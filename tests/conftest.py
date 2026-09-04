from __future__ import annotations

import gzip
from pathlib import Path

import imageio.v3 as iio
import numpy as np
import pandas as pd
import pytest
from skimage.segmentation import find_boundaries

BATCH = "2020_01_01_TEST"
PLATES = ("BR00000001", "BR00000002")
PLATE = PLATES[0]
OVERLAY_PLATE = PLATES[1]
CHANNELS = ("DNA", "RNA")
SHAPE = (32, 32)
PIXEL_SIZE = 1e-6


def _truth(kind: str) -> np.ndarray:
    """Two objects; nuclei sit inside cells, as CellProfiler grows them."""
    labels = np.zeros(SHAPE, np.uint32)
    if kind == "Cells":
        labels[3:16, 3:29] = 1
        labels[17:29, 3:29] = 2
    else:
        labels[7:13, 10:22] = 1
        labels[20:26, 10:22] = 2
    return labels


def _outlines(labels: np.ndarray) -> np.ndarray:
    return find_boundaries(labels, mode="inner")


def _write_site(directory: Path, well: str, site: int, *, overlay: bool) -> None:
    (directory / "outlines").mkdir(parents=True)
    for kind, name in (("Cells", "cell"), ("Nuclei", "nuclei")):
        truth = _truth(kind)
        path = directory / "outlines" / f"{well}_s{site}--{name}_outlines"
        if not overlay:
            iio.imwrite(path.with_suffix(".png"), (_outlines(truth) * 255).astype(np.uint8))
            continue
        # the outlines drawn over the greyscale image, and the cell image carries the nuclei outlines too
        image = np.full((*SHAPE, 3), 40, np.uint8)
        image[_outlines(truth)] = (0, 255, 255) if kind == "Cells" else (0, 128, 0)
        if kind == "Cells":
            image[_outlines(_truth("Nuclei"))] = (0, 128, 0)
        iio.imwrite(path.with_suffix(".tiff"), image)

    for kind in ("Cells", "Nuclei"):
        truth = _truth(kind)
        pd.DataFrame(
            {
                "ImageNumber": [1, 1],
                "ObjectNumber": [1, 2],
                "AreaShape_Center_X": [16.0, 16.0],
                "AreaShape_Center_Y": [
                    float(np.argwhere(truth == 1)[:, 0].mean()),
                    float(np.argwhere(truth == 2)[:, 0].mean()),
                ],
                "AreaShape_Area": [float((truth == 1).sum()), float((truth == 2).sum())],
                "Intensity_MeanIntensity_DNA": [0.5, 0.25],
            }
        ).to_csv(directory / f"{kind}.csv", index=False)


def _write_plate(root: Path, plate: str, *, overlay: bool) -> None:
    folder = f"{plate}__2020-01-01T00_00_00-Measurement1"
    images = root / "images" / BATCH / "images" / folder / "Images"
    images.mkdir(parents=True)
    rows = []
    for well, (row, column) in (("A01", (1, 1)), ("B02", (2, 2))):
        for site, (x, y) in enumerate([(-8e-6, 8e-6), (8e-6, -8e-6)], start=1):
            urls = {}
            for index, channel in enumerate(CHANNELS, start=1):
                name = f"r{row:02d}c{column:02d}f{site:02d}p01-ch{index}.tiff"
                iio.imwrite(images / name, np.full(SHAPE, index * 100, np.uint16))
                urls[f"URL_Orig{channel}"] = f"s3://bucket/acc/source_1/images/{BATCH}/images/{folder}/Images/{name}"
            rows.append(
                {
                    **urls,
                    "URL_IllumDNA": "ignored",
                    "Metadata_Plate": plate,
                    "Metadata_Well": well,
                    "Metadata_Site": site,
                    "Metadata_PositionX": x,
                    "Metadata_PositionY": y,
                    "Metadata_ImageResolutionX": PIXEL_SIZE,
                    "Metadata_ImageResolutionY": PIXEL_SIZE,
                    "Metadata_ImageSizeX": SHAPE[1],
                    "Metadata_ImageSizeY": SHAPE[0],
                }
            )
    load_data = root / "workspace/load_data_csv" / BATCH / plate
    load_data.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(load_data / "load_data.csv", index=False)

    analysis = root / "workspace/analysis" / BATCH / plate / "analysis"
    for site in (1, 2):
        _write_site(analysis / f"{plate}-A01-{site}", "A01", site, overlay=overlay)

    profiles = root / "workspace/profiles" / BATCH / plate
    profiles.mkdir(parents=True)
    plate_map = pd.DataFrame(
        {
            "Metadata_Plate": [plate] * 384,
            "Metadata_Well": [f"{chr(ord('A') + r)}{c + 1:02d}" for r in range(16) for c in range(24)],
            "Cells_AreaShape_Area": np.arange(384, dtype=float),
        }
    )
    with gzip.open(profiles / f"{plate}_test.csv.gz", "wt") as fh:
        plate_map.to_csv(fh, index=False)


@pytest.fixture
def gallery(tmp_path: Path) -> Path:
    """A minimal gallery source: two plates, two wells each, one well analysed, the second plate using overlays."""
    for plate in PLATES:
        _write_plate(tmp_path, plate, overlay=plate == OVERLAY_PLATE)
    return tmp_path
