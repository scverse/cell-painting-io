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
CHANNELS = ("DNA", "RNA")
SHAPE = (16, 16)
PIXEL_SIZE = 1e-6


def _truth() -> np.ndarray:
    labels = np.zeros(SHAPE, np.uint32)
    labels[2:8, 2:13] = 1
    labels[9:14, 2:13] = 2
    return labels


def _write_site(directory: Path) -> None:
    truth = _truth()
    (directory / "outlines").mkdir(parents=True)
    well, site = directory.name.split("-")[1:]
    for name, mask in (("cell", truth), ("nuclei", truth)):
        outlines = (find_boundaries(mask, mode="inner") * 255).astype(np.uint8)
        iio.imwrite(directory / "outlines" / f"{well}_s{site}--{name}_outlines.png", outlines)
    objects = pd.DataFrame(
        {
            "ImageNumber": [1, 1],
            "ObjectNumber": [1, 2],
            "Location_Center_X": [7.0, 7.0],
            "Location_Center_Y": [4.0, 11.0],
            "AreaShape_Area": [66.0, 55.0],
            "Intensity_MeanIntensity_DNA": [0.5, 0.25],
        }
    )
    for name in ("Cells", "Nuclei"):
        objects.to_csv(directory / f"{name}.csv", index=False)


def _write_plate(root: Path, plate: str) -> None:
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
        _write_site(analysis / f"{plate}-A01-{site}")

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
    """A minimal Cell Painting Gallery source tree: two plates, two wells each, one well of them analysed."""
    for plate in PLATES:
        _write_plate(tmp_path, plate)
    return tmp_path
