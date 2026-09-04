from __future__ import annotations

import re
from collections.abc import Mapping

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import ndimage as ndi
from skimage.segmentation import expand_labels

PLATE_FORMATS: Mapping[int, tuple[int, int, float]] = {
    96: (8, 12, 9.0e-3),
    384: (16, 24, 4.5e-3),
    1536: (32, 48, 2.25e-3),
}

_WELL = re.compile(r"([A-Za-z]+)(\d+)")


def parse_well(well: str) -> tuple[int, int]:
    """Split a well name into its zero-based row and column.

    Rows count the way spreadsheet columns do, so `A` is 0, `Z` is 25 and `AA` is 26, which is how 1536-well plates are named.

    Args:
        well: A well name such as `A01`, `P24` or `AF48`.

    Returns:
        The zero-based row and column.

    Raises:
        ValueError: The name is not letters followed by digits.
    """
    match = _WELL.fullmatch(well.strip())
    if match is None:
        raise ValueError(f"not a well name: {well!r}")
    letters, digits = match.groups()
    row = 0
    for letter in letters.upper():
        row = row * 26 + ord(letter) - ord("A") + 1
    return row - 1, int(digits) - 1


def fov_offsets(positions: pd.DataFrame, *, pixel_size: float, plate_format: int | None = 384) -> pd.DataFrame:
    """Pixel offsets of the top-left corner of each field of view, within its well and within the plate.

    The offsets have `y` pointing down, unlike the stage coordinates they come from, so that they can be used
    directly as `(y, x)` translations of an image whose row index grows downwards.
    Their origin is the top-left corner of the top-left field of a well, which is shared across wells, so every
    well of a plate is laid out in one frame.

    Args:
        positions: One row per field of view, with columns `well`, and `x` and `y` holding the stage coordinates
            of the field centre in metres, relative to the centre of its well and with `y` pointing up.
        pixel_size: Size of a pixel in metres, as `Metadata_ImageResolutionX` records it.
        plate_format: Number of wells on the plate, used to place the wells on their nominal grid.
            This assumes the standard well pitch of that format rather than anything measured; pass `None` to
            leave the plate offsets out.

    Returns:
        The offsets in pixels, indexed like `positions`, with columns `well_y` and `well_x`, and `plate_y` and
        `plate_x` unless `plate_format` is `None`.

    Raises:
        ValueError: The plate format is unknown, or a well falls outside a plate of that format.
    """
    x = positions["x"].to_numpy(float) / pixel_size
    y = positions["y"].to_numpy(float) / pixel_size
    offsets = pd.DataFrame({"well_y": y.max() - y, "well_x": x - x.min()}, index=positions.index)
    if plate_format is not None:
        if plate_format not in PLATE_FORMATS:
            raise ValueError(f"unknown plate format {plate_format}; known: {sorted(PLATE_FORMATS)}")
        n_rows, n_cols, pitch = PLATE_FORMATS[plate_format]
        grid = np.array([parse_well(well) for well in positions["well"]])
        if (grid < 0).any() or (grid >= [n_rows, n_cols]).any():
            raise ValueError(f"wells fall outside a {plate_format}-well plate")
        offsets[["plate_y", "plate_x"]] = offsets[["well_y", "well_x"]].to_numpy() + grid * (pitch / pixel_size)
    return offsets


def labels_from_outlines(
    outlines: npt.ArrayLike,
    centres: pd.DataFrame,
    *,
    x_column: str = "Location_Center_X",
    y_column: str = "Location_Center_Y",
    object_column: str = "ObjectNumber",
) -> npt.NDArray[np.uint32]:
    """Turn a CellProfiler outline image back into a label image carrying the CellProfiler object numbers.

    The Cell Painting Gallery publishes outlines rather than segmentation masks.
    Outlines are one pixel wide and shared between touching objects, so filling them and labelling connected
    components recovers the object interiors separately.
    Each component takes the object number of the centroid that falls inside it, and components without one -
    the background, and anything else the outlines close off - are dropped.
    The boundary is then grown back one pixel, which reproduces the CellProfiler areas to well under a percent.

    A component holding several centroids, which happens where an outline failed to close, keeps the lowest
    object number, so the objects that lost it are absent from the result.
    Compare the number of labels against the number of centroids to catch that.

    Args:
        outlines: A 2D outline image, anything above zero being outline.
        centres: One row per object, as CellProfiler wrote it, with the centroid and object number columns below.
        x_column: Column of `centres` holding the centroid column index.
        y_column: Column of `centres` holding the centroid row index.
        object_column: Column of `centres` holding the object number, which becomes the label value.

    Returns:
        A label image the shape of `outlines`, zero outside objects.

    Raises:
        ValueError: `outlines` is not 2D.
    """
    mask = np.asarray(outlines) > 0
    if mask.ndim != 2:
        raise ValueError(f"expected a 2D outline image, got shape {mask.shape}")
    components, n_components = ndi.label(ndi.binary_fill_holes(~mask) & ~mask)

    rows = np.rint(centres[y_column].to_numpy(float)).astype(int).clip(0, mask.shape[0] - 1)
    columns = np.rint(centres[x_column].to_numpy(float)).astype(int).clip(0, mask.shape[1] - 1)
    numbers = centres[object_column].to_numpy()
    order = np.argsort(numbers)[::-1]
    lookup = np.zeros(n_components + 1, dtype=np.uint32)
    lookup[components[rows, columns][order]] = numbers[order]
    lookup[0] = 0
    return expand_labels(lookup[components], distance=1).astype(np.uint32)
