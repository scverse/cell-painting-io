from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from skimage.segmentation import find_boundaries

from cell_painting_io.spatial import fov_offsets, labels_from_outlines, parse_well


@pytest.mark.parametrize(("well", "expected"), [("A01", (0, 0)), ("B02", (1, 1)), ("P24", (15, 23)), ("AA01", (26, 0))])
def test_parse_well(well: str, expected: tuple[int, int]) -> None:
    assert parse_well(well) == expected


def test_parse_well_rejects_nonsense() -> None:
    with pytest.raises(ValueError, match="not a well name"):
        parse_well("well-1")


def _positions() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "well": ["A01", "A01", "A01", "A01", "B02"],
            "x": [-5e-4, 5e-4, -5e-4, 5e-4, -5e-4],
            "y": [5e-4, 5e-4, -5e-4, -5e-4, 5e-4],
        }
    )


def test_fov_offsets_flip_y_and_share_a_well_origin() -> None:
    offsets = fov_offsets(_positions(), pixel_size=1e-6, plate_format=384)
    expected = pd.DataFrame({"well_y": [0.0, 0.0, 1000.0, 1000.0, 0.0], "well_x": [0.0, 1000.0, 0.0, 1000.0, 0.0]})
    pd.testing.assert_frame_equal(offsets[["well_y", "well_x"]], expected)


def test_fov_offsets_place_wells_on_the_plate_grid() -> None:
    offsets = fov_offsets(_positions(), pixel_size=1e-6, plate_format=384)
    expected = pd.DataFrame(
        {"plate_y": [0.0, 0.0, 1000.0, 1000.0, 4500.0], "plate_x": [0.0, 1000.0, 0.0, 1000.0, 4500.0]}
    )
    pd.testing.assert_frame_equal(offsets[["plate_y", "plate_x"]], expected)


def test_fov_offsets_without_a_plate_format() -> None:
    assert "plate_y" not in fov_offsets(_positions(), pixel_size=1e-6, plate_format=None)


def test_fov_offsets_reject_wells_off_the_plate() -> None:
    positions = _positions().assign(well="Z99")
    with pytest.raises(ValueError, match="outside a 96-well plate"):
        fov_offsets(positions, pixel_size=1e-6, plate_format=96)


def _outlined(labels: np.ndarray) -> np.ndarray:
    return find_boundaries(labels, mode="inner")


def test_labels_from_outlines_recovers_touching_objects() -> None:
    truth = np.zeros((40, 60), np.uint32)
    truth[5:35, 5:30] = 7
    truth[5:35, 30:55] = 3
    outlines = _outlined(truth)
    centres = pd.DataFrame(
        {"ObjectNumber": [7, 3], "Location_Center_X": [17.0, 42.0], "Location_Center_Y": [20.0, 20.0]}
    )

    labels = labels_from_outlines(outlines, centres)

    assert set(np.unique(labels)) == {0, 3, 7}
    for number in (3, 7):
        assert ((labels == number) & (truth == number)).sum() / (truth == number).sum() > 0.98


def test_labels_from_outlines_drops_unmatched_components() -> None:
    truth = np.zeros((40, 60), np.uint32)
    truth[5:35, 5:30] = 7
    truth[5:35, 30:55] = 3
    centres = pd.DataFrame({"ObjectNumber": [7], "Location_Center_X": [17.0], "Location_Center_Y": [20.0]})

    labels = labels_from_outlines(_outlined(truth), centres)

    assert set(np.unique(labels)) == {0, 7}


def test_labels_from_outlines_rejects_a_stack() -> None:
    centres = pd.DataFrame({"ObjectNumber": [1], "Location_Center_X": [1.0], "Location_Center_Y": [1.0]})
    with pytest.raises(ValueError, match="2D outline image"):
        labels_from_outlines(np.zeros((2, 10, 10)), centres)


def _centres(numbers: list[int], x: list[float], y: list[float], area: list[float]) -> pd.DataFrame:
    return pd.DataFrame(
        {"ObjectNumber": numbers, "Location_Center_X": x, "Location_Center_Y": y, "AreaShape_Area": area}
    )


def test_labels_from_outlines_rejects_the_background() -> None:
    truth = np.zeros((40, 60), np.uint32)
    truth[5:35, 5:30] = 7
    outlines = _outlined(truth)
    # object 3's outline never closed, so its centroid sits in the background component
    centres = _centres([7, 3], [17.0, 45.0], [20.0, 20.0], [750.0, 600.0])

    labels = labels_from_outlines(outlines, centres)

    assert set(np.unique(labels)) == {0, 7}


def test_labels_from_outlines_rejects_a_component_two_objects_claim() -> None:
    truth = np.zeros((40, 60), np.uint32)
    truth[5:35, 5:55] = 1
    outlines = _outlined(truth)
    # the boundary between the two objects is missing, so one component holds both centroids
    centres = _centres([1, 2], [17.0, 42.0], [20.0, 20.0], [750.0, 750.0])

    labels = labels_from_outlines(outlines, centres)

    assert not labels.any()


def test_labels_from_outlines_area_check_can_be_turned_off() -> None:
    truth = np.zeros((40, 60), np.uint32)
    truth[5:35, 5:30] = 7
    centres = _centres([7], [17.0], [20.0], [1.0])

    assert not labels_from_outlines(_outlined(truth), centres).any()
    assert labels_from_outlines(_outlined(truth), centres, area_column=None).any()


def test_labels_from_outlines_with_no_objects() -> None:
    centres = _centres([], [], [], [])
    labels = labels_from_outlines(np.zeros((8, 8), np.uint8), centres)
    assert labels.shape == (8, 8)
    assert not labels.any()
