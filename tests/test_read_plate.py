from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import spatialdata as sd
from conftest import OVERLAY_PLATE

from cell_painting_io import read_plate

BATCH = "2020_01_01_TEST"
PLATE = "BR00000001"
OTHER = "BR00000002"


@pytest.fixture
def sdata(gallery: Path) -> sd.SpatialData:
    return read_plate(gallery, BATCH, PLATE, profile="test")


def test_reads_every_downloaded_well(sdata: sd.SpatialData) -> None:
    assert set(sdata.images) == {f"{PLATE}_{w}_s{s}_image" for w in ("A01", "B02") for s in (1, 2)}


def test_labels_only_where_cellprofiler_ran(sdata: sd.SpatialData) -> None:
    assert set(sdata.labels) == {f"{PLATE}_A01_s{s}_{o}" for s in (1, 2) for o in ("nuclei", "cells", "cytoplasm")}


def test_channels_come_from_load_data(sdata: sd.SpatialData) -> None:
    image = sdata[f"{PLATE}_A01_s1_image"]["scale0"]["image"]
    assert list(image.coords["c"].to_numpy()) == ["DNA", "RNA"]


def test_labels_carry_cellprofiler_object_numbers(sdata: sd.SpatialData) -> None:
    mask = sdata[f"{PLATE}_A01_s1_cells"].to_numpy()
    assert set(np.unique(mask)) == {0, 1, 2}


def test_cytoplasm_is_cells_minus_nuclei(sdata: sd.SpatialData) -> None:
    cells = sdata[f"{PLATE}_A01_s1_cells"].to_numpy()
    nuclei = sdata[f"{PLATE}_A01_s1_nuclei"].to_numpy()
    cytoplasm = sdata[f"{PLATE}_A01_s1_cytoplasm"].to_numpy()

    assert np.array_equal(cytoplasm, np.where(nuclei > 0, 0, cells))
    assert cytoplasm.any()


def test_tables_annotate_their_elements(sdata: sd.SpatialData) -> None:
    cells, wells = sdata.tables["cells"], sdata.tables["wells"]
    assert cells.n_obs == 4
    assert set(cells.obs["region"]) == {f"{PLATE}_A01_s{s}_cells" for s in (1, 2)}
    assert wells.n_obs == 384
    assert set(wells.obs["region"]) == {f"{PLATE}_wells"}
    assert len(sdata.shapes[f"{PLATE}_wells"]) == 384


def test_three_coordinate_systems_per_element(sdata: sd.SpatialData) -> None:
    assert set(sdata.coordinate_systems) >= {PLATE, f"{PLATE}_A01", f"{PLATE}_A01_s1"}
    extent = sd.get_extent(sdata[f"{PLATE}_B02_s1_image"], coordinate_system=PLATE)
    assert extent["y"][0] == pytest.approx(4500.0)
    assert extent["x"][0] == pytest.approx(4500.0)


def test_selecting_wells(gallery: Path) -> None:
    sdata = read_plate(gallery, BATCH, PLATE, wells=["B02"], profile="test")
    assert set(sdata.images) == {f"{PLATE}_B02_s{s}_image" for s in (1, 2)}
    assert not sdata.labels
    assert set(sdata.tables) == {"wells"}


def test_without_a_profile_there_is_no_well_table(gallery: Path) -> None:
    sdata = read_plate(gallery, BATCH, PLATE, profile=None)
    assert set(sdata.tables) == {"cells"}
    assert not sdata.shapes


def test_two_plates_concatenate(sdata: sd.SpatialData, gallery: Path) -> None:
    other = read_plate(gallery, BATCH, OTHER, profile="test")
    assert not set(sdata.images) & set(other.images)

    merged = sd.concatenate([sdata, other], concatenate_tables=True)

    assert len(merged.images) == 8
    assert merged.tables["wells"].n_obs == 768
    assert merged.tables["cells"].n_obs == 8
    assert set(merged.coordinate_systems) >= {PLATE, OTHER}


def test_outlines_drawn_in_colour_over_the_image(gallery: Path) -> None:
    """The second plate publishes an overlay, and its cell image carries the nuclei outlines in a second colour."""
    sdata = read_plate(gallery, BATCH, OVERLAY_PLATE, profile="test")

    for kind in ("cells", "nuclei"):
        mask = sdata[f"{OVERLAY_PLATE}_A01_s1_{kind}"].to_numpy()
        assert set(np.unique(mask)) == {0, 1, 2}
    cells = sdata[f"{OVERLAY_PLATE}_A01_s1_cells"].to_numpy()
    nuclei = sdata[f"{OVERLAY_PLATE}_A01_s1_nuclei"].to_numpy()
    assert (cells == 1).sum() > (nuclei == 1).sum()
    assert ((nuclei == 1) & (cells == 1)).sum() == (nuclei == 1).sum()


def test_channels_ignore_illumination_columns(sdata: sd.SpatialData) -> None:
    image = sdata[f"{PLATE}_A01_s1_image"]["scale0"]["image"]
    assert list(image.coords["c"].to_numpy()) == ["DNA", "RNA"]
