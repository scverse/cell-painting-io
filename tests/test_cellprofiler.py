from __future__ import annotations

from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
import spatialdata as sd

from cell_painting_io.cellprofiler import cellprofiler_export_plates, read_cellprofiler_export

PREFIX = "Testrun"
PLATE = "BR00000001"
CHANNELS = ("DNA", "RNA")
OBJECTS = ("Nuclei", "Cells")
SHAPE = (8, 16)
FIELDS = ("A01_01_img1", "A01_02_img2")
ELEMENT_COLUMNS = (
    "sample_key",
    "image_number",
    "element_type",
    "element_name",
    "path",
    "shape",
    "element_dtype",
    "region_key_value",
    "status",
    "error",
)


def _write(path: Path, array: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as handle:
        handle.create_dataset("data", data=array, compression="gzip", compression_opts=1)


def _labels(n: int, dtype: str = "int8") -> np.ndarray:
    labels = np.zeros(SHAPE, dtype)
    for i in range(n):
        labels[2 * i + 1, 1 : 4 + i] = i + 1
    return labels


def _export(root: Path, *, fields: tuple[str, ...] = FIELDS, failed: str | None = None, label_dtype: str = "int8"):
    """One plate folder as ExportForSpatialData writes it, including its manifest.

    `failed` names a field whose arrays are not written, the way a cycle that raised leaves them, so its manifest
    rows carry a status of `failed` and an error.
    """
    plate = root / f"{PREFIX}_export" / PLATE
    rows, obs = [], []
    for number, field in enumerate(fields, start=1):
        broken = field == failed
        if not broken:
            _write(plate / "images" / f"{field}.h5", np.full((len(CHANNELS), *SHAPE), 0.5, "float32"))
        rows.append(
            {
                "sample_key": field,
                "image_number": number,
                "element_type": "image",
                "element_name": "",
                "path": f"images/{field}.h5",
                "shape": "" if broken else f"{len(CHANNELS)},{SHAPE[0]},{SHAPE[1]}",
                "element_dtype": "" if broken else "float32",
                "region_key_value": "",
                "status": "failed" if broken else "ok",
                "error": "RuntimeError: the cycle failed" if broken else "",
            }
        )
        for obj in OBJECTS:
            if not broken:
                _write(plate / "labels" / field / f"{obj}.h5", _labels(2, label_dtype))
            rows.append(
                {
                    "sample_key": field,
                    "image_number": number,
                    "element_type": "labels",
                    "element_name": obj,
                    "path": f"labels/{field}/{obj}.h5",
                    "shape": "" if broken else f"{SHAPE[0]},{SHAPE[1]}",
                    "element_dtype": "" if broken else label_dtype,
                    "region_key_value": f"{field}__{obj}",
                    "status": "failed" if broken else "ok",
                    "error": "RuntimeError: the cycle failed" if broken else "",
                }
            )
        obs += [{"region_key": f"{field}__Cells", "label_id": i + 1, "ImageNumber": number} for i in range(2)]

    frame = pd.DataFrame(obs)
    adata = ad.AnnData(
        np.arange(len(obs) * 3, dtype="float32").reshape(len(obs), 3),
        obs=frame.astype({"label_id": "int32", "ImageNumber": "int32"}),
        var=pd.DataFrame(index=["Cells__AreaShape_Area", "Cells__Intensity_MeanIntensity_DNA", "Nuclei__A"]),
    )
    adata.obs_names = [f"{r['region_key'].rsplit('__', 1)[0]}_{r['label_id']}" for r in obs]
    adata.uns["cellprofiler_mapping"] = {
        "elements": pd.DataFrame(rows, columns=list(ELEMENT_COLUMNS)),
        "image_channels": pd.DataFrame({"channel": list(CHANNELS), "stack_index": range(len(CHANNELS))}),
    }
    adata.uns["spatialdata_attrs"] = {
        "region": sorted(set(adata.obs["region_key"])),
        "region_key": "region_key",
        "instance_key": "label_id",
    }
    (plate / "tables").mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(plate / "tables" / f"{PREFIX}.h5ad")
    return plate


def test_reads_images_labels_and_the_table(tmp_path: Path) -> None:
    sdata = read_cellprofiler_export(_export(tmp_path))
    assert set(sdata.images) == {f"{field}_image" for field in FIELDS}
    assert set(sdata.labels) == {f"{field}__{obj}" for field in FIELDS for obj in OBJECTS}
    assert set(sdata.coordinate_systems) == set(FIELDS)
    assert sdata.tables["cells"].shape == (4, 3)


def test_channels_keep_the_order_the_module_stacked_them_in(tmp_path: Path) -> None:
    sdata = read_cellprofiler_export(_export(tmp_path))
    image = sdata.images[f"{FIELDS[0]}_image"]
    assert list(image.coords["c"].values) == list(CHANNELS)
    assert image.shape == (len(CHANNELS), *SHAPE)


def test_every_element_of_a_field_shares_that_field_s_coordinate_system(tmp_path: Path) -> None:
    """An image and the labels from the same field are the same pixel grid, which is what makes overlay work."""
    sdata = read_cellprofiler_export(_export(tmp_path))
    field = FIELDS[0]
    named = {name for _, name, _ in sdata.filter_by_coordinate_system(field).gen_spatial_elements()}
    assert named == {f"{field}_image", f"{field}__Nuclei", f"{field}__Cells"}


def test_rows_join_onto_the_label_arrays(tmp_path: Path) -> None:
    sdata = read_cellprofiler_export(_export(tmp_path))
    table = sdata.tables["cells"]
    for name in (f"{field}__Cells" for field in FIELDS):
        in_array = set(np.unique(np.asarray(sdata.labels[name])).tolist()) - {0}
        in_table = set(table.obs.loc[table.obs["region_key"] == name, "label_id"].tolist())
        assert in_array == in_table, name


def test_the_table_annotates_only_the_label_elements_that_exist(tmp_path: Path) -> None:
    """Nuclei elements are read but no row annotates them, so `region` must not claim them."""
    sdata = read_cellprofiler_export(_export(tmp_path))
    regions = sdata.tables["cells"].uns["spatialdata_attrs"]["region"]
    assert set(regions) == {f"{field}__Cells" for field in FIELDS}
    assert set(regions) <= set(sdata.labels)


def test_a_failed_field_is_left_out_along_with_its_rows(tmp_path: Path) -> None:
    """A cycle that failed wrote no arrays. The reader has to skip them and drop the rows that annotate them,
    rather than crash opening a file the manifest names."""
    sdata = read_cellprofiler_export(_export(tmp_path, failed=FIELDS[1]))
    assert set(sdata.images) == {f"{FIELDS[0]}_image"}
    assert set(sdata.labels) == {f"{FIELDS[0]}__{obj}" for obj in OBJECTS}
    table = sdata.tables["cells"]
    assert table.n_obs == 2
    assert set(table.obs["region_key"]) == {f"{FIELDS[0]}__Cells"}
    assert set(table.uns["spatialdata_attrs"]["region"]) <= set(sdata.labels)


@pytest.mark.parametrize("label_dtype", ["int8", "int16", "int32"])
def test_labels_arrive_as_one_type_whatever_width_cellprofiler_used(tmp_path: Path, label_dtype: str) -> None:
    """CellProfiler narrows its label arrays to fit the object count, so two fields of one plate can disagree."""
    sdata = read_cellprofiler_export(_export(tmp_path, label_dtype=label_dtype))
    assert sdata.labels[f"{FIELDS[0]}__Cells"].dtype == np.uint32


@pytest.mark.parametrize("lazy", [True, False])
def test_lazy_and_eager_reads_give_the_same_pixels(tmp_path: Path, lazy: bool) -> None:
    plate = _export(tmp_path)
    sdata = read_cellprofiler_export(plate, lazy=lazy)
    assert np.asarray(sdata.images[f"{FIELDS[0]}_image"]).shape == (len(CHANNELS), *SHAPE)
    assert np.unique(np.asarray(sdata.labels[f"{FIELDS[0]}__Cells"])).tolist() == [0, 1, 2]


def test_the_object_writes_to_zarr(tmp_path: Path) -> None:
    """The manifest travels in the table's `uns`, so anything unwritable in it makes the whole object
    unwritable. A manifest column named `dtype` did exactly that: pandas resolves `frame.dtype` to the column
    and anndata dispatches on `elem.dtype.kind`."""
    sdata = read_cellprofiler_export(_export(tmp_path))
    sdata.write(tmp_path / "plate.zarr")
    back = sd.read_zarr(tmp_path / "plate.zarr")
    assert set(back.labels) == set(sdata.labels)
    assert np.array_equal(np.asarray(back.tables["cells"].X), np.asarray(sdata.tables["cells"].X))
    assert "cellprofiler_mapping" in back.tables["cells"].uns


def test_a_table_without_a_manifest_says_which_module_writes_one(tmp_path: Path) -> None:
    plate = _export(tmp_path)
    path = plate / "tables" / f"{PREFIX}.h5ad"
    adata = ad.read_h5ad(path)
    del adata.uns["cellprofiler_mapping"]["elements"]
    adata.write_h5ad(path)
    with pytest.raises(ValueError, match="ExportForSpatialData writes the manifest"):
        read_cellprofiler_export(plate)


def test_a_table_without_region_key_is_refused(tmp_path: Path) -> None:
    plate = _export(tmp_path)
    path = plate / "tables" / f"{PREFIX}.h5ad"
    adata = ad.read_h5ad(path)
    del adata.obs["region_key"]
    del adata.uns["spatialdata_attrs"]
    adata.write_h5ad(path)
    with pytest.raises(ValueError, match="region_key"):
        read_cellprofiler_export(plate)


def test_a_folder_that_is_not_an_export_says_so(tmp_path: Path) -> None:
    (tmp_path / "tables").mkdir()
    with pytest.raises(FileNotFoundError, match="no table under"):
        read_cellprofiler_export(tmp_path)


def test_two_tables_in_one_plate_folder_are_refused(tmp_path: Path) -> None:
    plate = _export(tmp_path)
    (plate / "tables" / "second.h5ad").write_bytes(b"")
    with pytest.raises(ValueError, match="expected one table"):
        read_cellprofiler_export(plate)


def test_export_plates_lists_the_plate_folders(tmp_path: Path) -> None:
    plate = _export(tmp_path)
    other = plate.parent / "BR00000002"
    (other / "tables").mkdir(parents=True)
    (plate.parent / "notes.txt").write_text("not a plate")
    assert cellprofiler_export_plates(plate.parent) == [plate, other]
