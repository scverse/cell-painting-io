from __future__ import annotations

from pathlib import Path

import anndata as ad
import dask.array as da
import h5py
import numpy as np
import numpy.typing as npt
import pandas as pd
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, Labels2DModel, TableModel
from spatialdata.transformations import Identity

# Names the ExportForSpatialData CellProfiler module writes. Kept here rather than imported, because
# the exporter is a CellProfiler plugin and is not installable alongside this package.
DATASET = "data"
MANIFEST = "cellprofiler_mapping"
ELEMENTS = "elements"
IMAGE_CHANNELS = "image_channels"
ELEMENT_IMAGE = "image"
ELEMENT_LABELS = "labels"
STATUS_OK = "ok"
REGION_KEY = "region_key"
INSTANCE_KEY = "label_id"

# Labels are cast to one type on the way in. CellProfiler narrows its own label arrays to fit the
# object count, so a field with 102 objects arrives as int8 and a busier one as int16, which would
# otherwise make two fields of the same plate disagree on dtype.
LABEL_DTYPE = np.uint32


def _table_path(plate: Path) -> Path:
    tables = sorted((plate / "tables").glob("*.h5ad"))
    if not tables:
        raise FileNotFoundError(f"no table under {plate / 'tables'}; is {plate} a plate folder of an export?")
    if len(tables) > 1:
        raise ValueError(f"expected one table under {plate / 'tables'}, found {[p.name for p in tables]}")
    return tables[0]


def _manifest(adata: ad.AnnData, plate: Path) -> tuple[pd.DataFrame, list[str]]:
    mapping = adata.uns.get(MANIFEST, {})
    if ELEMENTS not in mapping:
        raise ValueError(
            f"{plate} has no uns['{MANIFEST}']['{ELEMENTS}'], so it holds no images or segmentations. "
            "ExportToAnnData writes a table alone; ExportForSpatialData writes the manifest this reader needs."
        )
    channels = mapping[IMAGE_CHANNELS].sort_values("stack_index")["channel"].astype(str).tolist()
    return mapping[ELEMENTS], channels


def _read_array(path: Path, *, lazy: bool) -> npt.NDArray | da.Array:
    if not lazy:
        with h5py.File(path, "r") as handle:
            return handle[DATASET][()]
    # The file stays open for as long as the dask array holds the dataset. `lock` serialises the
    # reads, which HDF5 needs when dask runs them on several threads.
    dataset = h5py.File(path, "r")[DATASET]
    return da.from_array(dataset, chunks=dataset.chunks or "auto", lock=True)


def read_cellprofiler_export(path: Path | str, *, lazy: bool = True) -> SpatialData:
    """Read one plate folder written by the ExportForSpatialData CellProfiler module.

    The module writes one folder per plate, holding an image stack and a label array per field of view and one
    table for the plate, with a manifest in the table's `uns` listing every array it wrote.
    This reader builds only what the manifest names and resolves each path relative to `path`, so it never walks
    the folder or reads a file name, and a folder that was moved or renamed still reads.

    Each field of view becomes one Image, named `{field}_image`, with its channels in the order the module stacked
    them. Each segmented object becomes one Labels element per field, named as the table's `region_key` column
    names it, so the rows join onto it. Labels arrive as `uint32` whatever width CellProfiler used.

    Every element of a field sits in one coordinate system named after that field, and nothing places the fields
    relative to each other: the module does not export stage coordinates yet, so a plate reads as unstitched
    fields rather than a well or plate mosaic.

    Element names come from the exporter as they are, which does not yet match the `{plate}_{well}_s{site}_cells`
    convention `read_plate` uses. The two are to be reconciled, and until they are, an object read here and one
    read from the gallery carry different element names for the same thing.

    Args:
        path: A plate folder of an export, the directory holding `images/`, `labels/` and `tables/`.
        lazy: Read arrays through dask, one HDF5 dataset per element, instead of loading them into memory.

    Returns:
        The plate, with the Images and Labels the manifest lists as written and a `cells` Table annotating them.
        Arrays the module recorded as failed are left out, and so are the table rows that annotate them.

    Raises:
        FileNotFoundError: No table under `path/tables`, or the manifest names an array that is not there.
        ValueError: Several tables under `path/tables`, or the table carries no element manifest, or it has no
            `region_key` column to join its rows onto the label arrays.
    """
    plate = Path(path)
    adata = ad.read_h5ad(_table_path(plate))
    elements, channels = _manifest(adata, plate)
    written = elements[elements["status"] == STATUS_OK]

    images: dict[str, object] = {}
    labels: dict[str, object] = {}
    for row in written.itertuples():
        transformations = {row.sample_key: Identity()}
        array = _read_array(plate / row.path, lazy=lazy)
        if row.element_type == ELEMENT_IMAGE:
            images[f"{row.sample_key}_image"] = Image2DModel.parse(
                array, dims=("c", "y", "x"), c_coords=channels, transformations=transformations
            )
        elif row.element_type == ELEMENT_LABELS:
            labels[row.region_key_value] = Labels2DModel.parse(
                array.astype(LABEL_DTYPE), dims=("y", "x"), transformations=transformations
            )

    tables = {}
    if labels:
        if REGION_KEY not in adata.obs:
            raise ValueError(
                f"the table in {plate} has no obs['{REGION_KEY}'] column, so its rows cannot be joined onto the "
                "label arrays. ExportForSpatialData adds it; a table from ExportToAnnData does not carry it."
            )
        regions = sorted(set(adata.obs[REGION_KEY].astype(str)) & set(labels))
        table = adata[adata.obs[REGION_KEY].astype(str).isin(regions)].copy()
        table.obs[REGION_KEY] = pd.Categorical(table.obs[REGION_KEY].astype(str), categories=regions)
        # The exporter names every region it wrote, including the ones whose arrays failed, which this reader
        # leaves out. The region list therefore has to be rebuilt from the elements that exist, and
        # TableModel.parse refuses to run while the key is still set.
        table.uns.pop(TableModel.ATTRS_KEY, None)
        tables["cells"] = TableModel.parse(table, region=regions, region_key=REGION_KEY, instance_key=INSTANCE_KEY)
    return SpatialData(images=images, labels=labels, tables=tables)


def cellprofiler_export_plates(root: Path | str) -> list[Path]:
    """List the plate folders of one export root, for reading a run that covered several plates.

    Args:
        root: The `<prefix>_export` directory the module wrote.

    Returns:
        The plate folders, sorted by name. Read each with `read_cellprofiler_export` and combine them with
        `spatialdata.concatenate`; element names carry the field-of-view key, which is unique across plates when
        the pipeline defines a plate metadata tag.
    """
    return sorted(p for p in Path(root).iterdir() if (p / "tables").is_dir())
