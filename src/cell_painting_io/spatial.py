from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path

import anndata as ad
import imageio.v3 as iio
import numpy as np
import numpy.typing as npt
import pandas as pd
from geopandas import GeoDataFrame
from scipy import ndimage as ndi
from shapely import Point
from skimage.segmentation import expand_labels
from spatialdata import SpatialData, sanitize_table
from spatialdata.models import Image2DModel, Labels2DModel, ShapesModel, TableModel
from spatialdata.transformations import Identity, Translation

from cell_painting_io.profiles import annotate_features
from cell_painting_io.reader import read_profiles

PLATE_FORMATS: Mapping[int, tuple[int, int, float]] = {
    96: (8, 12, 9.0e-3),
    384: (16, 24, 4.5e-3),
    1536: (32, 48, 2.25e-3),
}

WELL_RADIUS: float = 1.65e-3

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
    area_column: str | None = "AreaShape_Area",
    max_area_ratio: float = 1.25,
) -> npt.NDArray[np.uint32]:
    """Turn a CellProfiler outline image back into a label image carrying the CellProfiler object numbers.

    The Cell Painting Gallery publishes outlines rather than segmentation masks.
    Outlines are one pixel wide and shared between touching objects, so filling them and labelling connected
    components recovers the object interiors separately.
    Each component takes the object number of the centroid that falls inside it, and components without one -
    the background, and anything else the outlines close off - are dropped.
    The boundary is then grown back one pixel, which reproduces the CellProfiler areas to well under a percent.

    A component is only taken when it unambiguously belongs to one object: exactly one centroid falls in it, and
    its area is within `max_area_ratio` of the area CellProfiler measured.
    That rejects the two ways an unclosed outline goes wrong - two objects merging into one component, and a
    centroid landing in the background or in a fragment of its object - so every label in the result stands for
    one CellProfiler object rather than for a guess.
    Objects whose component was rejected are absent; compare the number of labels against the number of centroids.

    Args:
        outlines: A 2D outline image, anything above zero being outline.
        centres: One row per object, as CellProfiler wrote it, with the centroid and object number columns below.
        x_column: Column of `centres` holding the centroid column index.
        y_column: Column of `centres` holding the centroid row index.
        object_column: Column of `centres` holding the object number, which becomes the label value.
        area_column: Column of `centres` holding the area CellProfiler measured, used to reject components that
            cannot be the object. Pass `None`, or leave the column out of `centres`, to skip that check.
        max_area_ratio: How far a component's area may differ from the measured area, either way, and still be
            accepted. Reconstruction is exact to a fraction of a percent when it works, so the default is tight.

    Returns:
        A label image the shape of `outlines`, zero outside objects.

    Raises:
        ValueError: `outlines` is not 2D.
    """
    mask = np.asarray(outlines) > 0
    if mask.ndim != 2:
        raise ValueError(f"expected a 2D outline image, got shape {mask.shape}")
    if centres.empty:
        return np.zeros(mask.shape, np.uint32)
    components, n_components = ndi.label(ndi.binary_fill_holes(~mask) & ~mask)

    rows = np.rint(centres[y_column].to_numpy(float)).astype(int).clip(0, mask.shape[0] - 1)
    columns = np.rint(centres[x_column].to_numpy(float)).astype(int).clip(0, mask.shape[1] - 1)
    numbers = centres[object_column].to_numpy()
    areas = centres[area_column].to_numpy(float) if area_column is not None and area_column in centres else None
    hit = components[rows, columns]
    if areas is not None:
        # an interior is always smaller than the object, so only the upper bound means anything here; it is what
        # catches a centroid that landed in the background
        sizes = np.bincount(components.ravel(), minlength=n_components + 1)
        hit = np.where(sizes[hit] > max_area_ratio * areas, 0, hit)
    claims = np.bincount(hit, minlength=n_components + 1)
    hit = np.where(claims[hit] > 1, 0, hit)

    lookup = np.zeros(n_components + 1, dtype=np.uint32)
    lookup[hit] = numbers
    lookup[0] = 0
    labels = expand_labels(lookup[components], distance=1).astype(np.uint32)

    if areas is not None:
        grown = np.bincount(labels.ravel(), minlength=int(numbers.max()) + 1)[numbers] / areas
        rejected = numbers[(grown > max_area_ratio) | (grown < 1 / max_area_ratio)]
        labels[np.isin(labels, rejected)] = 0
    return labels


def _load_data(root: Path, batch: str, plate: str) -> pd.DataFrame:
    frame = pd.read_csv(root / "workspace/load_data_csv" / batch / plate / "load_data.csv")
    if "Metadata_Well" not in frame.columns:
        raise ValueError("load_data.csv does not name the well of each field")
    return frame.set_index(["Metadata_Well", "Metadata_Site"]).rename_axis(["well", "site"]).sort_index()


def _pixel_size(load_data: pd.DataFrame) -> float:
    sizes = np.unique(load_data[["Metadata_ImageResolutionX", "Metadata_ImageResolutionY"]].to_numpy().round(12))
    if sizes.size != 1:
        raise ValueError(f"the plate mixes pixel sizes: {sizes}")
    return float(sizes.item())


def _channels(load_data: pd.DataFrame) -> list[str]:
    for prefix in ("URL_Orig", "FileName_Orig"):
        names = [c.removeprefix(prefix) for c in load_data.columns if c.startswith(prefix)]
        if names:
            return names
    raise ValueError("load_data.csv has neither URL_Orig nor FileName_Orig columns")


def _image_path(root: Path, batch: str, row: pd.Series, channel: str) -> Path:
    """Locate a channel of one field under `root`.

    Sources record images either as an `URL_Orig*` S3 URI or as a `FileName_Orig*`/`PathName_Orig*` pair whose
    path is wherever the images sat when CellProfiler ran.
    Both end in the gallery's own `<batch>/images/<acquisition>/Images/` layout, which is what is matched on.
    """
    if f"URL_Orig{channel}" in row:
        location = str(row[f"URL_Orig{channel}"])
    else:
        location = f"{row[f'PathName_Orig{channel}']}/{row[f'FileName_Orig{channel}']}"
    marker = f"/{batch}/images/"
    if marker not in location:
        raise ValueError(f"image location does not follow the gallery layout: {location}")
    return root / "images" / location[location.index(marker) + 1 :]


def _read_fov(root: Path, batch: str, row: pd.Series, channels: Sequence[str]) -> npt.NDArray:
    return np.stack([iio.imread(_image_path(root, batch, row, channel)) for channel in channels])


def _site_dir(root: Path, batch: str, plate: str, well: str, site: int) -> Path:
    return root / "workspace/analysis" / batch / plate / "analysis" / f"{plate}-{well}-{site}"


def _outline_file(directory: Path, well: str, site: int, kind: str) -> Path | None:
    """Find one outline image, whatever the source called it.

    Seen across the gallery: `outlines/A01_s1--cell_outlines.png`, `outlines/a01_1--cell_outlines.png` and, in a
    directory named after the plate, `A01_s1_cell_outlines.tiff`.
    """
    for stem in (
        f"{well}_s{site}--{kind}_outlines",
        f"{well}_{site}--{kind}_outlines",
        f"{well}_s{site}_{kind}_outlines",
    ):
        matches = sorted(directory.glob(f"{stem}.*")) + sorted(directory.glob(f"*/{stem}.*"))
        if matches:
            return matches[0]
    return None


def _centre_columns(objects: pd.DataFrame) -> tuple[str, str]:
    for prefix in ("Location_Center", "AreaShape_Center"):
        if {f"{prefix}_X", f"{prefix}_Y"}.issubset(objects.columns):
            return f"{prefix}_X", f"{prefix}_Y"
    raise ValueError("object table has neither Location_Center nor AreaShape_Center columns")


def _site_labels(directory: Path, well: str, site: int) -> dict[str, npt.NDArray[np.uint32]]:
    masks = {}
    for name, kind, csv in (("nuclei", "nuclei", "Nuclei"), ("cells", "cell", "Cells")):
        path = _outline_file(directory, well, site, kind)
        if path is None or not (directory / f"{csv}.csv").exists():
            return {}
        outlines = np.squeeze(iio.imread(path))
        if outlines.ndim != 2:
            # some sources publish a colour overlay of the outlines on the image rather than the outlines alone
            return {}
        objects = pd.read_csv(directory / f"{csv}.csv")
        x_column, y_column = _centre_columns(objects)
        masks[name] = labels_from_outlines(outlines, objects, x_column=x_column, y_column=y_column)
    masks["cytoplasm"] = np.where(masks["nuclei"] > 0, 0, masks["cells"])
    return masks


def _well_shapes(load_data: pd.DataFrame, *, pixel_size: float, plate_format: int, system: str) -> GeoDataFrame:
    x = load_data["Metadata_PositionX"].to_numpy(float) / pixel_size
    y = load_data["Metadata_PositionY"].to_numpy(float) / pixel_size
    height, width = int(load_data["Metadata_ImageSizeY"].iloc[0]), int(load_data["Metadata_ImageSizeX"].iloc[0])
    centre = (y.max() + height / 2, -x.min() + width / 2)
    n_rows, n_cols, pitch = PLATE_FORMATS[plate_format]
    rows, columns = np.divmod(np.arange(n_rows * n_cols), n_cols)
    points = [
        Point(centre[1] + c * pitch / pixel_size, centre[0] + r * pitch / pixel_size)
        for r, c in zip(rows, columns, strict=True)
    ]
    frame = GeoDataFrame({"radius": WELL_RADIUS / pixel_size}, geometry=points, index=np.arange(n_rows * n_cols))
    return ShapesModel.parse(frame, transformations={system: Identity()})


def _well_table(path: Path, *, region: str | None, plate_format: int) -> ad.AnnData:
    adata = read_profiles(path, index_columns=("Plate", "Well"))
    annotate_features(adata.var)
    # some sources annotate wells with column names SpatialData rejects, "Common Name" among them
    sanitize_table(adata)
    if region is None:
        return TableModel.parse(adata)
    grid = np.array([parse_well(well) for well in adata.obs["Well"]])
    adata.obs["well_index"] = grid[:, 0] * PLATE_FORMATS[plate_format][1] + grid[:, 1]
    adata.obs["region"] = pd.Categorical([region] * adata.n_obs)
    return TableModel.parse(adata, region=region, region_key="region", instance_key="well_index")


def _cell_table(files: Sequence[Path], masks: Mapping[str, npt.NDArray]) -> ad.AnnData:
    adata = read_profiles(files, metadata_columns=("ImageNumber", "ObjectNumber"), path_columns={"Metadata_Key": 1})
    keys = adata.obs.pop("Key").str.rsplit("-", n=2, expand=True)
    plates, wells, sites = (keys[i].astype(str) for i in range(3))
    adata.obs["Well"] = pd.Categorical(wells)
    adata.obs["Site"] = sites.astype(int)
    adata.obs["region"] = pd.Categorical(plates + "_" + wells + "_s" + sites + "_cells")
    adata.obs_names = pd.Index(adata.obs["region"].astype(str) + ":" + adata.obs["ObjectNumber"].astype(str))
    annotate_features(adata.var)

    numbers = adata.obs["ObjectNumber"].to_numpy()
    keep = np.zeros(adata.n_obs, bool)
    for region, mask in masks.items():
        rows = (adata.obs["region"] == region).to_numpy()
        keep[rows] = np.isin(numbers[rows], np.unique(mask))
    adata = adata[keep].copy()
    adata.obs["region"] = adata.obs["region"].cat.remove_unused_categories()
    sanitize_table(adata)
    return TableModel.parse(
        adata, region=sorted(adata.obs["region"].cat.categories), region_key="region", instance_key="ObjectNumber"
    )


def read_plate(
    root: Path | str,
    batch: str,
    plate: str,
    *,
    wells: Sequence[str] | None = None,
    plane: int | None = None,
    profile: str | Path | None = "normalized_feature_select_negcon_batch",
    plate_format: int = 384,
) -> SpatialData:
    """Read one plate of a Cell Painting Gallery source into a SpatialData object.

    Fields of view become Images, one element per field with a channel per `URL_Orig*` or `FileName_Orig*` column
    of `load_data.csv`.
    The CellProfiler Nuclei, Cells and Cytoplasm segmentations become Labels, reconstructed from the published
    outlines by `labels_from_outlines` and carrying CellProfiler's own object numbers.
    The wells of the plate become Shapes, and the well- and cell-level measurements become Tables named
    `wells` and `cells`, each annotating the elements above.
    Cytoplasm is the cell mask minus the nucleus mask, as CellProfiler defines it, so it needs no files of its own.

    Where the source recorded stage coordinates, every element is placed in three coordinate systems, named
    `{plate}_{well}_s{site}`, `{plate}_{well}` and `{plate}`, so the fields of a well lay out as a mosaic and the
    wells as a plate map.
    Where it did not, or where it left out the pixel size those coordinates would be converted with, each field
    can only sit in its own frame, and the well shapes are left out along with the plate frame they would live in.
    Element names are prefixed with the plate barcode either way, so two plates concatenate without renaming.

    The gallery is uneven about what it publishes.
    A field contributes an image whether or not CellProfiler output exists for it, and labels only when that
    output includes outlines it can use; sources that publish a colour overlay instead of the outlines alone, or
    that pool a whole well into one analysis directory, yield images and tables but no labels.
    Rows of the cell table whose object did not survive the outline reconstruction are dropped, so that every row
    points at a label that exists.

    Args:
        root: Directory holding the `images/` and `workspace/` trees of one source of one accession.
        batch: Batch name, the directory below `images/` and `workspace/analysis/`.
        plate: Plate barcode.
        wells: Wells to read images and labels for. Defaults to every well whose images are present under `root`,
            so that a partial download reads back as itself. The well table always covers the whole plate.
        plane: Which `Metadata_PlaneID` to read where a source imaged a z stack. Required in that case, since
            there is no reason to prefer one plane over another and taking one silently would hide the rest.
        profile: Variant of the well-level profile, read as `workspace/profiles/{batch}/{plate}/{plate}_{profile}.csv.gz`.
            Sources that publish the profile under another name, `{plate}.parquet` among them, take a path instead.
            Pass `None` to leave out the well table and the well shapes.
        plate_format: Number of wells on the plate, used to place the wells on their nominal grid.

    Returns:
        The plate, with Images, Labels, Shapes and the `wells` and `cells` Tables.
        A table or element group is left out when nothing it would hold was read.

    Raises:
        ValueError: The plate mixes pixel sizes, `load_data.csv` names its images in an unknown way, or an image
            does not sit in the gallery layout.
        FileNotFoundError: `load_data.csv` or the requested profile is missing.
    """
    root = Path(root)
    load_data = _load_data(root, batch, plate)
    if load_data.index.duplicated().any():
        planes = sorted(load_data["Metadata_PlaneID"].unique()) if "Metadata_PlaneID" in load_data else []
        if plane is None:
            raise ValueError(f"{plate} has several rows per field; pass plane= to choose one of {planes}")
        load_data = load_data[load_data["Metadata_PlaneID"] == plane]
    channels = _channels(load_data)
    located = {
        "Metadata_PositionX",
        "Metadata_PositionY",
        "Metadata_ImageResolutionX",
        "Metadata_ImageResolutionY",
    }.issubset(load_data.columns)
    if located:
        pixel_size = _pixel_size(load_data)
        offsets = fov_offsets(
            pd.DataFrame(
                {
                    "well": load_data.index.get_level_values("well"),
                    "x": load_data["Metadata_PositionX"].to_numpy(float),
                    "y": load_data["Metadata_PositionY"].to_numpy(float),
                },
                index=load_data.index,
            ),
            pixel_size=pixel_size,
            plate_format=plate_format,
        )

    if wells is None:
        first = load_data.groupby(level="well").head(1)
        wells = [key[0] for key, row in first.iterrows() if _image_path(root, batch, row, channels[0]).exists()]

    images, labels, masks, analysed = {}, {}, {}, []
    for well in wells:
        for site in sorted(load_data.loc[well].index):
            fov = f"{plate}_{well}_s{site}"
            transformations: dict[str, Identity | Translation] = {fov: Identity()}
            if located:
                offset = offsets.loc[(well, site)]
                transformations[f"{plate}_{well}"] = Translation([offset["well_y"], offset["well_x"]], axes=("y", "x"))
                transformations[plate] = Translation([offset["plate_y"], offset["plate_x"]], axes=("y", "x"))
            images[f"{fov}_image"] = Image2DModel.parse(
                _read_fov(root, batch, load_data.loc[(well, site)], channels),
                dims=("c", "y", "x"),
                c_coords=channels,
                transformations=transformations,
                scale_factors=[2, 2],
            )
            directory = _site_dir(root, batch, plate, well, site)
            site_labels = _site_labels(directory, well, site) if directory.is_dir() else {}
            if not site_labels:
                continue
            analysed.append(directory / "Cells.csv")
            for name, mask in site_labels.items():
                masks[f"{fov}_{name}"] = mask
                labels[f"{fov}_{name}"] = Labels2DModel.parse(mask, dims=("y", "x"), transformations=transformations)

    tables, shapes = {}, {}
    if analysed:
        tables["cells"] = _cell_table(analysed, {k: v for k, v in masks.items() if k.endswith("_cells")})
    if profile is not None:
        path = (
            Path(profile)
            if isinstance(profile, Path)
            else root / "workspace/profiles" / batch / plate / f"{plate}_{profile}.csv.gz"
        )
        tables["wells"] = _well_table(path, region=f"{plate}_wells" if located else None, plate_format=plate_format)
        if located:
            shapes[f"{plate}_wells"] = _well_shapes(
                load_data, pixel_size=pixel_size, plate_format=plate_format, system=plate
            )
    return SpatialData(images=images, labels=labels, shapes=shapes, tables=tables)
