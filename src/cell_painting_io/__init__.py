from __future__ import annotations

from typing import TYPE_CHECKING, Any

from cell_painting_io.profiles import (
    CHANNEL_ALIASES,
    annotate_features,
    drop_constant_features,
    drop_extreme_features,
    drop_incomplete_features,
    drop_incomplete_wells,
    neighbour_enrichment,
)
from cell_painting_io.reader import METADATA_PREFIXES, read_profiles

if TYPE_CHECKING:
    from cell_painting_io.cellprofiler import cellprofiler_export_plates, read_cellprofiler_export
    from cell_painting_io.spatial import read_plate

__all__ = [
    "CHANNEL_ALIASES",
    "METADATA_PREFIXES",
    "annotate_features",
    "drop_constant_features",
    "drop_extreme_features",
    "drop_incomplete_features",
    "drop_incomplete_wells",
    "cellprofiler_export_plates",
    "neighbour_enrichment",
    "read_cellprofiler_export",
    "read_plate",
    "read_profiles",
]


# The readers below pull in spatialdata, which the profile readers do not need
_SPATIAL = {
    "read_plate": "cell_painting_io.spatial",
    "read_cellprofiler_export": "cell_painting_io.cellprofiler",
    "cellprofiler_export_plates": "cell_painting_io.cellprofiler",
}


def __getattr__(name: str) -> Any:
    if name in _SPATIAL:
        from importlib import import_module

        return getattr(import_module(_SPATIAL[name]), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
