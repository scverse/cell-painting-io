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
    from cell_painting_io.spatial import read_plate

__all__ = [
    "CHANNEL_ALIASES",
    "METADATA_PREFIXES",
    "annotate_features",
    "drop_constant_features",
    "drop_extreme_features",
    "drop_incomplete_features",
    "drop_incomplete_wells",
    "neighbour_enrichment",
    "read_plate",
    "read_profiles",
]


def __getattr__(name: str) -> Any:
    # read_plate pulls in spatialdata, which the profile readers do not need
    if name == "read_plate":
        from cell_painting_io.spatial import read_plate

        return read_plate
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
