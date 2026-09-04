from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence
from pathlib import Path
from typing import Literal

import anndata as ad
import numpy as np
import pandas as pd

METADATA_PREFIXES: Sequence[str] = ("Image_Metadata_", "Metadata_", "metadata_", "meta_")


def _read_frame(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path)


def _strip_prefix(name: str, prefixes: Sequence[str]) -> str:
    for prefix in sorted(prefixes, key=len, reverse=True):
        if name.startswith(prefix):
            return name[len(prefix) :]
    return name


def read_profiles(
    paths: Path | str | Sequence[Path | str],
    *,
    metadata_prefixes: Sequence[str] = METADATA_PREFIXES,
    metadata_columns: Sequence[str] = (),
    sentinels: float | Collection[float] | None = None,
    index_columns: Sequence[str] | None = None,
    on_column_mismatch: Literal["raise", "intersect"] = "raise",
    path_columns: Mapping[str, int] | None = None,
) -> ad.AnnData:
    """Read CellProfiler profiles into an AnnData of observations x features.

    Numeric columns become the feature matrix, everything else becomes `obs`.

    Args:
        paths: One profile file, or several to stack row-wise. `.parquet` is read as parquet, anything else as CSV.
        metadata_prefixes: Column-name prefixes marking metadata. The matching prefix is stripped from the `obs` column name.
        metadata_columns: Columns to treat as metadata even though they are numeric and unprefixed.
        sentinels: Values in the feature matrix standing for missing, replaced with NaN.
        index_columns: Metadata columns, named as they are after prefix stripping, joined with `:` into the observation index.
        on_column_mismatch: What to do when the files disagree on columns, either raise or keep the intersection in the column order of the first file.
        path_columns: Metadata columns taken from the file path, mapping each column name to how many directories up from the file to read the name of.

    Returns:
        An AnnData whose `X` is float32, `obs` holds the metadata, and `var` is indexed by feature name.

    Raises:
        ValueError: No files were given, the files share no columns, they disagree on columns while `on_column_mismatch` is `"raise"`, or `index_columns` do not identify observations uniquely.
        KeyError: A name in `metadata_columns` or `index_columns` is not in the data.
    """
    files = [Path(paths)] if isinstance(paths, str | Path) else [Path(p) for p in paths]
    if not files:
        raise ValueError("no profile files given")

    frames = [_read_frame(path) for path in files]
    # a file with no rows is read as all-object and would drag the dtypes of the others with it through concat
    kept = [index for index, frame in enumerate(frames) if len(frame)]
    if kept and len(kept) < len(frames):
        files = [files[index] for index in kept]
        frames = [frames[index] for index in kept]
    if len({tuple(frame.columns) for frame in frames}) > 1:
        if on_column_mismatch == "raise":
            raise ValueError("files disagree on columns; pass on_column_mismatch='intersect' to keep the shared ones")
        shared = set.intersection(*(set(frame.columns) for frame in frames))
        if not shared:
            raise ValueError("files share no columns")
        order = [c for c in frames[0].columns if c in shared]
        frames = [frame[order] for frame in frames]
    df = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
    if path_columns:
        # one concat rather than a per-file insert, which fragments a wide frame badly
        lengths = [len(frame) for frame in frames]
        derived = {
            name: np.repeat([path.parents[depth - 1].name for path in files], lengths)
            for name, depth in path_columns.items()
        }
        df = pd.concat([df, pd.DataFrame(derived, index=df.index)], axis=1)

    prefixes = tuple(metadata_prefixes)
    named = set(metadata_columns)
    if missing_named := named - set(df.columns):
        raise KeyError(f"metadata columns not in data: {sorted(missing_named)}")
    prefixed = {c for c in df.columns if c.startswith(prefixes)}
    numeric = set(df.select_dtypes("number").columns)
    meta_columns = [c for c in df.columns if c in prefixed or c in named or c not in numeric]
    feature_columns = [c for c in df.columns if c not in set(meta_columns)]

    x = df[feature_columns].to_numpy(np.float32)
    if sentinels is not None:
        values = [sentinels] if isinstance(sentinels, int | float) else list(sentinels)
        x[np.isin(x, np.asarray(values, dtype=x.dtype))] = np.nan

    obs = df[meta_columns].rename(columns=lambda c: _strip_prefix(c, prefixes))
    for column in obs.select_dtypes("object").columns:
        obs[column] = obs[column].astype("string")
    if index_columns:
        missing = [c for c in index_columns if c not in obs.columns]
        if missing:
            raise KeyError(f"index columns not in metadata: {missing}")
        index = obs[list(index_columns)].astype(str).agg(":".join, axis=1)
        if index.duplicated().any():
            raise ValueError(f"{index_columns} do not identify wells uniquely")
        obs.index = pd.Index(index, name="well")

    var = pd.DataFrame(index=pd.Index(feature_columns, name="feature"))
    return ad.AnnData(X=x, obs=obs, var=var)
