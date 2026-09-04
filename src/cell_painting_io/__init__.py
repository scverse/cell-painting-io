from __future__ import annotations

from collections.abc import Collection, Iterable, Mapping, Sequence
from pathlib import Path
from typing import Literal

import anndata as ad
import numpy as np
import pandas as pd

METADATA_PREFIXES: Sequence[str] = ("Image_Metadata_", "Metadata_", "metadata_", "meta_")

CHANNEL_ALIASES: Mapping[str, str] = {
    "dna": "dna",
    "hoechst": "dna",
    "rna": "rna",
    "agp": "agp",
    "er": "er",
    "mito": "mito",
    "brightfield": "brightfield",
    "lowzbf": "brightfield_low",
    "bflow": "brightfield_low",
    "highzbf": "brightfield_high",
    "bfhigh": "brightfield_high",
}


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


def drop_incomplete_features(adata: ad.AnnData, *, max_missing: float = 0.0) -> ad.AnnData:
    """Drop features that are missing in too large a fraction of observations.

    Args:
        adata: The profiles to filter.
        max_missing: Largest fraction of non-finite values a feature may have and be kept.

    Returns:
        A copy holding the features that passed.
    """
    missing = (~np.isfinite(adata.X)).mean(axis=0)
    return adata[:, missing <= max_missing].copy()


def drop_incomplete_wells(adata: ad.AnnData, *, max_missing: float = 0.0) -> ad.AnnData:
    """Drop observations that are missing in too large a fraction of features.

    Args:
        adata: The profiles to filter.
        max_missing: Largest fraction of non-finite values an observation may have and be kept.

    Returns:
        A copy holding the observations that passed.
    """
    missing = (~np.isfinite(adata.X)).mean(axis=1)
    return adata[missing <= max_missing].copy()


def drop_extreme_features(adata: ad.AnnData, *, max_abs: float = 1e6) -> ad.AnnData:
    """Drop features whose magnitude has blown up, which happens to ratios with a near-zero denominator.

    Args:
        adata: The profiles to filter.
        max_abs: Largest absolute value a feature may reach and be kept.

    Returns:
        A copy holding the features that passed.
    """
    largest = np.nanmax(np.abs(adata.X), axis=0)
    return adata[:, largest <= max_abs].copy()


def drop_constant_features(adata: ad.AnnData) -> ad.AnnData:
    """Drop features with zero variance, which carry no information and break scaling.

    Args:
        adata: The profiles to filter.

    Returns:
        A copy holding the features that vary.
    """
    return adata[:, np.nanstd(adata.X, axis=0) > 0].copy()


def annotate_features(var: pd.DataFrame, *, aliases: Mapping[str, str] = CHANNEL_ALIASES) -> None:
    """Annotate a feature table in place with what each CellProfiler feature name encodes.

    Names follow `<compartment>_<family>_<measurement>[_<channel>][_<parameters>]`.
    Compartment and family are read positionally, channels are the name tokens that match a known channel.
    `Correlation` features measure colocalization between a pair of channels, hence two channel columns; `AreaShape` and `Neighbors` are geometry and have none.
    The measurement itself and its parameters are not extracted.

    Args:
        var: The feature table to annotate, indexed by feature name. Gains the columns `compartment`, `family`, `channel`, `channel_2` and `n_channels`.
        aliases: Maps a lower-cased name token to the channel it stands for, so that datasets naming a channel differently stay comparable.
    """
    tokens = [name.lower().split("_") for name in var.index]
    matched = [[aliases[t] for t in token if t in aliases] for token in tokens]
    var["compartment"] = pd.Categorical([t[0] for t in tokens])
    var["family"] = pd.Categorical([t[1] if len(t) > 1 else "" for t in tokens])
    var["channel"] = pd.Categorical([c[0] if c else None for c in matched])
    var["channel_2"] = pd.Categorical([c[1] if len(c) > 1 else None for c in matched])
    var["n_channels"] = np.array([len(c) for c in matched], dtype=np.int8)


def neighbour_enrichment(adata: ad.AnnData, keys: Iterable[str]) -> pd.DataFrame:
    """Measure how much more often neighbours in the kNN graph share a label than chance would give.

    An embedding of Cell Painting profiles is easy to over-read, so this quantifies the structure before it is looked at.
    The baseline is the rate expected if the labels were shuffled, which is the sum of squared label frequencies.

    Args:
        adata: Profiles that have been through `scanpy.pp.neighbors`, so that `obsp["connectivities"]` exists.
        keys: Names of the `obs` columns to score. Observations with a missing label are ignored for that column.

    Returns:
        One row per covariate, with the observed rate, the shuffled baseline, and their ratio.
    """
    graph = adata.obsp["connectivities"].tocoo()
    rows = []
    for key in keys:
        labels = adata.obs[key]
        codes = labels.astype("category").cat.codes.to_numpy()
        known = labels.notna().to_numpy()
        edges = known[graph.row] & known[graph.col]
        observed = float((codes[graph.row][edges] == codes[graph.col][edges]).mean())
        baseline = float((labels.value_counts(normalize=True).to_numpy() ** 2).sum())
        rows.append({"covariate": key, "observed": observed, "baseline": baseline, "ratio": observed / baseline})
    if not rows:
        return pd.DataFrame(columns=["observed", "baseline", "ratio"], index=pd.Index([], name="covariate"))
    return pd.DataFrame(rows).set_index("covariate").round(3)
