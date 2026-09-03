from __future__ import annotations

from collections.abc import Collection, Iterable, Mapping, Sequence
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

METADATA_PREFIXES: Sequence[str] = ("Metadata_", "metadata_", "meta_")

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
    sentinels: float | Collection[float] | None = None,
    index_columns: Sequence[str] | None = None,
) -> ad.AnnData:
    files = [Path(paths)] if isinstance(paths, str | Path) else [Path(p) for p in paths]
    if not files:
        raise ValueError("no profile files given")

    frames = [_read_frame(path) for path in files]
    if len({tuple(frame.columns) for frame in frames}) > 1:
        raise ValueError("files disagree on columns; intersect features before concatenating")
    df = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]

    prefixes = tuple(metadata_prefixes)
    meta_columns = [c for c in df.columns if c.startswith(prefixes)]
    feature_columns = [c for c in df.columns if c not in meta_columns]

    x = df[feature_columns].to_numpy(np.float32)
    if sentinels is not None:
        values = [sentinels] if isinstance(sentinels, int | float) else list(sentinels)
        x[np.isin(x, np.asarray(values, dtype=x.dtype))] = np.nan

    obs = df[meta_columns].rename(columns=lambda c: _strip_prefix(c, prefixes))
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
    missing = np.isnan(adata.X).mean(axis=0)
    return adata[:, missing <= max_missing].copy()


def annotate_features(var: pd.DataFrame, *, aliases: Mapping[str, str] = CHANNEL_ALIASES) -> None:
    tokens = [name.lower().split("_") for name in var.index]
    matched = [[aliases[t] for t in token if t in aliases] for token in tokens]
    var["compartment"] = pd.Categorical([t[0] for t in tokens])
    var["family"] = pd.Categorical([t[1] if len(t) > 1 else "" for t in tokens])
    var["channel"] = pd.Categorical([c[0] if c else None for c in matched])
    var["channel_2"] = pd.Categorical([c[1] if len(c) > 1 else None for c in matched])
    var["n_channels"] = np.array([len(c) for c in matched], dtype=np.int8)


def categorize(adata: ad.AnnData, columns: Iterable[str]) -> None:
    for column in columns:
        adata.obs[column] = adata.obs[column].astype("category")


def neighbour_enrichment(adata: ad.AnnData, keys: Iterable[str]) -> pd.DataFrame:
    graph = adata.obsp["connectivities"].tocoo()
    rows = []
    for key in keys:
        labels = adata.obs[key].to_numpy()
        observed = float((labels[graph.row] == labels[graph.col]).mean())
        baseline = float((adata.obs[key].value_counts(normalize=True).to_numpy() ** 2).sum())
        rows.append({"covariate": key, "observed": observed, "baseline": baseline, "ratio": observed / baseline})
    return pd.DataFrame(rows).set_index("covariate").round(3)
