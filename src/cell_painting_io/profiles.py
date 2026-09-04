from __future__ import annotations

from collections.abc import Iterable, Mapping

import anndata as ad
import numpy as np
import pandas as pd

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
