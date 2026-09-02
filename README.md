# cell-painting-io

Experimental! A collection of scripts &amp; utilities to get Cell Painting datasets into scverse data structures.

## Notebooks

### `cpjump1_to_anndata.ipynb`

Builds an `AnnData` object from the CellProfiler well-level profiles of the JUMP-Cell Painting pilot dataset ([Chandrasekaran et al., *Nature Methods* 2024](https://www.nature.com/articles/s41592-024-02241-6)), and runs a short UMAP analysis on it.
One observation is a well, one variable is a CellProfiler feature.

Get the source data — note that a plain `git clone` leaves the profiles as Git LFS pointer stubs, so `git lfs pull` is required:

```bash
git clone https://github.com/jump-cellpainting/2024_Chandrasekaran_NatureMethods_CPJUMP1
cd 2024_Chandrasekaran_NatureMethods_CPJUMP1 && git lfs pull
```

Then point `DATA_ROOT` in the notebook at that checkout.

Environment:

```bash
mamba create -n cell-painting-io -c conda-forge python=3.14 anndata scanpy pandas numpy umap-learn leidenalg matplotlib jupyterlab
```

## Development

```bash
prek install
```
