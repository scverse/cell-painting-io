# cell-painting-io

Experimental! A collection of scripts &amp; utilities to get Cell Painting datasets into scverse data structures.

`cell_painting_io.py` holds the dataset-agnostic parts — reading profiles into an `AnnData`, masking missing-value sentinels, parsing CellProfiler feature names into `var`, and measuring how much of an embedding is explained by experimental design.
Each notebook adds only what its dataset needs on top.

## Notebooks

### `cpjump1_to_anndata.ipynb`

The JUMP-Cell Painting pilot dataset ([Chandrasekaran et al., *Nature Methods* 2024](https://www.nature.com/articles/s41592-024-02241-6)), batch `2020_11_04_CPJUMP1`: 19,498 wells x 903 features over 51 plates.
Profiles are per-plate `.csv.gz` with a `Metadata_` prefix, already pycytominer-normalized and feature selected.

Get the source data — a plain `git clone` leaves the profiles as Git LFS pointer stubs, so `git lfs pull` is required:

```bash
git clone https://github.com/jump-cellpainting/2024_Chandrasekaran_NatureMethods_CPJUMP1
cd 2024_Chandrasekaran_NatureMethods_CPJUMP1 && git lfs pull
```

### `eubopen_to_anndata.ipynb`

The EUbOPEN compound dataset ([Zenodo 10894238](https://zenodo.org/records/10894238), CC-BY-4.0): 39,206 wells over 105 plates in four runs, U2OS, 575 perturbagen ids at 8 doses.
A single flat parquet with a `meta_`/`metadata_` prefix, unnormalized, and with missing values written as `-9999` rather than NaN — 24% of the matrix, leaving 4,943 complete features of 7,703.

```bash
curl -L -o hcs_cellpainting_eubopen.parquet \
  "https://zenodo.org/api/records/10894238/files/hcs_cellpainting_eubopen(1).parquet/content"
```

Point the path constant at the top of each notebook at your download.

## Environment

```bash
mamba create -n cell-painting-io -c conda-forge python=3.14 anndata scanpy pandas numpy umap-learn leidenalg matplotlib jupyterlab
```

## Development

```bash
prek install
```
