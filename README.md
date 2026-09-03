# cell-painting-io

Experimental! A collection of scripts &amp; utilities to get Cell Painting datasets into scverse data structures.

`cell_painting_io.py` holds the dataset-agnostic parts — reading profiles into an `AnnData`, masking missing-value sentinels, parsing CellProfiler feature names into `var`, and measuring how much of an embedding is explained by experimental design.
Each notebook adds only what its dataset needs on top.

Everything that has differed between datasets so far is a parameter rather than an assumption: the missing-value sentinel (`-9999` in EUbOPEN, absent elsewhere), the metadata column prefix (`Metadata_`, `meta_`, `metadata_`), and channel naming (`LowZBF` and `bflow` canonicalise to the same channel).
`read_profiles(..., on_column_mismatch="intersect")` handles collections of plates that disagree on columns, and `drop_incomplete_features` / `drop_extreme_features` remove features that are missing or numerically blown up.

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

### `lincs_to_anndata.ipynb`

LINCS ([`cpg0004-lincs`](https://broadinstitute.github.io/cellpainting-gallery/), idr0125): A549, 1,571 compounds at six doses, batch `2016_04_01_a549_48hr_batch1`.
52,223 wells over 136 plates.
Uses `_normalized_dmso` rather than the feature-selected variants, because LINCS ran feature selection per *plate* and the intersection over 270 plates is 30 features.

```bash
aws s3 cp --no-sign-request --recursive --exclude "*" --include "*_normalized_dmso.csv.gz" \
  s3://cellpainting-gallery/cpg0004-lincs/broad/workspace/profiles/2016_04_01_a549_48hr_batch1/ .
```

### `rohban_to_anndata.ipynb`

Rohban pathways ([`cpg0017-rohban-pathways`](https://broadinstitute.github.io/cellpainting-gallery/), TA-ORF, BBBC037, idr0033): U2OS, 323 genes overexpressed by ORF, five plates.

```bash
aws s3 cp --no-sign-request --recursive --exclude "*" --include "*_normalized_feature_select_batch.csv.gz" \
  s3://cellpainting-gallery/cpg0017-rohban-pathways/broad/workspace/profiles/ .
```

### `luad_to_anndata.ipynb`

Caicedo CMVIP ([`cpg0031-caicedo-cmvip`](https://broadinstitute.github.io/cellpainting-gallery/), LUAD, BBBC043): A549, 596 alleles across 53 lung adenocarcinoma genes, sixteen plates, wild-type and mutant alleles side by side.

```bash
aws s3 cp --no-sign-request --recursive --exclude "*" --include "*_normalized_feature_select_batch.csv.gz" \
  s3://cellpainting-gallery/cpg0031-caicedo-cmvip/broad/workspace/profiles/ .
```

Point the path constant at the top of each notebook at your download.
None of these notebooks need the images, only the well-level profiles.

## Environment

```bash
mamba create -n cell-painting-io -c conda-forge python=3.14 anndata scanpy pandas numpy umap-learn leidenalg matplotlib jupyterlab
```

## Development

```bash
prek install
```
