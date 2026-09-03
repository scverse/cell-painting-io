# cell-painting-io

Experimental! A collection of scripts &amp; utilities to get Cell Painting datasets into scverse data structures.

`cell_painting_io.py` holds the dataset-agnostic parts — reading profiles into an `AnnData`, masking missing-value sentinels, parsing CellProfiler feature names into `var`, and measuring how much of an embedding is explained by experimental design.
Each notebook adds only what its dataset needs on top.

Everything that has differed between datasets so far is a parameter rather than an assumption: the missing-value sentinel (`-9999` in EUbOPEN, absent elsewhere), the metadata column prefix (`Metadata_`, `meta_`, `metadata_`), and channel naming (`LowZBF` and `bflow` canonicalise to the same channel).
`read_profiles(..., on_column_mismatch="intersect")` handles collections of plates that disagree on columns, and `drop_incomplete_features` / `drop_extreme_features` remove features that are missing or numerically blown up.

[COVERAGE.md](COVERAGE.md) lists every accession in the Cell Painting Gallery with its status, and a reason for each one that is not covered.

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

### The rest

One notebook per remaining dataset, each following the same shape: read, annotate, drop unusable features, save, embed, then measure what the embedding actually encodes before looking at it.

| Notebook | Dataset | Wells x features |
| --- | --- | --- |
| `jump_crispr_to_anndata.ipynb` | JUMP CRISPR, assembled (`cpg0016-jump-assembled`) | 51,185 x 599 |
| `pki_to_anndata.ipynb` | Kinase inhibitors (`cpg0008`) | 3,072 x 1,061 |
| `varchamp_to_anndata.ipynb` | Coding variants (`cpg0020`) | 7,177 x 853 |
| `rare_diseases_to_anndata.ipynb` | Rare disease variants (`cpg0026`) | 5,634 x 1,795 |
| `jump_adipocyte_to_anndata.ipynb` | Adipocytes (`cpg0014`) | 9,216 x 848 |
| `oasis_pilot_to_anndata.ipynb` | OASIS pilot (`cpg0033`) | 2,303 x 718 |
| `bortezomib_to_anndata.ipynb` | Bortezomib resistance (`cpg0024`) | 768 x 4,695 |
| `miami_to_anndata.ipynb` | MIAMI (`cpg0006`) | 1,387 x 264 |
| `chroma_to_anndata.ipynb` | Alternative dyes (`cpg0029`) | 767 x 675 |
| `amish_to_anndata.ipynb` | Amish cohort (`cpg0047`) | 468 x 790 |
| `neuropainting_to_anndata.ipynb` | Brain cell types (`cpg0038`) | 188 x 563 |
| `wawer_to_anndata.ipynb` | CDRP bioactives (`cpg0012`) | 153,022 x 781 |
| `periscope_to_anndata.ipynb` | PERISCOPE pooled screen (`cpg0021`), genes not wells | 20,384 x 3,640 |
| `molglue_to_anndata.ipynb` | Molecular glues (`cpg0009`), from `backend/` | 1,920 x 7,587 |
| `cmqtl_to_anndata.ipynb` | iPSC lines (`cpg0022`) | 4,148 x 1,327 |
| `cellpainting_protocol_to_anndata.ipynb` | Protocol variants (`cpg0001`) | 5,794 x 187 |
| `jump_scope_to_anndata.ipynb` | JUMP-SCOPE microscopes (`cpg0002`) | 1,529 x 819 |
| `kelley_to_anndata.ipynb` | Bortezomib resistance (`cpg0028`) | 1,920 x 221 |
| `caie_to_anndata.ipynb` | Caie drug response (`cpg0010`) | 632 x 516 |
| `garcia_fossa_live_to_anndata.ipynb` | Live Cell Painting (`cpg0039`) | 536 x 82 |
| `pooled_rare_to_anndata.ipynb` | Pooled rare variants (`cpg0032`), barcodes | 290 x 4,046 |
| `garcia_fossa_agnp_to_anndata.ipynb` | Silver nanoparticles (`cpg0040`) | 180 x 112 |

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
