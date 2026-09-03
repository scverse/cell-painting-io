# Cell Painting Gallery coverage

Every accession in `s3://cellpainting-gallery/` as of this survey — 45, against the ~27 listed on the [dataset page](https://broadinstitute.github.io/cellpainting-gallery/complete_datasets.html).
Sizes are the profile prefixes only, never images.

Where a dataset is not covered the reason is stated.
"Too large" means the profiles alone exceed what is reasonable to pull onto one machine for a demo notebook, not that the reader cannot handle it.

## Covered

| Accession | Notebook | Notes |
| --- | --- | --- |
| `cpg0000-jump-pilot` | `cpjump1_to_anndata.ipynb` | CPJUMP1, taken from the paper's GitHub mirror of the same profiles |
| `cpg0004-lincs` | `lincs_to_anndata.ipynb` | Feature selection ran per plate, so `_normalized_dmso` is used |
| `cpg0006-miami` | `miami_to_anndata.ipynb` | |
| `cpg0008-pki` | `pki_to_anndata.ipynb` | Kinase inhibitors with dose and MoA |
| `cpg0014-jump-adipocyte` | `jump_adipocyte_to_anndata.ipynb` | Differentiated adipocytes |
| `cpg0016-jump-assembled` | `jump_crispr_to_anndata.ipynb` | The assembled CRISPR arm, 13 imaging sites |
| `cpg0017-rohban-pathways` | `rohban_to_anndata.ipynb` | Gene identity outranks plate |
| `cpg0020-varchamp` | `varchamp_to_anndata.ipynb` | 20 wells are duplicated under one plate/well id |
| `cpg0024-bortezomib` | `bortezomib_to_anndata.ipynb` | Parquet, unnormalized, so it is scaled before PCA |
| `cpg0026-lacoste_haghighi-rare-diseases` | `rare_diseases_to_anndata.ipynb` | Published as `population_profiles`; cancer mutations batch |
| `cpg0029-chroma-pilot` | `chroma_to_anndata.ipynb` | Alternative dye sets differ per batch |
| `cpg0031-caicedo-cmvip` | `luad_to_anndata.ipynb` | LUAD alleles, gene outranks plate |
| `cpg0033-oasis-pilot` | `oasis_pilot_to_anndata.ipynb` | U2OS and HepaRG assay development |
| `cpg0038-tegtmeyer-neuropainting` | `neuropainting_to_anndata.ipynb` | Astrocytes and neurons, 20x and 63x |
| `cpg0047-amish` | `amish_to_anndata.ipynb` | Cell line, seeding density and timepoint |

[EUbOPEN](https://zenodo.org/records/10894238) is also covered by `eubopen_to_anndata.ipynb`; it is on Zenodo, not in the gallery.

## Has profiles, not covered

| Accession | Profiles | Why not |
| --- | --- | --- |
| `cpg0005-gerry-bioactivity` | 64 MB | Loads as 4,608 x 1,475, but the file carries **no metadata columns at all** — `obs` is empty, so there is nothing to annotate wells with or colour a UMAP by. Would need the platemap joined from `metadata/`. |
| `cpg0011-lipocyteprofiler` | 16 MB | **Not well-level.** 225 rows with `patientID` and `cellType` only, so it is a different aggregation unit from every other dataset here. |
| `cpg0037-oasis` | 2.6 MB – 3.3 GB | The small `hepatopac` source loads (192 x 624) but `obs` holds only `Site_Count`; the informative sources (`axiom`, `xellar`, `insphero`) are 1.2–3.3 GB. Covered in spirit by `cpg0033-oasis-pilot`. |
| `cpg0001-cellpainting-protocol` | 5.1 GB | Too large |
| `cpg0002-jump-scope` | 4.5 GB | Too large |
| `cpg0012-wawer-bioactivecompoundprofiling` | 8.6 GB | Too large |
| `cpg0016-jump` (13 sources) | 1.7–10.5 GB each | Superseded by `cpg0016-jump-assembled`, which is the table people actually use |
| `cpg0022-cmqtl` | 4.7 GB | Too large |
| `cpg0032-pooled-rare` | 23 GB | Too large, and pooled optical screening rather than arrayed Cell Painting |
| `cpg0039-garcia-fossa-livecellpainting` | 10.9 GB | Too large |
| `cpg0040-garcia-fossa-AgNP` | 12.7 GB | Too large |
| `cpg0021-periscope` | 758 GB | Far too large, and a pooled in-situ-sequencing variant, not well-level Cell Painting |

## No well-level profiles published

Checked by listing each `*/workspace/` prefix.

| Accession | What is published | Why not |
| --- | --- | --- |
| `cpg0003-rosetta` | `preprocessed_data/` | Re-processed copies of CDRP, LINCS-Pilot1, LUAD and TA-ORF — datasets already covered here |
| `cpg0042-chandrasekaran-jump` | `profiles_assembled/` | Assembled CPJUMP1; duplicates `cpg0000-jump-pilot`, and the compound table is 2.8 GB |
| `cpg0046-microrna` | `profiles_assembled/*.h5ad` | **Already distributed as AnnData** (2.2–5.9 GB per cell line), so there is nothing for this repo to convert |
| `cpg0009-molglue` | `backend/` | Only per-plate backend databases, no aggregated profiles |
| `cpg0010-caie-drugresponse` | `backend/` | Only backend databases |
| `cpg0025-dactyloscopy` | `backend/` | Only backend databases |
| `cpg0028-kelley-resistance` | `backend/` | Only backend databases |
| `cpg0015-heterogeneity` | `supplementary/` | No profile tables |
| `cpg0034-arevalo-su-motive` | `publication_data/` | A graph/network dataset, not well-level profiles |
| `cpg0018-singh-seedseq` | nothing under `workspace/` | Images only |
| `cpg0019-moshkov-deepprofiler` | nothing under `workspace/` | Images and DeepProfiler artefacts only |
| `cpg0023-mpi` | `scratch/` | Nothing published |
| `cpg0030-gustafsdottir-cellpainting` | `load_data_csv/`, `metadata/` | Image metadata only, no profiles |
| `cpg0036-EU-OS-bioactives` | `load_data_csv/`, `metadata/` | Image metadata only |
| `cpg0043-segmentation` | `analysis/`, `pipelines/` | Segmentation outputs, not profiles |
| `cpg0045-ncats-mito` | `load_data_csv/` | Image metadata only |
| `cpg0049-ipsc-diff-pancreatic-progenitor` | `load_data_csv/`, `metadata/` | Image metadata only |
