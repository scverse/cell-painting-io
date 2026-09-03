# Cell Painting Gallery coverage

Every accession in `s3://cellpainting-gallery/` — 44, against the ~27 listed on the [dataset page](https://broadinstitute.github.io/cellpainting-gallery/complete_datasets.html).
Only well-level (or, for the pooled screens, barcode-level) profiles are ever downloaded; no images.

Sizes below are what a notebook actually pulls, not the whole profile tree.
An earlier version of this file excluded several datasets as "too large" by measuring every profile variant; measuring only the variant in use collapsed most of those, and they are now covered.

## Covered

| Accession | Notebook | Downloaded | Notes |
| --- | --- | --- | --- |
| `cpg0000-jump-pilot` | `cpjump1_to_anndata.ipynb` | 5.5 GB | CPJUMP1, via the paper's GitHub mirror |
| `cpg0001-cellpainting-protocol` | `cellpainting_protocol_to_anndata.ipynb` | 118 MB | Batches share *no* features; one batch used |
| `cpg0002-jump-scope` | `jump_scope_to_anndata.ipynb` | 158 MB | Batches are imaging configurations |
| `cpg0004-lincs` | `lincs_to_anndata.ipynb` | 789 MB | Feature selection per plate, so `_normalized_dmso` |
| `cpg0006-miami` | `miami_to_anndata.ipynb` | 6.7 MB | |
| `cpg0008-pki` | `pki_to_anndata.ipynb` | 11 MB | Kinase inhibitors, dose and MoA |
| `cpg0009-molglue` | `molglue_to_anndata.ipynb` | 267 MB | No `profiles/`; aggregated CSVs under `backend/` |
| `cpg0010-caie-drugresponse` | `caie_to_anndata.ipynb` | 6 MB | Uses the `Image_Metadata_` prefix |
| `cpg0012-wawer-bioactivecompoundprofiling` | `wawer_to_anndata.ipynb` | 386 MB | CDRP, 153,022 wells, the largest per-plate set |
| `cpg0014-jump-adipocyte` | `jump_adipocyte_to_anndata.ipynb` | 120 MB | Differentiated adipocytes |
| `cpg0016-jump-assembled` | `jump_crispr_to_anndata.ipynb` | 177 MB | Assembled CRISPR arm |
| `cpg0017-rohban-pathways` | `rohban_to_anndata.ipynb` | 12 MB | Gene identity outranks plate |
| `cpg0020-varchamp` | `varchamp_to_anndata.ipynb` | ~90 MB | 20 duplicated plate/well rows |
| `cpg0021-periscope` | `periscope_to_anndata.ipynb` | 216 MB | Pooled optical screen; **an observation is a gene**, one gene-level file of 402 GB |
| `cpg0022-cmqtl` | `cmqtl_to_anndata.ipynb` | ~300 MB | iPSC lines |
| `cpg0024-bortezomib` | `bortezomib_to_anndata.ipynb` | 75 MB | Parquet, unnormalized, scaled before PCA |
| `cpg0026-lacoste_haghighi-rare-diseases` | `rare_diseases_to_anndata.ipynb` | 91 MB | `population_profiles/`; cancer mutations batch |
| `cpg0028-kelley-resistance` | `kelley_to_anndata.ipynb` | 809 MB | `backend/` CSVs; duplicated plate/well rows |
| `cpg0029-chroma-pilot` | `chroma_to_anndata.ipynb` | 7.7 MB | Dye sets differ per batch |
| `cpg0031-caicedo-cmvip` | `luad_to_anndata.ipynb` | 26 MB | LUAD alleles, gene outranks plate |
| `cpg0032-pooled-rare` | `pooled_rare_to_anndata.ipynb` | 3.7 MB | Pooled; observation is a barcode |
| `cpg0033-oasis-pilot` | `oasis_pilot_to_anndata.ipynb` | 11 MB | U2OS and HepaRG |
| `cpg0038-tegtmeyer-neuropainting` | `neuropainting_to_anndata.ipynb` | 4.5 MB | Astrocytes and neurons, 20x and 63x |
| `cpg0039-garcia-fossa-livecellpainting` | `garcia_fossa_live_to_anndata.ipynb` | 6.1 MB | Live Cell Painting |
| `cpg0040-garcia-fossa-AgNP` | `garcia_fossa_agnp_to_anndata.ipynb` | 288 KB | Nanoparticle size, dose, time |
| `cpg0047-amish` | `amish_to_anndata.ipynb` | 1.2 MB | Cell line, density, timepoint |

[EUbOPEN](https://zenodo.org/records/10894238) is covered by `eubopen_to_anndata.ipynb`; it is on Zenodo, not the gallery.

## Not covered, and why

Each of these was checked by listing the accession's `workspace/` prefixes, not assumed.

| Accession | What is published | Why not |
| --- | --- | --- |
| `cpg0003-rosetta` | `preprocessed_data/` | Re-processed copies of CDRP, LINCS-Pilot1, LUAD and TA-ORF, all already covered here from their own accessions |
| `cpg0005-gerry-bioactivity` | one 64 MB CSV | Loads as 4,608 x 1,475 but carries **no metadata columns at all** — `obs` is empty, so wells cannot be annotated and no covariate can be tested |
| `cpg0011-lipocyteprofiler` | `profiles/` | **Not well-level**: 225 rows keyed by `patientID` and `cellType`, a different aggregation unit from every other dataset |
| `cpg0015-heterogeneity` | `supplementary/` | No profile tables published |
| `cpg0016-jump` (13 sources) | per-source `profiles/` | Superseded by `cpg0016-jump-assembled`, which is the table the field actually uses and which is covered |
| `cpg0018-singh-seedseq` | nothing under `workspace/` | Images only |
| `cpg0019-moshkov-deepprofiler` | nothing under `workspace/` | Images and DeepProfiler artefacts only |
| `cpg0023-mpi` | `scratch/` | Nothing published |
| `cpg0025-dactyloscopy` | `backend/` | 83 aggregated CSVs totalling 45 GB, roughly 543 MB each; the only dataset still excluded purely on size |
| `cpg0030-gustafsdottir-cellpainting` | `load_data_csv/`, `metadata/` | Image metadata only, no profiles |
| `cpg0034-arevalo-su-motive` | `publication_data/` | A graph dataset, not well-level profiles |
| `cpg0036-EU-OS-bioactives` | `load_data_csv/`, `metadata/` | Image metadata only |
| `cpg0037-oasis` | `profiles/` per source | The `hepatopac` source loads but `obs` holds only `Site_Count`; the assay is covered by `cpg0033-oasis-pilot` |
| `cpg0042-chandrasekaran-jump` | `profiles_assembled/` | Assembled CPJUMP1, duplicating `cpg0000-jump-pilot` |
| `cpg0043-segmentation` | `analysis/`, `pipelines/` | Segmentation outputs, not profiles |
| `cpg0045-ncats-mito` | `load_data_csv/` | Image metadata only |
| `cpg0046-microrna` | `profiles_assembled/*.h5ad` | **Already distributed as AnnData**, so there is nothing here to convert |
| `cpg0049-ipsc-diff-pancreatic-progenitor` | `load_data_csv/`, `metadata/` | Image metadata only |
