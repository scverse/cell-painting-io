# Coverage

## Cell Painting Gallery

Every accession in `s3://cellpainting-gallery/` — 44, against the ~27 listed on the [dataset page](https://broadinstitute.github.io/cellpainting-gallery/complete_datasets.html).
Only well-level (or, for the pooled screens, barcode-level) profiles are downloaded, with one exception: `cpg0000-jump-pilot` also has a notebook that reads images.

Sizes below are what a notebook actually pulls, not the whole profile tree.
An earlier version of this file excluded several datasets as "too large" by measuring every profile variant rather than the one a notebook reads; every one of those exclusions turned out to be wrong and they are now covered.
The remaining exclusions are structural: the accession publishes no profiles, duplicates one already covered, or already ships AnnData.

## Covered

| Accession | Notebook | Downloaded | Notes |
| --- | --- | --- | --- |
| `cpg0000-jump-pilot` | `notebooks/cpjump1_to_anndata.ipynb` | 5.5 GB | CPJUMP1, via the paper's GitHub mirror |
| `cpg0000-jump-pilot` | `notebooks/cpjump1_to_spatialdata.ipynb` | 700 MB | Images and segmentations of two wells on each of two plates |
| `cpg0001-cellpainting-protocol` | `notebooks/cellpainting_protocol_to_anndata.ipynb` | 118 MB | Batches share *no* features; one batch used |
| `cpg0002-jump-scope` | `notebooks/jump_scope_to_anndata.ipynb` | 158 MB | Batches are imaging configurations |
| `cpg0004-lincs` | `notebooks/lincs_to_anndata.ipynb` | 789 MB | Feature selection per plate, so `_normalized_dmso` |
| `cpg0006-miami` | `notebooks/miami_to_anndata.ipynb` | 6.7 MB | |
| `cpg0008-pki` | `notebooks/pki_to_anndata.ipynb` | 11 MB | Kinase inhibitors, dose and MoA |
| `cpg0009-molglue` | `notebooks/molglue_to_anndata.ipynb` | 267 MB | No `profiles/`; aggregated CSVs under `backend/` |
| `cpg0010-caie-drugresponse` | `notebooks/caie_to_anndata.ipynb` | 6 MB | Uses the `Image_Metadata_` prefix |
| `cpg0012-wawer-bioactivecompoundprofiling` | `notebooks/wawer_to_anndata.ipynb` | 386 MB | CDRP, 153,022 wells, the largest per-plate set |
| `cpg0014-jump-adipocyte` | `notebooks/jump_adipocyte_to_anndata.ipynb` | 120 MB | Differentiated adipocytes |
| `cpg0016-jump-assembled` | `notebooks/jump_crispr_to_anndata.ipynb` | 177 MB | Assembled CRISPR arm |
| `cpg0017-rohban-pathways` | `notebooks/rohban_to_anndata.ipynb` | 12 MB | Gene identity outranks plate |
| `cpg0020-varchamp` | `notebooks/varchamp_to_anndata.ipynb` | ~90 MB | 20 duplicated plate/well rows |
| `cpg0021-periscope` | `notebooks/periscope_to_anndata.ipynb` | 216 MB | Pooled optical screen; **an observation is a gene**, one gene-level file of 402 GB |
| `cpg0022-cmqtl` | `notebooks/cmqtl_to_anndata.ipynb` | ~300 MB | iPSC lines |
| `cpg0024-bortezomib` | `notebooks/bortezomib_to_anndata.ipynb` | 75 MB | Parquet, unnormalized, scaled before PCA |
| `cpg0026-lacoste_haghighi-rare-diseases` | `notebooks/rare_diseases_to_anndata.ipynb` | 91 MB | `population_profiles/`; cancer mutations batch |
| `cpg0028-kelley-resistance` | `notebooks/kelley_to_anndata.ipynb` | 809 MB | `backend/` CSVs; duplicated plate/well rows |
| `cpg0029-chroma-pilot` | `notebooks/chroma_to_anndata.ipynb` | 7.7 MB | Dye sets differ per batch |
| `cpg0031-caicedo-cmvip` | `notebooks/luad_to_anndata.ipynb` | 26 MB | LUAD alleles, gene outranks plate |
| `cpg0032-pooled-rare` | `notebooks/pooled_rare_to_anndata.ipynb` | 3.7 MB | Pooled; observation is a barcode |
| `cpg0033-oasis-pilot` | `notebooks/oasis_pilot_to_anndata.ipynb` | 11 MB | U2OS and HepaRG |
| `cpg0005-gerry-bioactivity` | `notebooks/gerry_to_anndata.ipynb` | 62 MB | One flat CSV; the only non-feature column is an opaque `well_profile_id`, so nothing can be tested against the embedding |
| `cpg0011-lipocyteprofiler` | `notebooks/lipocyte_to_anndata.ipynb` | 16 MB | **Observations are not wells** — aggregated per patient and cell type, 225 rows |
| `cpg0025-dactyloscopy` | `notebooks/dactyloscopy_to_anndata.ipynb` | 264 MB | `backend/` CSVs; five cell lines |
| `cpg0037-oasis` | `notebooks/oasis_to_anndata.ipynb` | 52 MB | `axiom` source; 299 duplicated plate/well rows |
| `cpg0038-tegtmeyer-neuropainting` | `notebooks/neuropainting_to_anndata.ipynb` | 4.5 MB | Astrocytes and neurons, 20x and 63x |
| `cpg0039-garcia-fossa-livecellpainting` | `notebooks/garcia_fossa_live_to_anndata.ipynb` | 6.1 MB | Live Cell Painting |
| `cpg0040-garcia-fossa-AgNP` | `notebooks/garcia_fossa_agnp_to_anndata.ipynb` | 288 KB | Nanoparticle size, dose, time |
| `cpg0047-amish` | `notebooks/amish_to_anndata.ipynb` | 1.2 MB | Cell line, density, timepoint |

[EUbOPEN](https://zenodo.org/records/10894238) is covered by `notebooks/eubopen_to_anndata.ipynb`; it is on Zenodo, not the gallery.

## Not covered, and why

Each of these was checked by listing the accession's `workspace/` prefixes, not assumed.

| Accession | What is published | Why not |
| --- | --- | --- |
| `cpg0003-rosetta` | `preprocessed_data/` | Re-processed copies of CDRP, LINCS-Pilot1, LUAD and TA-ORF, all already covered here from their own accessions |
| `cpg0015-heterogeneity` | `supplementary/` | No profile tables published |
| `cpg0016-jump` (13 sources) | per-source `profiles/` | Superseded by `cpg0016-jump-assembled`, which is the table the field actually uses and which is covered |
| `cpg0018-singh-seedseq` | nothing under `workspace/` | Images only |
| `cpg0019-moshkov-deepprofiler` | nothing under `workspace/` | Images and DeepProfiler artefacts only |
| `cpg0023-mpi` | `scratch/` | Nothing published |
| `cpg0030-gustafsdottir-cellpainting` | `load_data_csv/`, `metadata/` | Image metadata only, no profiles |
| `cpg0034-arevalo-su-motive` | `publication_data/` | A graph dataset, not well-level profiles |
| `cpg0036-EU-OS-bioactives` | `load_data_csv/`, `metadata/` | Image metadata only |
| `cpg0042-chandrasekaran-jump` | `profiles_assembled/` | Assembled CPJUMP1, duplicating `cpg0000-jump-pilot` |
| `cpg0043-segmentation` | `analysis/`, `pipelines/` | Segmentation outputs, not profiles |
| `cpg0045-ncats-mito` | `load_data_csv/` | Image metadata only |
| `cpg0046-microrna` | `profiles_assembled/*.h5ad` | **Already distributed as AnnData**, so there is nothing here to convert |
| `cpg0049-ipsc-diff-pancreatic-progenitor` | `load_data_csv/`, `metadata/` | Image metadata only |


## Image Data Resource

Seven studies in the [IDR](https://idr.openmicroscopy.org/) use the Cell Painting assay, identified from the `Study Type`, keywords and dye lists in each study's metadata record rather than from the short API descriptions, which omit the term for several of them.

IDR is an image repository: each study publishes one ISA-Tab annotation table per screen, and only one puts features in it.
The rest carry experimental annotation only; where their features exist at all they are in the Cell Painting Gallery, already covered here.

Annotation tables are read from each study's GitHub repository rather than the IDR API, to keep load off a shared public server.

| Study | Annotation columns | Features | Status |
| --- | --- | --- | --- |
| `idr0133-dahlin-cellpainting` | 403 | 372 | `notebooks/idr0133_to_anndata.ipynb` |
| `idr0016-wawer-bioactivecompoundprofiling` | 24 | 0 | annotation only; features are `cpg0012` |
| `idr0033-rohban-pathways` | 77 | 0 | annotation only; features are `cpg0017` |
| `idr0088-cox-phenomicprofiling` | 21 | 0 | annotation only |
| `idr0080-way-perturbation` | 17 | 0 | annotation only |
| `idr0160-lippincott-pyroptosis` | 38 | 0 | annotation only |
| `idr0036-gustafsdottir-cellpainting` | — | — | no annotation repository on GitHub; publishes no profiles here or in the gallery (`cpg0030`) |

`idr0035-caie-drugresponse` is BBBC021 rather than Cell Painting proper and is covered from `cpg0010`.

`idr0093-mueller-perturbation` publishes ~380 per-well features, but under its own naming (`nuclei_area_mean`, `frac_G1`) rather than CellProfiler's, and it is an EU/nascent-RNA assay rather than Cell Painting, so it is out of scope.
Detection here keys on CellProfiler naming, so features published under another convention would be missed.


## BioImage Archive

The [BioImage Archive](https://www.ebi.ac.uk/bioimage-archive/) hosts raw Cell Painting images, and a few studies also deposit the CellProfiler feature tables.
Where those tables appear nowhere else — no gallery accession, no Zenodo record — the study is covered here.
Candidates were found by searching the BioImages collection for the Cell Painting assay and walking each study's public file tree for well- or object-level `.parquet` / `.csv` matrices.

| Accession | What it publishes | Status |
| --- | --- | --- |
| `S-BIAD2254` | 3D colorectal-spheroid Cell Painting (HCT116, HT29), per-well and per-object profiles plus optically-cleared z-stacks | `notebooks/spheroid_3d_to_anndata.ipynb` |
| `S-BIAD2262` | EUbOPEN compounds Cell Painting, per-plate well profiles plus images | Profiles duplicate `notebooks/eubopen_to_anndata.ipynb`, the Zenodo deposit of the same Servier U2OS screen |
| `S-BIAD1094` | EUbOPEN Wave 1 compounds Cell Painting | Images only, no feature tables |
| `S-BIAD0847`, `S-BIAD0848`, `S-BIAD0851`, `S-BIAD0855` | OME-NGFF mirrors of IDR Cell Painting screens | Images only; where features exist they are in the gallery |
