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
| `cpg0000-jump-pilot` | `notebooks/anndata/cpjump1.ipynb` | 5.5 GB | CPJUMP1, via the paper's GitHub mirror |
| `cpg0000-jump-pilot` | `notebooks/spatialdata/cpjump1.ipynb` | 700 MB | Images and segmentations of two wells on each of two plates |
| `cpg0001-cellpainting-protocol` | `notebooks/anndata/cellpainting_protocol.ipynb` | 118 MB | Batches share *no* features; one batch used |
| `cpg0001-cellpainting-protocol` | `notebooks/spatialdata/cellpainting_protocol.ipynb` | one well | Only some fields of a well were analysed |
| `cpg0002-jump-scope` | `notebooks/anndata/jump_scope.ipynb` | 158 MB | Batches are imaging configurations |
| `cpg0002-jump-scope` | `notebooks/spatialdata/jump_scope.ipynb` | one well | Profile under `backend/`; no stage coordinates |
| `cpg0004-lincs` | `notebooks/anndata/lincs.ipynb` | 789 MB | Feature selection per plate, so `_normalized_dmso` |
| `cpg0006-miami` | `notebooks/anndata/miami.ipynb` | 6.7 MB | |
| `cpg0006-miami` | `notebooks/spatialdata/miami.ipynb` | one well | Images and segmentations |
| `cpg0008-pki` | `notebooks/anndata/pki.ipynb` | 11 MB | Kinase inhibitors, dose and MoA |
| `cpg0008-pki` | `notebooks/spatialdata/pki.ipynb` | one well | Images and segmentations |
| `cpg0009-molglue` | `notebooks/anndata/molglue.ipynb` | 267 MB | No `profiles/`; aggregated CSVs under `backend/` |
| `cpg0010-caie-drugresponse` | `notebooks/anndata/caie.ipynb` | 6 MB | Uses the `Image_Metadata_` prefix |
| `cpg0012-wawer-bioactivecompoundprofiling` | `notebooks/anndata/wawer.ipynb` | 386 MB | CDRP, 153,022 wells, the largest per-plate set |
| `cpg0012-wawer-bioactivecompoundprofiling` | `notebooks/spatialdata/wawer.ipynb` | one well | Lower-case wells; no stage coordinates |
| `cpg0014-jump-adipocyte` | `notebooks/anndata/jump_adipocyte.ipynb` | 120 MB | Differentiated adipocytes |
| `cpg0014-jump-adipocyte` | `notebooks/spatialdata/jump_adipocyte.ipynb` | one well | Profile under `backend/` |
| `cpg0016-jump-assembled` | `notebooks/anndata/jump_crispr.ipynb` | 177 MB | Assembled CRISPR arm |
| `cpg0016-jump` | `notebooks/spatialdata/jump.ipynb`, `notebooks/spatialdata/jump_source1.ipynb` | one well each | `source_4` and the 1536-well `source_1` |
| `cpg0017-rohban-pathways` | `notebooks/anndata/rohban.ipynb` | 12 MB | Gene identity outranks plate |
| `cpg0017-rohban-pathways` | `notebooks/spatialdata/rohban.ipynb` | one well | No stage coordinates, so fields only |
| `cpg0020-varchamp` | `notebooks/anndata/varchamp.ipynb` | ~90 MB | 20 duplicated plate/well rows |
| `cpg0020-varchamp` | `notebooks/spatialdata/varchamp.ipynb` | one well | Outlines drawn in colour, two object types in one image |
| `cpg0021-periscope` | `notebooks/anndata/periscope.ipynb` | 216 MB | Pooled optical screen; **an observation is a gene**, one gene-level file of 402 GB |
| `cpg0022-cmqtl` | `notebooks/anndata/cmqtl.ipynb` | ~300 MB | iPSC lines |
| `cpg0022-cmqtl` | `notebooks/spatialdata/cmqtl.ipynb` | one well | Profile is an uncompressed CSV |
| `cpg0024-bortezomib` | `notebooks/anndata/bortezomib.ipynb` | 75 MB | Parquet, unnormalized, scaled before PCA |
| `cpg0024-bortezomib` | `notebooks/spatialdata/bortezomib.ipynb` | one well | No stage coordinates |
| `cpg0026-lacoste_haghighi-rare-diseases` | `notebooks/anndata/rare_diseases.ipynb` | 91 MB | `population_profiles/`; cancer mutations batch |
| `cpg0026-lacoste_haghighi-rare-diseases` | `notebooks/spatialdata/rare_diseases.ipynb` | one well | Colour outlines; no well-level profile published |
| `cpg0028-kelley-resistance` | `notebooks/anndata/kelley.ipynb` | 809 MB | `backend/` CSVs; duplicated plate/well rows |
| `cpg0028-kelley-resistance` | *(no spatialdata notebook)* | — | `load_data.csv` gives `Metadata_Row`/`Col` transposed against the well names its own analysis directories use, so the well of a field cannot be read off it |
| `cpg0029-chroma-pilot` | `notebooks/anndata/chroma.ipynb` | 7.7 MB | Dye sets differ per batch |
| `cpg0029-chroma-pilot` | `notebooks/spatialdata/chroma.ipynb` | one well | Eight channels; images named by `FileName_`/`PathName_` |
| `cpg0031-caicedo-cmvip` | `notebooks/anndata/luad.ipynb` | 26 MB | LUAD alleles, gene outranks plate |
| `cpg0031-caicedo-cmvip` | `notebooks/spatialdata/caicedo_cmvip.ipynb` | one well | Lower-case wells; no stage coordinates |
| `cpg0032-pooled-rare` | `notebooks/anndata/pooled_rare.ipynb` | 3.7 MB | Pooled; observation is a barcode |
| `cpg0033-oasis-pilot` | `notebooks/anndata/oasis_pilot.ipynb` | 11 MB | U2OS and HepaRG |
| `cpg0033-oasis-pilot` | `notebooks/spatialdata/oasis_pilot.ipynb` | one well | Images and segmentations |
| `cpg0005-gerry-bioactivity` | `notebooks/anndata/gerry.ipynb` | 62 MB | One flat CSV; the only non-feature column is an opaque `well_profile_id`, so nothing can be tested against the embedding |
| `cpg0011-lipocyteprofiler` | `notebooks/anndata/lipocyte.ipynb` | 16 MB | **Observations are not wells** — aggregated per patient and cell type, 225 rows |
| `cpg0025-dactyloscopy` | `notebooks/anndata/dactyloscopy.ipynb` | 264 MB | `backend/` CSVs; five cell lines |
| `cpg0025-dactyloscopy` | `notebooks/spatialdata/dactyloscopy.ipynb` | one well | 2160x2160 fields; outlines carry a singleton axis |
| `cpg0037-oasis` | `notebooks/anndata/oasis.ipynb` | 52 MB | `axiom` source; 299 duplicated plate/well rows |
| `cpg0037-oasis` | `notebooks/spatialdata/oasis.ipynb` | one well | Stage coordinates but no pixel size, so fields only |
| `cpg0038-tegtmeyer-neuropainting` | `notebooks/anndata/neuropainting.ipynb` | 4.5 MB | Astrocytes and neurons, 20x and 63x |
| `cpg0038-tegtmeyer-neuropainting` | `notebooks/spatialdata/neuropainting.ipynb` | one well | z stacks, so `plane` picks one |
| `cpg0039-garcia-fossa-livecellpainting` | `notebooks/anndata/garcia_fossa_live.ipynb` | 6.1 MB | Live Cell Painting |
| `cpg0040-garcia-fossa-AgNP` | `notebooks/anndata/garcia_fossa_agnp.ipynb` | 288 KB | Nanoparticle size, dose, time |
| `cpg0040-garcia-fossa-AgNP` | `notebooks/spatialdata/garcia_fossa_agnp.ipynb` | one well | Two channels named without `Orig`; colour outlines |
| `cpg0047-amish` | `notebooks/anndata/amish.ipynb` | 1.2 MB | Cell line, density, timepoint |
| `cpg0047-amish` | `notebooks/spatialdata/amish.ipynb` | one well | Images and segmentations |

[EUbOPEN](https://zenodo.org/records/10894238) is covered by `notebooks/anndata/eubopen.ipynb`; it is on Zenodo, not the gallery.

## Not covered, and why

Each of these was checked by listing the accession's `workspace/` prefixes, not assumed.

| Accession | What is published | Why not |
| --- | --- | --- |
| `cpg0003-rosetta` | `preprocessed_data/` | Re-processed copies of CDRP, LINCS-Pilot1, LUAD and TA-ORF, all already covered here from their own accessions |
| `cpg0015-heterogeneity` | `supplementary/` | No profile tables published |
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
| `idr0133-dahlin-cellpainting` | 403 | 372 | `notebooks/anndata/idr0133.ipynb` |
| `idr0016-wawer-bioactivecompoundprofiling` | 24 | 0 | annotation only; features are `cpg0012` |
| `idr0033-rohban-pathways` | 77 | 0 | annotation only; features are `cpg0017` |
| `idr0088-cox-phenomicprofiling` | 21 | 0 | annotation only |
| `idr0080-way-perturbation` | 17 | 0 | annotation only |
| `idr0160-lippincott-pyroptosis` | 38 | 0 | annotation only |
| `idr0036-gustafsdottir-cellpainting` | — | — | no annotation repository on GitHub; publishes no profiles here or in the gallery (`cpg0030`) |

`idr0035-caie-drugresponse` is BBBC021 rather than Cell Painting proper and is covered from `cpg0010`.

`idr0093-mueller-perturbation` publishes ~380 per-well features, but under its own naming (`nuclei_area_mean`, `frac_G1`) rather than CellProfiler's, and it is an EU/nascent-RNA assay rather than Cell Painting, so it is out of scope.
Detection here keys on CellProfiler naming, so features published under another convention would be missed.
