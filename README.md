# cell-painting-io

Experimental! A collection of scripts &amp; utilities to get Cell Painting datasets into scverse data structures.

`src/cell_painting_io/` holds the dataset-agnostic parts; each notebook in `notebooks/` adds only what its dataset needs.
Everything that has differed between datasets is a parameter rather than an assumption — missing-value sentinels, metadata prefixes, channel naming, columns that disagree across plates, features that are blown up or constant.

[COVERAGE.md](COVERAGE.md) lists every Cell Painting Gallery and IDR accession with its status, and a reason for each one not covered.

Each notebook starts with the command that fetches its data and a path constant to point at the download.
All but one of them read profiles only and never touch the images.

| Notebook (in `notebooks/`) | Dataset | Accession | Wells x features |
| --- | --- | --- | --- |
| `wawer` | CDRP bioactives, 30k compounds | `cpg0012` | 153,022 x 781 |
| `lincs` | LINCS, 1,571 compounds at six doses | `cpg0004` | 52,223 x 1,603 |
| `jump_crispr` | JUMP CRISPR, assembled | `cpg0016-jump-assembled` | 51,185 x 599 |
| `eubopen` | EUbOPEN compounds | Zenodo 10894238 | 39,206 x 7,703 |
| `idr0133` | Dahlin bioactives, from IDR | `idr0133` | 23,108 x 372 |
| `oasis` | OASIS, axiom site | `cpg0037` | 22,042 x 155 |
| `periscope` | PERISCOPE pooled screen, genes not wells | `cpg0021` | 20,384 x 3,640 |
| `cpjump1` | JUMP pilot, batch `2020_11_04_CPJUMP1` | `cpg0000` | 19,498 x 903 |
| `jump_adipocyte` | Adipocytes | `cpg0014` | 9,216 x 848 |
| `varchamp` | Coding variants | `cpg0020` | 7,177 x 853 |
| `dactyloscopy` | Dactyloscopy, five cell lines | `cpg0025` | 6,420 x 2,121 |
| `luad` | LUAD alleles, wild-type and mutant | `cpg0031` | 6,144 x 564 |
| `cellpainting_protocol` | Protocol variants | `cpg0001` | 5,794 x 187 |
| `rare_diseases` | Rare disease variants | `cpg0026` | 5,634 x 1,795 |
| `gerry` | Gerry bioactivity, no annotation | `cpg0005` | 4,608 x 1,467 |
| `cmqtl` | iPSC lines | `cpg0022` | 4,148 x 1,327 |
| `pki` | Kinase inhibitors, dose and MoA | `cpg0008` | 3,072 x 1,061 |
| `oasis_pilot` | OASIS pilot, U2OS and HepaRG | `cpg0033` | 2,303 x 718 |
| `molglue` | Molecular glues | `cpg0009` | 1,920 x 7,587 |
| `kelley` | Bortezomib resistance | `cpg0028` | 1,920 x 221 |
| `rohban` | TA-ORF, 323 genes overexpressed | `cpg0017` | 1,918 x 847 |
| `jump_scope` | JUMP-SCOPE, microscopes and settings | `cpg0002` | 1,529 x 819 |
| `miami` | MIAMI | `cpg0006` | 1,387 x 264 |
| `bortezomib` | Bortezomib clones, unnormalized parquet | `cpg0024` | 768 x 4,695 |
| `chroma` | Alternative dyes | `cpg0029` | 767 x 675 |
| `caie` | Caie drug response | `cpg0010` | 632 x 516 |
| `garcia_fossa_live` | Live Cell Painting | `cpg0039` | 536 x 82 |
| `amish` | Amish cohort, density and timepoint | `cpg0047` | 468 x 790 |
| `pooled_rare` | Pooled rare variants, barcodes | `cpg0032` | 290 x 4,046 |
| `lipocyte` | Lipocyte Profiler, rows are patient x cell type | `cpg0011` | 225 x 2,870 |
| `neuropainting` | Astrocytes and neurons, 20x and 63x | `cpg0038` | 188 x 563 |
| `garcia_fossa_agnp` | Silver nanoparticles, size and dose | `cpg0040` | 180 x 112 |

## Images

`cell_painting_io.read_plate` reads one plate of a Cell Painting Gallery source into a [SpatialData](https://spatialdata.scverse.org) object: fields of view as Images, the CellProfiler Nuclei, Cells and Cytoplasm segmentations as Labels, the wells as Shapes, and the well- and cell-level measurements as Tables.
Every element sits in a field, a well and a plate coordinate system, built from the stage coordinates the microscope recorded, and carries the plate barcode in its name so two plates concatenate.

The gallery publishes one-pixel object outlines rather than segmentation masks, so the Labels are reconstructed from them, keeping CellProfiler's object numbers so the per-object tables join onto them.

`notebooks/cpjump1_to_spatialdata.ipynb` runs it on two wells of two JUMP pilot plates, roughly 700 MB against the tens of terabytes of images for that accession.

## CellProfiler exports

`cell_painting_io.read_cellprofiler_export` reads one plate folder written by the `ExportForSpatialData` CellProfiler module, which exports images, segmentation masks and the per-cell table straight out of a pipeline run.
Nothing is reconstructed, unlike the gallery path: the module writes real label arrays rather than outlines, and a per-cell table already joined across compartments.
The reader builds Images, Labels and a `cells` Table from the manifest the module puts in the table's `uns`, so it never walks the folder or parses a file name, and a folder that was moved still reads.
Arrays the module recorded as failed are skipped, along with the rows that annotate them.
`cellprofiler_export_plates` lists the plate folders of one run, for a run that covered several plates.

The module does not export stage coordinates yet, so each element sits in its own field-of-view coordinate system and the fields are not placed relative to each other.
Element names are the exporter's own and do not yet match the `{plate}_{well}_s{site}_cells` names `read_plate` uses. Reconciling the two conventions is still open.

## Environment

```bash
mamba create -n cell-painting-io -c conda-forge python=3.14 anndata scanpy pandas numpy umap-learn leidenalg matplotlib jupyterlab
pip install -e '.[notebooks,spatial,test]'
prek install
```

The notebooks import `cell_painting_io`, so the package needs installing before they run.
