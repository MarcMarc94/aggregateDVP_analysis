# Deep Visual Proteomics of α-synuclein inclusions: data analysis

Analysis code for the manuscript "Deep Visual Proteomics reveals selective organization
of α-synuclein inclusions" (Oeller et al., in preparation). Five Jupyter notebooks read
the protein and phosphosite quantities of the aggregate Deep Visual Proteomics (DVP)
experiments and make the data panels of Figures 2 to 6 and of the supplementary figures
listed below.

## Notebooks

The numbers follow the order of the figures in the paper; each notebook runs on its own.

| Notebook | Samples | Figures |
|---|---|---|
| `1_HEK.ipynb` | HEK293 cells expressing A53T α-synuclein, seeded with α-synuclein preformed fibrils (PFF) or treated with PBS, without and with bafilomycin A1 | Figure 2, Figures S2 and S3 |
| `2_iDA.ipynb` | iPSC-derived dopaminergic (iDA) neurons treated with PFF or PBS | Figure 3, Figure S4 |
| `3_BrainTissue.ipynb` | post-mortem brain tissue of donors with Parkinson's disease or Lewy body disease and of controls, from two brain banks, and a replication cohort | Figure 4, Figure S6 |
| `4_CrossSystem.ipynb` | the HEK293, iDA and brain datasets, compared | Figure 5, Figures S7 and S8 |
| `5_Nanophos.ipynb` | phosphoproteome and proteome of HEK293 cells, from the same lysates | Figure 6, Figures S9 and S10 |

The first cell of every notebook lists its inputs, its outputs and the abbreviations it
uses. Every path is set in the first code cell of section 2. Each notebook runs from top
to bottom without the others; the cross-system notebook loads and filters the three
datasets itself.

## Repository layout

```
README.md  LICENSE  environment.yml  GSEA_REVIEW_PROMPT.md
*.ipynb            the five notebooks; they import src/ by a path relative to their folder
src/               modules, gene-set caches, UniProt-to-gene table
data/<dataset>/    sample annotations (metadata.csv) and annotation resources; the search
                   outputs from PRIDE go here
gsea/<dataset>/    GSEA decision files (review.csv) and review notes; the GSEA results are
                   written next to them
results/           result tables, written by the notebooks
figures/           figure panels and assembled figures, written by the notebooks
```

Two panels are images made outside the notebooks and are not in this repository:
Figure 3h (immunofluorescence, `2_iDA`) and Figure S3f (microscopy, `1_HEK`). The
figure-assembly cell of each notebook states how to assemble the figure without it.

## Installation

```
git clone https://github.com/MarcMarc94/aggregateDVP_analysis.git
cd aggregateDVP_analysis
conda env create -f environment.yml
conda activate aggregateDVP_analysis
python -m ipykernel install --user --name aggregateDVP_analysis
```

`environment.yml` pins the package versions that produced the published results, with
alphapepttools at its released version 0.2.0. The figure assembly converts SVG files with
CairoSVG, which calls the cairo library; conda installs cairo from conda-forge.

## Data

The raw files and the search outputs are deposited in the ProteomeXchange Consortium via
the PRIDE partner repository; the dataset identifier will be added upon publication. The
search outputs are DIA-NN reports for the HEK293, iDA and brain datasets and Spectronaut
reports for the phosphoproteomics dataset.

Copy the search outputs of each dataset from PRIDE into its folder in `data/`. The sample
annotations (`metadata.csv`) and the annotation resources are part of this repository.

| Folder | Files from PRIDE | In this repository | Read by |
|---|---|---|---|
| `data/HEK/DIANN_OptimizedSettings/` | `report.pg_matrix.tsv`, `report.pr_matrix.tsv` | `metadata.csv` | `1_HEK`, `4_CrossSystem` |
| `data/iDA/DIANN_OptimizedSettings/` | `report.pg_matrix.tsv`, `report.pr_matrix.tsv` | `metadata.csv` | `2_iDA`, `4_CrossSystem` |
| `data/BrainTissue/ADVP04_011_Combined_withBannerLibrary/` | `report.pg_matrix.tsv`, `report.pr_matrix.tsv` (discovery cohort) | `metadata.csv` | `3_BrainTissue`, `4_CrossSystem` |
| `data/BrainTissue/ADVP04_011_17_CombinedwithBannerLibrary/` | `report.pg_matrix.tsv` (discovery and replication cohorts, searched together) | `metadata.csv` | `3_BrainTissue` |
| `data/Nanophos/` | the Spectronaut phosphopeptide report (`…_Report_Denys_ms2_level.tsv`), proteome report (`…_proteome_aggregates_new_fasta_Report.tsv`) and fragment report (`…_fragments_Report.tsv`, 3.1 GB), and the search FASTA (`260522_Human_CAN_20416entries_ASYN_A53T.fasta`) | `metadata.csv` | `5_Nanophos` |
| `data/CrossSystem/` | — | `LITERATURE_MASTER.xlsx` (published gene lists); `Organelles/`: Gene Ontology `go-basic.obo` (release 2026-07-26), human annotations `goa_human.gaf.gz` (generated 2026-05-21), organelle module list; `Complexes/`: Complex Portal `9606.tsv` and two UniProt tables (downloaded 2026-09-10); `Physicochemical_Properties/`: UniProt sequences and CamSol scores | `4_CrossSystem` |

`4_CrossSystem` reads `report.pg_matrix.tsv` and `metadata.csv` of the HEK, iDA and brain
discovery folders. The sources of the annotation resources are given in the markdown of
the sections that read them.

## Running a notebook

The notebooks are stored without outputs. Open a notebook in Jupyter or VS Code, select
the kernel `aggregateDVP_analysis` and run all cells. The working directory must be the
repository root: the notebooks import the modules in `src/` by a path relative to it.

The notebooks write result tables to `results/`, GSEA results to `gsea/` and figure
panels to `figures/`. The name of a panel file starts with its panel letter and the date
of the run; the last section of each notebook assembles the figures from the newest file
of each letter.

## Gene sets and GSEA decisions

The gene set enrichment analyses (GSEA) use MSigDB gene sets cached in `src/`: Reactome
and GO biological process (`gsea_genesets_cache.parquet`, retrieved 2026-07-31) and GO
cellular component (`gsea_genesets_gocc_cache.parquet`, retrieved 2026-09-08). Both were
identical to MSigDB on 2026-10-29. A notebook that does not find a cache downloads the
current MSigDB release, which can differ from the cached one.

Each GSEA contrast has a decision file `review.csv` (representative term, exclusion and
display label of each cluster of significant terms). The decisions were made by the
authors with the help of an AI review (Claude, Anthropic) that followed the rules in
`GSEA_REVIEW_PROMPT.md`. The decision files are inputs of the notebooks; a re-run of the
GSEA keeps every decision by term name.

## Modules in `src/`

| File | Content |
|---|---|
| `figure_style.py` | page grid, fonts, line widths and colours of all panels; panel saving and figure assembly |
| `gsea_pipeline.py`, `GSEA_PIPELINE.md` | GSEA prerank, clustering of significant terms, decision file, dot plots; method description |
| `gsea_families.yaml`, `gsea_label_genes.yaml` | gene-family dictionary for the review; gene list of each row label of the iDA multi-contrast figure |
| `gsea_genesets_cache.parquet`, `gsea_genesets_gocc_cache.parquet` | MSigDB gene sets (see above) |
| `preprocessing.py` | contaminant flags |
| `heatmap.py`, `boxplot.py`, `scatter.py`, `dotplot.py` | plot functions |
| `crosssystem_plots.py` | plot functions of `4_CrossSystem` |
| `nanophos_preprocessing.py`, `PeptideCollapse_v4.py` | preprocessing of the phosphoproteomics data; collapse of Spectronaut precursors to peptides and sites |
| `Uniprot_to_GeneSymbol.csv` | UniProt-to-gene table from the alphapepttools documentation |

## Citation

Oeller M. et al. Deep Visual Proteomics reveals selective organization of α-synuclein
inclusions. In preparation.

## License

Apache License 2.0; see `LICENSE`.
