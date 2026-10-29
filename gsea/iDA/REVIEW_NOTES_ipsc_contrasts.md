# GSEA review notes: iPSC contrasts (2026-10-03)

## Input and scope

- **Dataset:** the Welch t ranking was computed on `adata_strict_no_contaminants`: 2467 protein groups, complete in every sample, no imputation (`params.json`: `n_protein_rows` 2467, `pct_missing` 0.0). The imputed matrix (6069 protein groups) was not used. No cluster in this review is driven by imputed values. The imputation check in the review prompt therefore does not apply to this run.
- **Review status:** this is a fresh run. All 2,270 terms in the six decision files were `NEW` and none was `MOVED`. Every HIGH cluster was reviewed from its `core_genes`, or from its `consensus_genes` where `core_size` < 5.
- **System:** iPSC-derived dopaminergic neurons, day 56, seeded with PFF or PBS. Terms that require immune cells, vasculature, infection, muscle, embryonic germ layers or epithelia were excluded when their core was generic or incoherent. When the core was a recognisable neuronal module, the cluster was relabelled and kept.
- **Contrasts and figure columns:** Soma PFF/PBS, Neur PFF/PBS, SAgg/Soma, NAgg/Neur, NAgg/SAgg, PBS N/S. These form the six columns of `20261003_ipsc_contrasts_dotplot` (k = 6). The first two also form `20261003_pff_contrasts_dotplot` (k = 2).
- **Figure parameters:** `cell_terms="purity"`, `purity_floor` 0.45, `coverage_floor` 0.10, `cell_purity_floor` 0.30, `anchor_rule="coverage"` and `min_rows_per_contrast` 3. The six-contrast figure uses 25 up and 25 down rows. The PFF figure uses 20 + 20 rows, with `PFF_FIGSIZE` and `PFF_DOT_AREA` set as in the notebook.

## Review coverage (HIGH clusters at FDR < 0.05)

| contrast | up labelled / significant | down labelled / significant |
|---|---|---|
| Soma PFF/PBS | 11 / 12 | 11 / 20 |
| Neur PFF/PBS | 25 / 40 | 32 / 40 |
| SAgg/Soma | 30 / 40 | 29 / 40 |
| NAgg/Neur | 16 / 39 | 21 / 22 |
| NAgg/SAgg | 21 / 40 | 22 / 30 |
| PBS N/S | 24 / 40 | 30 / 40 |

Every block reaches the floor of 10 labelled clusters. The unlabelled clusters fall into three groups: excluded clusters (listed below); clusters whose core is a mixture of modules that already carry a label in the same contrast (rep blanked, no row); and clusters whose core is two modules linked through a chain.

The pass 2 replay assertion passed for all six files (22, 57, 59, 37, 43 and 54 representatives).

## Rows drawn: six-contrast figure (NES; * = FDR < 0.05; one gene set per row)

| Row | term used (all columns) | SomaPFF | NeurPFF | SAgg | NAgg | NAggSAgg | PBS_NS |
|---|---|---|---|---|---|---|---|
| Heterochromatin maintenance (DNMT/TRIM28) | GOBP_HETEROCHROMATIN_FORMATION | -1.68* | -1.66* | -1.79* | -2.03* | -1.54 | -1.97* |
| Spliceosome (U1/U2 snRNP) | REACTOME_MRNA_SPLICING_MINOR_PATHWAY | -0.49 | -1.85* | -1.92* | -1.95* | -1.69* | -1.97* |
| mRNA splicing regulation (SR/hnRNP) | GOBP_POSITIVE_REGULATION_OF_RNA_SPLICING | -0.90 | -1.80* | -1.94* | -1.40 | -1.71* | -1.93* |
| Cytoplasmic translation (ribosome) | REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION | +1.71* | -2.94* | -1.94* | -2.31* | -2.28* | -1.61* |
| Nucleosome assembly (FACT/SWI-SNF) | GOBP_PROTEIN_DNA_COMPLEX_ORGANIZATION | -1.39 | -1.53 | -1.66* | -1.43 | -1.37 | -1.82* |
| Histone deacetylation (NuRD) | REACTOME_CHROMATIN_ORGANIZATION | -1.79* | -1.66* | -1.63* | -1.27 | -0.99 | -1.71* |
| Chaperonin CCT/TRiC | GOBP_REGULATION_OF_DNA_BIOSYNTHETIC_PROCESS | -1.18 | -1.84* | -1.74* | -1.30 | -1.38 | -1.56* |
| Ribosome biogenesis | GOBP_RIBOSOMAL_LARGE_SUBUNIT_BIOGENESIS | -1.00 | -2.21* | -1.64* | -0.84 | -1.47 | -1.34 |
| DNA damage response (NHEJ) | GOBP_RESPONSE_TO_IONIZING_RADIATION | -0.98 | -1.81* | -1.58* | -1.50 | -0.71 | -1.18 |
| Nuclear pore and mRNA export | REACTOME_NUCLEAR_PORE_COMPLEX_NPC_DISASSEMBLY | +0.79 | -1.89* | -1.76* | -1.46 | -1.52 | -1.52* |
| mRNA stability and translational control | GOBP_POSITIVE_REGULATION_OF_MRNA_CATABOLIC_PROCESS | +1.51 | -1.65* | -1.66* | -1.95* | -1.82* | -1.72* |
| Translation initiation (eIF3/eIF4) | GOBP_CYTOPLASMIC_TRANSLATIONAL_INITIATION | +1.57 | -1.43 | -1.79* | -1.79* | -1.91* | -1.90* |
| Nuclear import (importins/RAN) | GOBP_IMPORT_INTO_NUCLEUS | -1.70* | -2.04* | -1.46 | +0.66 | -1.29 | -1.24 |
| tRNA aminoacylation | REACTOME_CYTOSOLIC_TRNA_AMINOACYLATION | +1.86* | -1.78* | -1.72* | -1.60 | -2.17* | -1.50* |
| Proteasome / UPS | REACTOME_TRANSCRIPTIONAL_REGULATION_BY_RUNX3 | -2.79* | -1.17 | -1.21 | +2.49* | -1.39 | -2.33* |
| mRNA deadenylation and decapping (CCR4-NOT/DCP) | REACTOME_DEADENYLATION_DEPENDENT_MRNA_DECAY | +1.34 | -1.24 | -1.58* | -1.74* | -1.49 | -1.61* |
| Chaperone-assisted UPS (VCP/UBQLN) | REACTOME_DNA_REPAIR | -0.78 | -1.69* | -1.61* | +1.19 | -1.10 | -1.68* |
| Nuclear envelope (lamins/LEM) | GOBP_NUCLEAR_ENVELOPE_ORGANIZATION | -1.34 | -1.08 | -1.27 | +1.56 | -1.07 | -1.65* |
| COPI retrograde transport | GOBP_RETROGRADE_VESICLE_MEDIATED_TRANSPORT_GOLGI_TO_ENDOPLASMIC_RETICULUM | +1.98* | -0.85 | -0.88 | -1.59 | -1.84* | -1.03 |
| ER to Golgi transport (COPII) | REACTOME_COPII_MEDIATED_VESICLE_TRANSPORT | +2.05* | +0.75 | -1.33 | -1.99* | -1.87* | -1.52* |
| Lysosomal glycosidases (sphingolipid catabolism) | GOBP_SPHINGOLIPID_METABOLIC_PROCESS | +1.54 | -0.72 | -1.15 | -0.85 | -1.59* | -0.97 |
| mRNA 3'-end processing | GOBP_PROTEIN_TETRAMERIZATION | +0.87 | +1.15 | -1.19 | -0.89 | -0.85 | -1.71* |
| Unfolded protein response (ER) | REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS | +1.07 | +2.02* | -1.01 | -1.36 | -1.88* | -1.33 |
| Cholesterol biosynthesis | GOBP_STEROL_BIOSYNTHETIC_PROCESS | -2.43* | -1.75* | +1.61* | +1.30 | -0.43 | +0.78 |
| Mitochondrial translation (mitoribosome) | REACTOME_MITOCHONDRIAL_TRANSLATION | +1.74* | -1.02 | -1.28 | -1.79* | +1.00 | +1.19 |
| Aldehyde and xenobiotic oxidoreductases | REACTOME_BIOLOGICAL_OXIDATIONS | -2.15* | -1.68* | +0.90 | +1.09 | +0.94 | +0.80 |
| Oxidative phosphorylation (ETC) | REACTOME_COMPLEX_I_BIOGENESIS | +2.27* | +2.59* | -1.25 | -1.75* | -0.49 | -1.07 |
| Integrin-ECM adhesion (ILK/FERMT2, laminins) | GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY | -2.04* | -1.67* | +2.14* | -0.65 | +1.65* | +2.39* |
| ESCRT and autophagosome maturation | GOBP_ENDOSOME_TO_LYSOSOME_TRANSPORT | +0.89 | +1.21 | -0.61 | +2.09* | -1.08 | -0.67 |
| Lysosomal transport and sorting | GOBP_LYSOSOMAL_TRANSPORT | +1.53 | +1.40 | -0.64 | +1.97* | -1.60* | -0.72 |
| Amino-acid transporters (SLC7A5/SLC3A2) | GOBP_ORGANIC_ACID_TRANSMEMBRANE_TRANSPORT | -1.42 | -1.50 | +1.92* | -1.27 | +2.20* | +2.07* |
| v-ATPase / lysosomal acidification | REACTOME_IRON_UPTAKE_AND_TRANSPORT | -1.25 | +1.34 | +2.41* | +1.94* | -1.01 | -1.13 |
| Purine nucleotide biosynthesis | GOBP_RIBONUCLEOSIDE_MONOPHOSPHATE_METABOLIC_PROCESS | -1.52 | -1.87* | +1.51 | +1.26 | +0.93 | +2.26* |
| Rab GTPase trafficking | REACTOME_RAB_GERANYLGERANYLATION | -1.76* | -0.86 | +2.52* | +1.49 | +1.13 | +2.31* |
| Ion transport (Na+/K+ ATPase, channels) | GOBP_POTASSIUM_ION_TRANSPORT | -1.33 | -0.49 | +2.29* | +1.49 | +1.22 | +2.06* |
| Glycolysis | GOBP_NUCLEOSIDE_DIPHOSPHATE_CATABOLIC_PROCESS | -1.21 | +0.79 | +1.70* | +2.00* | +1.08 | +0.93 |
| Synaptic vesicle endocytosis (SYNJ1/auxilin) | REACTOME_CARGO_RECOGNITION_FOR_CLATHRIN_MEDIATED_ENDOCYTOSIS | +0.55 | +1.72* | +1.73* | +2.06* | +1.03 | -0.76 |
| Heterotrimeric G-protein signalling | GOBP_ADENYLATE_CYCLASE_MODULATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY | +1.24 | +1.41 | +1.62* | -1.06 | +1.39 | +1.90* |
| Cytoskeletal adhesion dynamics (CLASP/ARF6) | GOBP_REGULATION_OF_EPITHELIAL_CELL_MIGRATION | +0.58 | -0.76 | +2.32* | +1.66 | +1.70* | +2.02* |
| Parkinson-associated proteins (PARK7/UCHL1/SNCA) | GOBP_ADULT_LOCOMOTORY_BEHAVIOR | +1.56 | +2.01* | +1.46 | +0.85 | +1.27 | +0.68 |
| Synaptic scaffolding (neurexin-CASK) | REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE | -0.92 | +1.80* | +2.48* | +0.73 | +2.02* | +1.98* |
| PKA and phosphodiesterase signalling | REACTOME_REGULATION_OF_INSULIN_SECRETION | +0.93 | +1.26 | +2.24* | -1.25 | +2.28* | +2.68* |
| Microtubule motors and tubulin | REACTOME_RECRUITMENT_OF_NUMA_TO_MITOTIC_CENTROSOMES | +0.70 | +1.95* | +2.10* | +1.32 | +1.49 | +1.24 |
| Semaphorin-plexin signalling | REACTOME_SEMAPHORIN_INTERACTIONS | -0.71 | +1.05 | +2.79* | +1.48 | +2.12* | +2.76* |
| MAPK/ERK signalling | REACTOME_SIGNALING_BY_MODERATE_KINASE_ACTIVITY_BRAF_MUTANTS | +1.16 | +1.37 | +1.87* | +1.36 | +1.80* | +1.98* |
| Calcium import and calcineurin | GOBP_CALCIUM_ION_TRANSMEMBRANE_TRANSPORT | +1.54 | +1.76* | +1.47 | +1.18 | +2.06* | +1.58* |
| Synaptic vesicle exocytosis (SNARE/SYT) | GOBP_EXOCYTIC_PROCESS | +0.98 | +2.04* | +2.82* | +1.31 | +2.06* | +2.74* |
| Axon development and neuron morphogenesis | GOBP_NEURON_PROJECTION_GUIDANCE | +1.08 | +2.01* | +2.92* | +1.26 | +2.43* | +3.04* |

The PFF two-contrast figure draws 31 rows. 29 of them are also rows of the six-contrast figure. For 23 of those 29, the PFF figure reads the same gene set. For 6 it chose a different one (`Proteasome / UPS`, `mRNA splicing regulation`, `Spliceosome`, `mRNA stability and translational control`, `Translation initiation`, `PKA and phosphodiesterase signalling`), because the term is chosen over the contrasts passed in. Two PFF rows are absent from the six-contrast figure: `Actomyosin organisation and cell polarity` and `Synapse assembly (neurexin-neuroligin)`. The two figures must therefore not be presented as the same measurement. The full per-cell record is in `ipsc_contrasts_row_table.csv` and in `20261003_pff_contrasts_dotplot_grid.csv`.

## Labels with no term clearing the purity floor in any contrast

These modules carry labelled clusters, but no tested gene set has a leading edge that is ≥ 45 % the label's genes. They are findings at protein level that do not appear as rows.

- Six-contrast figure: `Cell-cell adhesion (cadherins)` (SAgg/Soma, PBS N/S), `DNA replication and repair` (Neur PFF/PBS), `Eph-ephrin signalling (SRC/PAK)` (PBS N/S), `Nonsense-mediated mRNA decay` (Neur PFF/PBS), `Transcription elongation (Mediator)` (SAgg/Soma, PBS N/S), `ncRNA-mediated gene silencing` (Neur PFF/PBS, SAgg/Soma, NAgg/Neur), `rRNA synthesis and processing` (SAgg/Soma).
- PFF figure: `DNA replication and repair`, `Nonsense-mediated mRNA decay`, `Spectrin membrane skeleton (SPTAN1/SPTBN1)`, `mRNA 3'-end processing`, `ncRNA-mediated gene silencing`.
- The labels `Synapse assembly`, `Actin nucleation`, `Actomyosin organisation`, `Dense-core vesicle cargo (granins)`, `Spectrin membrane skeleton` and `Netrin/Robo guidance` are not rows of the six-contrast figure: they are absent from the label set that has a qualifying term and significance in at least one contrast after the per-direction trim. `Actomyosin organisation` and `Synapse assembly` are rows of the PFF figure. `Actin nucleation` was a PFF row only while `REACTOME_PARASITE_INFECTION`, an artifact-listed term, could fill it; after the artifact fix it is no longer drawn in either figure.

## Clusters excluded, and why

Wrong lineage. The term name requires a cell type or tissue that the culture lacks, and the core is generic:
- Soma PFF/PBS: `GOBP_CELLULAR_RESPONSE_TO_TYPE_II_INTERFERON` (GSN, ACTR2/3, FLNB, ZYX), `GOBP_ERYTHROCYTE_HOMEOSTASIS` (5 genes), `GOBP_SMOOTH_MUSCLE_CELL_PROLIFERATION`, `REACTOME_NEUTROPHIL_DEGRANULATION` (leading edge 68 genes; the proteasome part of it already has a row).
- Neur PFF/PBS: `GOBP_REGULATION_OF_MUSCLE_SYSTEM_PROCESS`, `GOBP_ERYTHROCYTE_HOMEOSTASIS`, `GOBP_BLASTOCYST_DEVELOPMENT`.
- SAgg/Soma: `GOBP_VIRAL_GENOME_REPLICATION`, `REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION`, `GOBP_DEFENSE_RESPONSE_TO_VIRUS`, `GOBP_NEGATIVE_REGULATION_OF_RESPONSE_TO_BIOTIC_STIMULUS`.
- NAgg/Neur: `REACTOME_NEUTROPHIL_DEGRANULATION` (contains S100A8/S100A9; see the open questions), `GOBP_BIOLOGICAL_PROCESS_INVOLVED_IN_SYMBIOTIC_INTERACTION`, `GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY`, `GOBP_NEGATIVE_REGULATION_OF_CYTOKINE_PRODUCTION`, `REACTOME_MUSCLE_CONTRACTION`, `REACTOME_BACTERIAL_INFECTION_PATHWAYS`, `GOBP_NEUROMUSCULAR_PROCESS`, `REACTOME_MYD88_INDEPENDENT_TLR4_CASCADE`, `GOBP_NEGATIVE_REGULATION_OF_RESPONSE_TO_EXTERNAL_STIMULUS`.
- NAgg/SAgg: `GOBP_PLATELET_ACTIVATION`, `GOBP_MORPHOGENESIS_OF_AN_EPITHELIUM`, `GOBP_POSITIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY`, `GOBP_RESPONSE_TO_VIRUS`, `GOBP_VIRAL_GENOME_REPLICATION`.
- PBS N/S: `GOBP_ENDOTHELIUM_DEVELOPMENT`, `GOBP_MORPHOGENESIS_OF_AN_EPITHELIUM`, `GOBP_EPITHELIAL_TUBE_MORPHOGENESIS`, `GOBP_DEFENSE_RESPONSE_TO_VIRUS`.

Incoherent core, with no single module among the core genes:
- Soma PFF/PBS: `GOBP_DENDRITIC_SPINE_DEVELOPMENT` (RNA granule plus synaptic genes), `GOBP_POSITIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION`, `REACTOME_CHROMATIN_ORGANIZATION`, `GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION`, `REACTOME_METABOLISM_OF_WATER_SOLUBLE_VITAMINS_AND_COFACTORS`, `GOBP_RESPONSE_TO_NUTRIENT`.
- Neur PFF/PBS: `GOBP_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY`, `GOBP_ORGANIC_HYDROXY_COMPOUND_TRANSPORT`, `GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_POTENTIAL`, `GOBP_PEPTIDE_TRANSPORT`, `GOBP_NERVOUS_SYSTEM_PROCESS`, `GOBP_REGULATION_OF_MEMBRANE_POTENTIAL`, `GOBP_OSSIFICATION`.
- SAgg/Soma: `GOBP_RESPONSE_TO_MECHANICAL_STIMULUS`, `GOBP_REGULATION_OF_CIRCADIAN_RHYTHM` (paraspeckle NONO/SFPQ/PSPC1 plus unrelated genes), `GOBP_POSITIVE_REGULATION_OF_TRANSCRIPTION_BY_RNA_POLYMERASE_II`.
- PBS N/S: `GOBP_POSITIVE_REGULATION_OF_TRANSCRIPTION_BY_RNA_POLYMERASE_II`.

Chained clusters (core ≤ 2 genes): NAgg/Neur `GOBP_REGULATION_OF_RECEPTOR_MEDIATED_ENDOCYTOSIS` and `REACTOME_CARGO_RECOGNITION_FOR_CLATHRIN_MEDIATED_ENDOCYTOSIS`; NAgg/SAgg `GOBP_CELL_SUBSTRATE_JUNCTION_ORGANIZATION` and `GOBP_EXOCYTOSIS` (core 1).

## Clusters relabelled because the database name was misleading

Each line gives the representative term, the core genes that decided the label, and the label assigned.

- `REACTOME_NUCLEAR_EVENTS_MEDIATED_BY_NFE2L2` (Soma PFF/PBS), `REACTOME_KEAP1_NFE2L2_PATHWAY` and `REACTOME_PCP_CE_PATHWAY` (NAgg/Neur), `REACTOME_PTEN_REGULATION` (NAgg/SAgg), and `REACTOME_RUNX1_REGULATES...`, `REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT` and `REACTOME_UB_SPECIFIC_PROCESSING_PROTEASES` (PBS N/S): core PSMA/PSMB/PSMC/PSMD subunits. Label: **Proteasome / UPS**. The Soma PFF/PBS leading edge also holds PRDX1, PGD, G6PD, IDH1, ME1 and TXNRD1, which are NRF2 target genes. That cell nevertheless reads as proteasome down (-2.79), not NRF2 activation.
- `GOBP_PLATELET_ACTIVATION` (Soma PFF/PBS): FLNA, VCL, ILK, ACTN1, CSRP1, MYH9. Label: **Integrin-ECM adhesion**.
- `REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION` (Neur PFF/PBS): DYNC1H1, DYNC1I1/I2, DYNC1LI2 and TUBB isoforms. Label: **Microtubule motors and tubulin**. The same applies to `REACTOME_RECRUITMENT_OF_NUMA_TO_MITOTIC_CENTROSOMES` (dynein and dynactin in post-mitotic neurons), `REACTOME_L1CAM_INTERACTIONS`, `REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4...` and `REACTOME_FACTORS_INVOLVED_IN_MEGAKARYOCYTE...`, whose cores are tubulin isoforms only.
- `REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT...` (Neur PFF/PBS): SCG2, SCG3, CHGB, VGF, APLP2, NUCB1. Label: **Dense-core vesicle cargo (granins)**, which is new. The same term in NAgg/SAgg has the core PDIA6, MIA3, CALU, NUCB1, KTN1, PRKCSH, CKAP4, which are ER luminal and export proteins with no granins. It was labelled **Unfolded protein response (ER)** (see the open questions).
- `REACTOME_SIGNALING_BY_VEGF`: in Neur PFF/PBS and NAgg/SAgg the core is the WAVE complex (ABI1/2, CYFIP2, BRK1, NCKAP1, BAIAP2), labelled **Actin nucleation (ARP2/3 complex)**. In SAgg/Soma it is SRC/FYN/RAC1/CDC42/RHOA/ABI/NRP, labelled **Netrin/Robo guidance and actin**. In PBS N/S it is RHOA/ROCK2/PAK1-3/CDC42/NCK2/SRC, labelled **Eph-ephrin signalling (SRC/PAK)**.
- `REACTOME_LEISHMANIA_INFECTION` (SAgg/Soma): ACTR2/3, ARPC2/3/4, ABI2. Label: **Actin nucleation (ARP2/3 complex)**.
- `REACTOME_EPH_EPHRIN_SIGNALING` (SAgg/Soma): MYH9/10, MYL6, ROCK2, RHOA, CFL1, PAK2. Label: **Actomyosin organisation and cell polarity**.
- `REACTOME_ION_CHANNEL_TRANSPORT`, `GOBP_DEVELOPMENTAL_MATURATION` (SAgg/Soma) and `GOBP_MACROAUTOPHAGY` (NAgg/Neur): ATP6V0/ATP6V1 subunits. Label: **v-ATPase / lysosomal acidification**.
- `GOBP_CIRCULATORY_SYSTEM_PROCESS` (SAgg/Soma) and `GOBP_CARDIAC_CONDUCTION` (NAgg/SAgg): ATP1A1/A3, ATP1B1, ATP2B1/4, SLC8A2. Label: **Ion transport (Na+/K+ ATPase, channels)**.
- `REACTOME_REGULATION_OF_INSULIN_SECRETION` (NAgg/SAgg, PBS N/S) and `REACTOME_SIGNALING_BY_GPCR` (PBS N/S): GNAI1/2, GNAQ, GNA11, GNAS, GNB1/2, GNG2. Label: **Heterotrimeric G-protein signalling**.
- `REACTOME_SIGNALING_BY_NTRKS` (PBS N/S): KRAS, HRAS, NRAS, MAP2K2, RAP1A, SRC. Label: **MAPK/ERK signalling**.
- `REACTOME_SUMOYLATION`, `REACTOME_TRANSCRIPTIONAL_REGULATION_BY_SMALL_RNAS`, `REACTOME_SLC_TRANSPORTER_DISORDERS` and `REACTOME_SARS_COV_2_ACTIVATES_MODULATES_INNATE...`: NUP core (NUP93/98/107/133/155/160/205/214, TPR, RANBP2). Label: **Nuclear pore and mRNA export**.
- `GOBP_DNA_BIOSYNTHETIC_PROCESS`, `GOBP_TELOMERE_ORGANIZATION`, `GOBP_TELOMERE_MAINTENANCE_VIA_TELOMERE_LENGTHENING` and `GOBP_POSITIVE_REGULATION_OF_CHROMOSOME_ORGANIZATION`: TCP1, CCT2-8. Label: **Chaperonin CCT/TRiC**.
- `GOBP_VIRAL_GENE_EXPRESSION` (Neur PFF/PBS, PBS N/S): EIF3A/B/D/F/L. Label: **Translation initiation (eIF3/eIF4)**.
- `REACTOME_SELENOAMINO_ACID_METABOLISM` (NAgg/SAgg, PBS N/S) and `REACTOME_TRANSCRIPTIONAL_AND_POST_TRANSLATIONAL_REGULATION_OF_MITF_M...` (Soma PFF/PBS): RARS1, MARS1, KARS1, IARS1, EPRS1, LARS1, QARS1, DARS1, AIMP1/2. Label: **tRNA aminoacylation**.
- `GOBP_NUCLEOSIDE_DIPHOSPHATE_CATABOLIC_PROCESS` (NAgg/Neur): PKM, GAPDH, ALDOC, ENO1, LDHA, GPI, PFKP, PGK1. Label: **Glycolysis**, which is new.
- `GOBP_REGULATION_OF_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY` (NAgg/Neur): UFD1, NPLOC4, HDAC6, HSPA8, HSPA1A, LAMP2. Label: **Chaperone-assisted UPS (VCP/UBQLN)**.
- `GOBP_MITOTIC_CYTOKINESIS` (NAgg/Neur): CHMP2A, CHMP4B, CHMP1A, VPS4A, IST1. Label: **ESCRT and autophagosome maturation**.
- `GOBP_LOCOMOTORY_BEHAVIOR` (Soma PFF/PBS, Neur PFF/PBS): SNCA, SNCG, UCHL1, PARK7 (and HTRA2). Label: **Parkinson-associated proteins**.
- `REACTOME_OPIOID_SIGNALLING` and `REACTOME_G_ALPHA_I_SIGNALLING_EVENTS`: PRKACA, PRKACB, PRKAR1A, PRKAR2A/B, CALM1. Label: **PKA and phosphodiesterase signalling**.

## New labels and the gene lists written for them

Each list was written from the label name before any leading edge was read, and was added to `src/gsea_label_genes.yaml`.

- **Mitochondrial translation (mitoribosome)**: every MRPL and MRPS subunit, plus DAP3, PTCD3, AURKAIP1, CHCHD1, GADD45GIP1, TUFM, TSFM, GFM1, GFM2, MTIF2, MTIF3, MRRF, MTRF1L, ERAL1, MTG1, MTG2, MALSU1, TACO1. Used for Soma PFF/PBS (up) and NAgg/Neur (down).
- **Dense-core vesicle cargo (granins)**: CHGA, CHGB, SCG2, SCG3, SCG5, VGF, PCSK1N, CPE, PCSK1, PCSK2, PAM, PTPRN, PTPRN2, CADPS, CADPS2, BAIAP3, NPY, TAC1, PENK, CARTPT, SST, CALCA. Used for Neur PFF/PBS (up). A qualifying term exists, but the row is not drawn in either figure because of the row budget.
- **Glycolysis**: HK1, HK2, GPI, PFKL, PFKM, PFKP, ALDOA, ALDOC, TPI1, GAPDH, PGK1, PGAM1, ENO1, ENO2, ENO3, PKM, LDHA, LDHB, BPGM, SLC2A1, SLC2A3. Used for NAgg/Neur (up).

All other labels reuse existing entries. `label_vocabulary.txt` did not exist and was written by `update_vocabulary`; it contains the 65 labels used.

## Labels with near-identical gene lists

`Stress granule assembly` shares 15 of its 31 genes with `mRNA stability and translational control` (ELAVL1, PABPC1, FMR1, FXR1/2, DDX6, CAPRIN1, CIRBP, CSDE1, PUM1, NUFIP2, ...). The two clusters I had first given the stress-granule label (NAgg/SAgg and PBS N/S, both `REGULATORY_NCRNA_MEDIATED_GENE_SILENCING` with cores of TIAL1/STAU1/ELAVL1/FMR1/PABPC1) were moved to `mRNA stability and translational control`. `Stress granule assembly` is therefore not used.

`Clathrin-mediated endocytosis` and `Synaptic vesicle endocytosis (SYNJ1/auxilin)` share clathrin, AP2, dynamin, EPS15/L1, ITSN, SNAP91, PICALM, FCHO and NECAP. Only the synaptic-vesicle label was used, because every endocytic cluster carried SH3GL2, DNAJC6, SYNJ1, AMPH or DNM1/3. This keeps the module as one row.

All other pairs of labels in use overlap by less than 0.30 (overlap coefficient on the lists).

## Redundancy (rows whose drawn terms share ≥ 50 % of a leading edge)

All pairs below were kept apart. In each pair the two rows differ in sign or significance in at least one column.

- `Nuclear import (importins/RAN)` vs `Nuclear pore and mRNA export` (4/6 columns). Soma PFF/PBS: import -1.70*, pore +0.79. NAgg/Neur: import +0.66, pore -1.46. PBS N/S: pore -1.52*, import -1.24.
- `Cytoplasmic translation (ribosome)` vs `Ribosome biogenesis` (2/6). Soma PFF/PBS: +1.71* vs -1.00. NAgg/Neur: -2.31* vs -0.84. This is the worked example in the review prompt.
- `Heterotrimeric G-protein signalling` vs `PKA and phosphodiesterase signalling` (2/6). NAgg/SAgg: PKA +2.28* vs G-protein +1.39 (not significant). PBS N/S: +2.68* vs +1.90*. They share significance status in 4 of 6 columns. This pair is the weakest separation in the figure; see the open questions.
- `Chaperone-assisted UPS` vs `DNA damage response` (1/6) and `ESCRT` vs `Lysosomal transport and sorting` (1/6). These overlap in one column only, and the rows differ elsewhere: NAgg/SAgg -1.08 vs -1.60*.
- An earlier draw also showed `Endolysosomal trafficking (VPS/RAB)` and `ESCRT and autophagosome maturation` filled from the same term (`GOBP_ENDOSOME_TO_LYSOSOME_TRANSPORT`), so the two rows were identical. The Endolysosomal label came from a single Neur PFF/PBS cluster (`GOBP_MEMBRANE_FUSION`, core SYT1/7/11, STX7/12, SNAP29, EEA1, ANKFY1, CHMP1A), whose core is SNARE fusion plus endosomal tethers. That label was removed and the cluster left without a row.

## Off-list cells

The six-contrast figure has 51 hollow cells for which another significant term covers the module. **No gene list was changed in response.**

- **31 cells: the other term's purity is 0.27–0.44, below the 0.45 floor.** These are borderline enrichments that the label's list recognises only in part, mostly broad GO parents such as `REACTOME_CHROMATIN_ORGANIZATION`, `GOBP_DENDRITIC_SPINE_DEVELOPMENT` and `GOBP_CALCIUM_MEDIATED_SIGNALING`. Widening a list to admit any of them would break the rule that a list is never widened for a particular term. They were left.
- **20 cells: the other term clears the floor (purity 0.47–1.00).** Here the gene list is adequate. The cell is hollow because the row reads one gene set in every column, chosen where the module is enriched. Three of these cells should be stated in the legend or text:
  - `Synaptic vesicle endocytosis`, PBS N/S: the row term reads -0.76, while `GOBP_PRESYNAPTIC_ENDOCYTOSIS` reads +2.22 (FDR 4e-4, purity 0.47). The sign is opposite.
  - `ESCRT and autophagosome maturation`, NAgg/SAgg: the row term reads -1.08, while `GOBP_MEMBRANE_FISSION` reads +1.61 (FDR 0.035, purity 0.55). The sign is opposite.
  - `Proteasome / UPS`, NAgg/SAgg: the row term reads -1.39 (FDR 0.16), while `REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA`, whose leading edge is 100 % proteasome subunits, reads -1.58 (FDR 0.037). Before the artifact fix this cell was filled at -1.63 by the artifact-listed `REACTOME_CROSS_PRESENTATION_OF_SOLUBLE_EXOGENOUS_ANTIGENS_ENDOSOMES`. The proteasome-down signal in NAgg/SAgg is therefore significant on a pure proteasome term, but not on the term the row uses.
- The full list of 51 cells is in the run output (`report_offlist_cells`).

## Clusters I was unsure about, and what would resolve each

- **Aggrephagy in NAgg/Neur** (`REACTOME_AGGREPHAGY`, core UBB, HDAC6, DYNC1I2, TUBB4B, TUBA1A). No row was given. The core is 5 genes and is shared with the tubulin and UPS rows. A list for "Aggrephagy (HDAC6/dynein)" written from the name would decide whether this is a module separate from the microtubule row.
- **NAgg/SAgg `REACTOME_EPH_EPHRIN_SIGNALING`** (AP2A1/A2/M1/S1 together with FYN, SRC, HRAS, RAC1). This core is two modules. No row was given.
- **NAgg/SAgg IGF-transport term labelled UPR.** The core proteins are ER luminal and export proteins (PDIA6, CALU, NUCB1, PRKCSH, MIA3, KTN1, CKAP4); none is a UPR sensor or effector. The UPR row in the figure is filled by the same IGF term. If the reader would take "UPR" as stress signalling, an "ER lumen and secretory cargo" label with its own list would state the finding more precisely.
- **S100A8/S100A9 in the NAgg/Neur neutrophil-degranulation leading edge.** These are myeloid proteins and should not be present in a neuronal culture. Check whether they are carry-over from handling, not culture protein, before any interpretation of that contrast's up block.
- **`GOBP_MEMBRANE_FUSION` (Neur PFF/PBS)**: described above. Splitting the cluster at a higher overlap threshold would show whether it holds an endosomal module.

## Edits beyond labelling

- **Low-priority representatives blanked:** 211 pre-filled `rep` values in low-priority clusters were set to blank (Neur PFF/PBS 25, SAgg/Soma 88, NAgg/SAgg 38, PBS N/S 60). Left in place, they would have reached the figure under their database names; `plot_combined` would have raised a `KeyError`, because those names have no gene list.
- **HIGH representatives blanked:** for HIGH clusters that were not labelled and not excluded, the pre-filled `rep` was blanked, with no row for that cluster.
- **Representatives not moved:** where a cluster was kept, the label sits on the pre-filled representative.
- **Columns and code:** only `display_choice` and `display_label` were changed. The run asserted that all other columns are byte-identical. `src/gsea_label_genes.yaml` gained the three new lists. No code was changed.
- **Code changed after the review (2026-10-03), in `src/gsea_pipeline.py`:**
  1. **Artifact terms.** `plot_combined` now removes `cfg.artifact_pathways` from each raw table before the purity selection. `report_offlist_cells` excludes them from the "other set" search. Before this change the list was applied only when clusters were built (`build_review`), so an artifact term could fill a row. Effects on this figure: `Amino-acid transporters` changed from `REACTOME_CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL` to `GOBP_ORGANIC_ACID_TRANSMEMBRANE_TRANSPORT` (same significant columns: SAgg/Soma, NAgg/SAgg, PBS N/S). `Proteasome / UPS` changed from `REACTOME_CROSS_PRESENTATION...` to `REACTOME_TRANSCRIPTIONAL_REGULATION_BY_RUNX3`, and its NAgg/SAgg cell went from filled to hollow (see Off-list cells). The PFF figure lost `Actin nucleation`, which only `REACTOME_PARASITE_INFECTION` had filled. No artifact term fills any cell in either figure.
  2. **Exact ties.** The term for a row was chosen with `max()` over a set of term names. Set order changes between Python processes, so exact ties were broken differently on each run. `REACTOME_NUCLEAR_PORE_COMPLEX_NPC_DISASSEMBLY` and `REACTOME_REGULATION_OF_GLUCOKINASE_BY_GLUCOKINASE_REGULATORY_PROTEIN` detect the same 20 genes in every contrast, so the nuclear pore row named one or the other from run to run, with identical values. The candidates are now sorted first, so the alphabetically first of a tie wins. Two draws under different `PYTHONHASHSEED` values produced byte-identical grids.
  3. **Tests.** Three tests were added to `src/test_gsea_pipeline.py`. Each was checked to fail on the code without the fix (the tie test under 5 hash seeds). The suite (28 tests) passes. pytest is not installed in any permitted environment, so the suite was run with a small runner that calls the test functions directly.

## Methods sentences required (from GSEA_PIPELINE.md §9)

- FDR is controlled within each contrast. "Significant in ≥ 1 of k contrasts" (k = 6, or k = 2 for the PFF figure) is a display rule, not an inferential claim.
- Each row is one gene set, chosen where the module is enriched, with a leading edge that is ≥ 45 % the label's genes there and ≥ 30 % in every column drawn filled.
- Row admission uses a per-contrast quota of 3 rows per direction.
