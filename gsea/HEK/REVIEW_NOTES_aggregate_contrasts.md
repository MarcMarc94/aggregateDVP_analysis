# GSEA Pass 2 review — aggregate contrasts

Reviewed 2026-09-24 by Claude, following `GSEA_REVIEW_PROMPT.md` and
`src/GSEA_PIPELINE.md` (§7 review, §9 multi-contrast).

Contrasts, all HEK293, Welch t ranking, gene-permutation null:

| folder | figure column |
|---|---|
| `PureAggs_vs_RestCells_minusBafA1`     | EA / RC −BafA1 |
| `PureAggs_vs_RestCells_plusBafA1`      | EA / RC +BafA1 |
| `CellsWith_vs_CellsWithout_minusBafA1` | CW / CWO −BafA1 |
| `CellsWith_vs_CellsWithout_plusBafA1`  | CW / CWO +BafA1 |

Files written by this review:

- `*/review.csv` — decisions, in `display_choice` and `display_label` only.
- `label_vocabulary.txt` — 95 labels (18 pre-existing reused, 73 added).
- `aggregate_contrasts_row_table.csv` — every candidate row: label, canonical
  term, per-contrast NES, FDR, number of terms collapsed and number of clusters
  merged. This is where the collapse counts live; see "Collapse counts" below.
- this file.

Replay check passes on all four files (80 / 79 / 76 / 79 representatives).

## Two interventions beyond labelling

Both change what the figure shows and are not part of the standard review.

**1. Unreviewed clusters removed from the candidate pool.** Pass 1 pre-fills
`display_choice="rep"` on *every* cluster, including low-priority ones that
carry no gene lists and are never shown to the reviewer. `plot_combined` admits
those with an auto-generated title-case name. Before the fix, 308 of 399
admitted labels came from such clusters and duplicated reviewed rows cell for
cell — "Proton Transmembrane Transport" (+3.42) beside "Oxidative
phosphorylation / respiratory chain" built from the same term. `display_choice`
was blanked on every row with `review_priority != HIGH`: 168, 143, 72 and 23
rows. Candidate labels fell from 399 to 91.

This is a property of the pipeline, not of this dataset. The cross-system figure
runs the same code path.

**2. Ten representatives moved so that no term serves two labels.** A row is
keyed by its label but drawn from one canonical term, so two labels sharing a
term draw two identical rows. Twelve terms had this problem. Where a member term
named the cluster's core, the representative was moved onto it:

| moved off | moved onto | why |
|---|---|---|
| `GOBP_MEMBRANE_LIPID_METABOLIC_PROCESS` | `GOBP_GLYCOLIPID_BIOSYNTHETIC_PROCESS`, `GOBP_PHOSPHATIDYLINOSITOL_BIOSYNTHETIC_PROCESS`, `GOBP_MEMBRANE_LIPID_BIOSYNTHETIC_PROCESS` | the term is sphingolipid-driven in EA/RC +BafA1 and GPI-driven in the CW contrasts |
| `GOBP_DOUBLE_STRAND_BREAK_REPAIR`, `GOBP_NUCLEOTIDE_EXCISION_REPAIR` | `REACTOME_ATP_DEPENDENT_CHROMATIN_REMODELERS` | SWI/SNF cores under repair names |
| `REACTOME_EPIGENETIC_REGULATION_OF_GENE_EXPRESSION` | `REACTOME_RNA_POLYMERASE_I_TRANSCRIPTION` | Pol I cores |
| `GOBP_RRNA_METABOLIC_PROCESS` | `GOBP_RRNA_TRANSCRIPTION` | separates Pol I transcription from ribosome biogenesis |
| `REACTOME_CHROMATIN_ORGANIZATION`, `REACTOME_REGULATION_OF_ENDOGENOUS_RETROELEMENTS` | `REACTOME_HDACS_DEACETYLATE_HISTONES`, `..._BY_KRAB_ZFP_PROTEINS` | NuRD cores |
| `GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS` | `GOBP_SPLICEOSOMAL_SNRNP_ASSEMBLY`, `GOBP_U2_TYPE_PRESPLICEOSOME_ASSEMBLY` | spliceosome cores under a generic RNP name |
| `REACTOME_SPHINGOLIPID_METABOLISM` | `GOBP_MEMBRANE_LIPID_CATABOLIC_PROCESS` | lysosomal hydrolase core (HEXA/HEXB/NEU1) |
| `GOBP_CELLULAR_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN` | `GOBP_CELLULAR_RESPONSE_TO_UNFOLDED_PROTEIN` | separates UPR from ERAD |
| `GOBP_DNA_TEMPLATED_TRANSCRIPTION_INITIATION` | `GOBP_POSITIVE_REGULATION_OF_GENE_EXPRESSION_EPIGENETIC` | COMPASS core |
| `GOBP_FATTY_ACID_METABOLIC_PROCESS` | `GOBP_MONOCARBOXYLIC_ACID_BIOSYNTHETIC_PROCESS` | separates synthesis from β-oxidation |

Where no member named the core, the labels were merged instead (below).

## Clusters excluded (7)

| contrast | representative | NES | reason |
|---|---|---|---|
| EA/RC +Baf | `GOBP_REGULATION_OF_CILIUM_ASSEMBLY` | +2.38 | incoherent core: 4 genes mixing centriolar (CCP110, CEP97) with vesicle proteins (SCFD1, SEC22B) |
| CW/CWO −Baf | `GOBP_REGULATION_OF_BLOOD_PRESSURE` | +2.12 | wrong lineage: no vasculature in a HEK293 monolayer |
| CW/CWO −Baf | `GOBP_PANCREAS_DEVELOPMENT` | +1.95 | wrong lineage: developmental programme absent; genes (GSK3B, CTNNB1, WFS1) are generic |
| CW/CWO −Baf | `REACTOME_SEALING_OF_THE_NUCLEAR_ENVELOPE_NE_BY_ESCRT_III` | −1.98 | core is pure α/β-tubulin and no member term names it; the same decrease is carried by the TRiC / CCT row |
| CW/CWO −Baf | `GOBP_TELOMERE_MAINTENANCE_VIA_TELOMERE_LENGTHENING` | −1.85 | chained: core mixes CCT chaperonin with H/ACA RNP (DKC1, NAF1) |
| CW/CWO +Baf | `GOBP_REGULATION_OF_LEUKOCYTE_APOPTOTIC_PROCESS` | +1.80 | wrong lineage: no leukocytes; genes are generic apoptosis regulators |

No cluster was excluded as a duplicate. Redundancy is handled by the shared
label, which merges clusters into one row.

## Relabelled rather than excluded

- Proteasome cores under signalling names (`DEFECTIVE_CFTR_CAUSES_CYSTIC_FIBROSIS`,
  `HEDGEHOG_LIGAND_BIOGENESIS`, `REGULATION_OF_ORNITHINE_DECARBOXYLASE_ODC`, the
  NF-κB / TCR / BCR set) → **Proteasome / UPS induction**.
- SARS-CoV-2 spike maturation = RPN1/2, STT3A, DDOST, MOGS, CANX →
  **N-linked glycosylation (ER)**.
- MHC-I peptide loading = calnexin / calreticulin / PDIA3 →
  **Calnexin / calreticulin glycoprotein QC**.
- `REACTOME_ADIPOGENESIS` = Mediator core →
  **RNA Pol II pre-initiation (Mediator / TFIID)**.
- `GOBP_NEURON_CELLULAR_HOMEOSTASIS` = V-ATPase →
  **Lysosome acidification (V-ATPase)**.
- `GOBP_VASCULAR_TRANSPORT` = SLC transporters → **Amino-acid transport (SLC)**.
- Immunoglobulin / B-cell terms = MSH2, MSH6, MLH1, EXO1 → **Mismatch repair**.
- HIV transcription terms = POLR2 subunits →
  **RNA Pol II transcription elongation**.

## Open questions

1. **Three clusters named for splicing whose core is Pol II** (EA/RC −Baf d01,
   EA/RC +Baf d22, CW/CWO −Baf d02). Every member term names splicing or FGFR
   signalling, so no Pol II-named member exists to move onto, and
   `REACTOME_MRNA_SPLICING` does contain POLR2 genes. Merged into
   **mRNA splicing / spliceosome**. Resolved by inspecting the leading edge of
   `REACTOME_MRNA_SPLICING` in these contrasts: if it is POLR2 with no SF3/SNRP
   genes, the row is transcription and the merge hides a separate Pol II decrease.
2. **`GOBP_NUCLEOSOME_ORGANIZATION`** carries a mixed FACT + SWI/SNF + histone
   chaperone core in all three contrasts where it represents a cluster. Labelled
   **Nucleosome / chromatin assembly**; **Chromatin remodelling (SWI/SNF)** is
   kept for clusters whose representative is a remodeller-specific term.
3. **`GOBP_INTRACELLULAR_PROTEIN_TRANSMEMBRANE_TRANSPORT`** (EA/RC +Baf) holds
   two modules: Sec61 translocon and peroxisomal PEX import. Labelled
   **Protein translocation (Sec61 / peroxisomal import)**.
4. **`GOBP_INTRACILIARY_TRANSPORT`** (CW/CWO +Baf, −1.70) has a clean IFT core.
   HEK293 cells form primary cilia, so it was kept. This is the one lineage call
   that could go either way.
5. **`GOBP_PROTEIN_QUALITY_CONTROL_FOR_MISFOLDED_PROTEINS`** (CW/CWO +Baf) mixes
   cytosolic (VCP, RNF126, ATXN3, UBR4) with mitochondrial proteases (CLPP,
   LONP1, AFG3L2). Labelled **Protein quality control (VCP / mitochondrial
   proteases)**; neither half alone is accurate.

## Collapse counts

§9 drops the "(N terms collapsed)" suffix the single-contrast figure carries,
because the count is per contrast. On these four contrasts it varies by up to a
factor of two on the same row — mRNA splicing collapses 18, 18, 12 and 8 terms;
oxidative phosphorylation 12, none, 23 and 16. Rows also merge more than one
cluster inside a contrast (DNA replication merges up to 6), so the number is
"terms represented by this row in this contrast", not a cluster size. All four
numbers per row are in `aggregate_contrasts_row_table.csv`.

## Result

20 rows, 10 per direction. Up: the secretory pathway and mitochondrial membrane
machinery. Down: the nucleus. The Golgi rows separate the two fraction
comparisons — intra-Golgi and endosome→TGN transport are strongly up in excised
aggregates and flat in whole cells with aggregates.

3 of 80 drawn cells (3.8 %) fall outside a contrast's detected-size window and
are drawn as a grey ×: COPII and N-linked glycosylation in EA/RC +BafA1,
O-linked glycosylation in EA/RC −BafA1.

Methods must state k = 4: FDR is controlled within each contrast, so
"significant in at least one of four" is a display rule, not an inferential
claim.
