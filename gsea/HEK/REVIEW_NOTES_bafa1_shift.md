# GSEA Pass 2 review — BafA1-shift ranking (section 19)

Reviewed 2026-09-25 by Claude, following `GSEA_REVIEW_PROMPT.md` and
`src/GSEA_PIPELINE.md` (§7 review, §8 single contrast).

Contrast folder: `BafA1_shift_PureAggs_vs_RestCells/`.
System: HEK293 monolayer.

What was ranked, from `params.json`: `bafa1_shift = log2FC(Pure Aggregates /
Rest Cells, +BafA1) − log2FC(Pure Aggregates / Rest Cells, −BafA1)`, symbols
collapsed by highest mean log2 intensity. 6868 genes, 2662 sets tested.
A term's direction says how lysosomal inhibition changes aggregate composition:

- **up** = enriched in aggregates MORE when the lysosome is blocked
- **down** = enriched LESS

Files written by this review:

- `BafA1_shift_PureAggs_vs_RestCells/review.csv` — decisions, in
  `display_choice` and `display_label` only.
- `bafa1_shift_row_table.csv` — all 26 kept clusters: term, label, direction,
  NES, FDR, terms collapsed, and whether the row reaches the figure.
- `label_vocabulary.txt` — 5 labels added, now 100.
- this file.

Replay check passes: all 26 representatives survive selection.

## The asymmetry is the finding

41 clusters cleared FDR < 0.05: **3 up and 38 down**. Blocking the lysosome does
not add a programme to the aggregates so much as it removes the relative excess
of everything nuclear, mitochondrial and ribosomal. The three up clusters are
the substantive positive result and all three are membrane or degradative:

| label | NES | FDR | terms | core |
|---|---|---|---|---|
| SNARE-mediated membrane fusion | +2.12 | 0.0056 | 2 | HOPS (VPS11/18/39/41), STX6/10/16, VAMP3, CLN3, ATP13A2, EEA1 |
| Peptide catabolism (peptidases) | +2.00 | 0.025 | 1 | IDE, NPEPPS, ERAP1, THOP1, DNPEP, CPE, LNPEP, PAM |
| Centrosome / MTOC assembly | +1.99 | 0.019 | 2 | CEP72/131/192/250/350, PCM1, AURKA, HOOK3, DCTN2 |

The figure draws 10 down rows and 3 up rows; `plot_contrast` prints a NOTE that
fewer than 10 were available in the up direction. That is a real property of the
ranking, not a review artefact.

## Rows in the figure

| label | NES | FDR | terms collapsed |
|---|---|---|---|
| mRNA splicing / spliceosome | −2.92 | 0 | 7 |
| Oxidative phosphorylation / respiratory chain | −2.44 | 0 | 9 |
| Ribosome biogenesis / rRNA processing | −2.44 | 0 | 5 |
| Splicing regulation (SR / hnRNP) | −2.38 | 0 | 5 |
| rRNA transcription (Pol I) | −2.37 | 0 | 6 |
| B-WICH complex (rDNA chromatin) | −2.35 | 0 | 3 |
| NuRD complex (CHD4 / MTA) | −2.27 | 2.4e-4 | 4 |
| Respiratory chain complex assembly | −2.20 | 3.4e-4 | 1 |
| RNA Pol II transcription elongation | −2.14 | 1.1e-3 | 7 |
| COMPASS / MLL histone methylation | −2.10 | 1.6e-3 | 2 |
| Centrosome / MTOC assembly | +1.99 | 0.019 | 2 |
| Peptide catabolism (peptidases) | +2.00 | 0.025 | 1 |
| SNARE-mediated membrane fusion | +2.12 | 0.0056 | 2 |

The other 13 kept clusters are in `bafa1_shift_row_table.csv` with
`in_figure = False`; they rank below the top ten of their direction.

## Clusters excluded (15 of 41, all on the down side)

Rows are clusters in §8, so a duplicate module must be excluded rather than
merged under a shared label — the opposite of the multi-contrast figures.

**Duplicates — the same core, weaker cluster dropped (11)**

| cluster | representative | NES | kept instead |
|---|---|---|---|
| d09 | `REACTOME_RRNA_PROCESSING` | −2.21 | d03, same UTP / EXOSC core |
| d23 | `GOBP_RIBOSOMAL_LARGE_SUBUNIT_BIOGENESIS` | −1.97 | d03 |
| d15 | `GOBP_PROTEIN_RNA_COMPLEX_ORGANIZATION` | −2.12 | d01, same SNRP core |
| d25 | `GOBP_REGULATION_OF_RNA_SPLICING` | −1.91 | d04, same SR / hnRNP core |
| d31 | `GOBP_MRNA_SPLICE_SITE_RECOGNITION` | −1.80 | d04 |
| d12 | `REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT` | −2.16 | d02, complex I versus complex V of one chain |
| d10 | `REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION` | −2.21 | composite: its leading edge is the union of d02, d17 and d19 |
| d28 | `GOBP_MITOCHONDRIAL_TRANSMEMBRANE_TRANSPORT` | −1.88 | d19, same TOM / TIM core |
| d20 | `GOBP_TRANSCRIPTION_BY_RNA_POLYMERASE_I` | −2.05 | d05 and d06 |
| d13 | `REACTOME_FORMATION_OF_THE_BETA_CATENIN_TCF_TRANSACTIVATING_COMPLEX` | −2.16 | d16, same COMPASS core |
| d38 | `GOBP_DNA_TEMPLATED_TRANSCRIPTION_INITIATION` | −1.70 | d30, same Mediator / TFIID core |

**Wrong lineage or incoherent (4)**

| cluster | representative | NES | reason |
|---|---|---|---|
| d08 | `REACTOME_MEIOTIC_RECOMBINATION` | −2.24 | meiosis does not occur in HEK293; the core is four replication-dependent histones plus NBN |
| d29 | `REACTOME_RUNX1_REGULATES_TRANSCRIPTION_OF_GENES_INVOLVED_IN_DIFFERENTIATION_OF_HSCS` | −1.83 | haematopoietic programme, and the leading edge mixes core histones with proteasome β subunits — two modules, no coherent core |
| d27 | `REACTOME_REGULATION_OF_MITF_M_DEPENDENT_GENES_INVOLVED_IN_PIGMENTATION` | −1.91 | melanocytes absent; core is pure SWI/SNF, already carried by d21 |
| d26 | `GOBP_CELLULAR_RESPONSE_TO_CALCIUM_ION` | −1.91 | core is 3 genes (MICU3, NLGN1, BRAF) and the members name neurotransmitter secretion — no synapses here |

## Relabelled because the database name is misleading

The core histone octamer genes (H2BC1, H2BC12, H2BC26, H2AX, H3C15) move down
with BafA1 and drag every term that contains them. Five clusters are named for
biology that is not present:

- **d16** `REACTOME_RUNX1_REGULATES_GENES_INVOLVED_IN_MEGAKARYOCYTE_DIFFERENTIATION`
  → **COMPASS / MLL histone methylation**. Core: WDR5, DPY30, ASH2L, RBBP5,
  KMT2B/C/D, EP300 plus the histones. Megakaryocytes are absent; the complex is real.
- **d07** `REACTOME_REGULATION_OF_ENDOGENOUS_RETROELEMENTS` → **NuRD complex
  (CHD4 / MTA)**. Core: CHD4, GATAD2A, MBD3, MTA2, RBBP4/7, HDAC2, CHD3.
- **d06** `REACTOME_B_WICH_COMPLEX_POSITIVELY_REGULATES_RRNA_EXPRESSION` →
  **B-WICH complex (rDNA chromatin)**. Kept its own row rather than folding into
  Pol I: SMARCA5, BAZ1B, MYBBP1A, DDX21 and DEK are a remodeller acting at the
  rDNA, not the polymerase.
- **d32** `REACTOME_BASE_EXCISION_REPAIR` → kept as **Base excision repair**. The
  leading edge carries genuine BER enzymes (SMUG1, MBD4, APEX1, FEN1, POLE3/4,
  POLD3) alongside the histones, so this one is real biology under a correct name.
- **d05** `REACTOME_NEGATIVE_EPIGENETIC_REGULATION_OF_RRNA_EXPRESSION` →
  **rRNA transcription (Pol I)**. Core is POLR1A/C/F/G, TBP, TAF1B, UBTF: the
  polymerase, not a repressive programme.

## Clusters I was not certain about

1. **d02 versus d12** — complex V (ATP5\*) and complex I (NDUF\*) formed separate
   clusters. I kept the stronger, d02, under the label
   *Oxidative phosphorylation / respiratory chain* and excluded d12 as a
   duplicate. Both are the same respiratory chain, but if the intent is to show
   that one complex moves and the other does not, they should be two rows.
   Resolved by checking whether the two leading edges overlap at all.
2. **d34 SUMOylation (PIAS / UBE2I)** — its leading edge also carries Polycomb
   (BMI1, PHC1/2, CBX8), which is d37's territory. Kept both; the SUMO ligases
   (PIAS1, PIAS4, UBE2I, SUMO2) are the part that is unique to d34.
3. **d33 Cytoplasmic ribosome (RPL / RPS)** — the existing vocabulary entry is
   *Cytoplasmic translation (40S)*, minted on the PFF/PBS contrast. Here the
   leading edge is large and small subunit proteins in equal measure, so the
   (40S) qualifier would be wrong and a new label was minted. If the PFF/PBS
   label was also 60S-containing, the two should be merged into one entry.
4. **d21 Nucleosome / chromatin assembly** — its core mixes FACT (SSRP1,
   SUPT16H), SWI/SNF (SMARCA4/B1/C1/C2/D1/E1) and the histones. It is the row
   that absorbed the SWI/SNF cluster d27, so it is broader than its label says.

## New vocabulary labels (5)

- B-WICH complex (rDNA chromatin)
- Centrosome / MTOC assembly
- Citric acid cycle (TCA)
- Cytoplasmic ribosome (RPL / RPS)
- Peptide catabolism (peptidases)

21 existing labels were reused, so this figure and the section 16 figure name the
same biology the same way.

## Edits beyond labelling

None. Every representative is the one Pass 1 suggested; no row was moved and no
decision was applied outside the HIGH block. All 41 clusters are HIGH in this
contrast, so the unreviewed-cluster problem that affects the four-contrast figure
does not arise here.

## For the caption

All rows shown are FDR < 0.05 by gene permutation (§8), so there is no
significance legend. State the sentinel terms dropped by the run, the ranking as
`params.json` records it, and that k = 1 — no multiplicity correction spans this
figure and the others in the notebook.
