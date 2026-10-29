# GSEA pipeline — specification

Gene Set Enrichment Analysis per contrast, with automated redundancy collapse, a
human/LLM review step, and single- or multi-contrast dot plots.

Implemented in `src/gsea_pipeline.py`, driven from the analysis notebooks.

**Vocabulary, used consistently throughout.** A **term** is one gene set
(`REACTOME_MRNA_SPLICING`). A **cluster** is a group of redundant terms. A **row**
is what appears in the figure. Every count states which granularity it refers to.

**Abbreviations.** ES = Enrichment Score. NES = Normalized Enrichment Score.
LE = leading edge. FDR = false discovery rate. ORA = over-representation analysis.

---

## 1. Ranking (per contrast)

GSEA is inherently two-class. Every comparison is a pairwise contrast; a
multi-condition analysis is several contrasts shown side by side.

- Subset to the two conditions being compared.
- **Collapse duplicate gene symbols** by fewest missing values, tie-broken by
  highest mean abundance. Explicitly *not* max-variance or max-|log2 fold change|
  — both are noise-biased, because the most variable row for a symbol is usually
  the least reliable one.
- **Rank by the Welch t-statistic** (unequal variance),
  `(mean_A - mean_B) / sqrt(s_A^2/n_A + s_B^2/n_B)`, computed **ignoring missing
  values** and requiring **at least 3 real values per group**. No imputation.
  Verified in `gseapy` source that `method='t_test'` is genuinely Welch, not pooled.
- **`run_contrast_from_rank`** takes a ranking computed in the notebook instead, for
  designs that are not two independent groups of samples: the signed z of the
  donor-level eBayes fit (BrainTissue), the BafA1 shift of the log2 fold changes
  (Autophagy), the eBayes moderated t of the brain contrasts (CrossSystem). The
  ranking is described in `params.json`.

### Two gene universes, deliberately not reconciled

Applies only when the phenotype null runs (§3). The primary analysis needs only one number per gene, so it runs on every gene with
≥ 3 real values per group. The phenotype-permutation check re-shuffles sample
labels and recomputes from scratch, so it needs a complete rectangular matrix and
runs on the smaller imputed set (≥ 10/12 observed, remaining gaps filled with the
gene's overall mean).

These are **not** forced to match. Shrinking the primary result to fit a check that
no longer appears in the figure would discard genes for nothing. Both counts are
recorded in `params.json` and the difference gets one methods sentence.

The imputation caveat — filling shrinks the between-class difference but also the
within-class standard deviation, so its net effect on the Welch t is not
one-directional — therefore applies **only** to the phenotype check, never to the
primary result.

---

## 2. Gene sets

- MSigDB: Reactome pathways + GO biological process; BrainTissue uses GO cellular
  component, cached separately (`src/gsea_genesets_gocc_cache.parquet`).
- **Artefact terms** (`artifact_pathways`: lineages absent from the system) are
  tested with all other terms, so they count in the FDR, and are removed before
  clustering and display.
- **Cached to disk in the repo** (`src/gsea_genesets_cache.parquet`), loaded from
  cache by default, remote fetch as a documented fallback. MSigDB is versioned and
  `dc.op.resource()` fetches live, so without a cache the figure cannot be
  reproduced later. The fetch is also slow.
- Gene set size is **detected-restricted**: `|annotation ∩ detected genes|`,
  keeping **20 ≤ size ≤ 500**. Used for both the size filter and the dot size, so
  dots encode evidence actually present in the data rather than annotation bulk.

---

## 3. Nulls — gene permutation in every notebook

**Gene permutation** (`gseapy.prerank`) on the ranking, giving NES, FDR, nominal p
and the leading edge. Powered in both directions. This selects, ranks and colours
everything, and it is the only null the notebooks of this repository run.

**Phenotype permutation** (`gseapy.gsea`, `permutation_type='phenotype'`) is
available in `run_contrast` with `GSEAConfig(run_phenotype=True)`; every notebook
sets `run_phenotype=False`. When run, it is stored as `fdr_pheno` and
`pheno_robust` and **never encoded in the figure**: at 6 vs 6 it discriminates
roughly 4 of 10 up rows and 0 of 10 down rows, and drawing it would imply the down
programmes are weak when the null cannot resolve them.

The two nulls ask different questions. Gene permutation asks "is this set enriched relative to random sets
of the same size", which co-regulated modules (ribosome, oxidative phosphorylation,
spliceosome) can answer trivially. Phenotype permutation preserves gene-gene
covariance and asks "would random relabelling produce this", which is the question
you actually care about — but it is severely underpowered at small n.

### Reporting language for the phenotype null (precise forms)

- **Permutations.** *1000 phenotype permutations drawn with replacement from the
  462 distinct relabellings that exist at 6 vs 6, so the null is saturated and the
  resolution ceiling is ≈ 2.2 × 10⁻³.* Not "1000 permutations" (overstates
  resolution) and not "all 462 unique relabellings" (`gseapy` samples randomly, it
  does not enumerate; exhaustive enumeration would have to be written separately).
- **Down direction.** *No down-regulated programme reached FDR < 0.05 under
  phenotype permutation, because the permuted negative NES distribution is as wide
  as the observed values at 6 vs 6.* They **were** tested and returned FDR ≈ 0.2–1.0;
  "not testable" overstates.
- **Robustness summary.** *n of m primary hits also survive phenotype permutation* —
  this is **term-level**, and must say so.

---

## 4. Clustering — Pass 1, per contrast, drops nothing

- Candidates: FDR < 0.05 on the primary null, minus the artefact list.
- Similarity: **overlap coefficient on the leading edge**, `|A∩B| / min(|A|,|B|)`.
  Not Jaccard — Jaccard dilutes with size asymmetry and misses containment. A
  68-gene leading edge sitting ~60 genes inside a 186-gene one gives Jaccard 0.31
  (escapes) but overlap 0.88 (correctly collapsed).
- **Threshold and linkage.** `distance = 1 − overlap coefficient`; **complete
  linkage**; cut at `1 − t` with **t = 0.5**. Because linkage is complete this is
  exactly equivalent to requiring *every* pairwise overlap within a cluster to
  exceed `t` — no chain of partial overlaps can merge two terms that do not
  themselves overlap.
- **Why not single linkage.** It chains catastrophically here: 377 of 402 up-terms
  merged into one cluster with an empty core. Raising the threshold does not help,
  because containment already scores 1.0.
- **Within NES direction only** — an up term can never collapse into a down term.
- Empty leading edge → dropped, counted, asserted. Never an unconditional keep.
- **Null-safe parsing.** `str(numpy.nan)` yields `'nan'`, which would produce the
  one-element set `{'nan'}` and give every missing term overlap 1.00 with every
  other missing term — one spurious mega-cluster. Nulls are checked before
  stringifying.
- `cluster_id` is **globally unique across directions**, because Pass 2 keys
  exclusions and de-duplication on it alone (§10).

**Per cluster:** core (intersection of all members' leading edges), consensus
(genes in ≥ 50% of members), `core_size`, `min_pairwise_overlap` (chaining
diagnostic), `nes_min`, `nes_max`.

**Representative:** `|NES| desc → FDR asc → smaller leading edge → term name`.
Deterministic, and it makes the plotted value the one that earned the row its slot.
(Best-FDR-first could plot a weak member of a strong cluster.)

**Gene list ordering:** by |Welch t| **in the term's own direction** — up-terms
highest t first, down-terms most negative first — truncated to the top 20. Never
alphabetical: alphabetical implies the first genes shown are the important ones,
and biases toward early-alphabet families (ACTB, AP2\*, ATP5\*).

---

## 5. Family labels — a review aid only

`family_label` is a **hint for the reviewer**. It never affects clustering,
selection, or plotting. It must still be right, because a wrong label is worse
than no label.

### Fallback for small cores

Complete linkage guarantees pairwise overlap but **not** a non-empty intersection.
As members accumulate the intersection shrinks monotonically, so a 28-member
cluster may have a 3-gene core while its consensus set stays rich. A family label
computed from a 3-gene core is fragile — one atypical gene flips the result.

**If `core_size < 5`, compute the label from consensus instead.** Both are already
computed; `core_size` is the diagnostic and this is the action attached to it.

### Three-tier resolution: exact names → prefixes → unlabelled

Exact names must be resolved **first**, because prefixes collide. The buckets, their
prefixes and their exact members are listed in `src/gsea_families.yaml`.

**Four collisions the exact tier exists to resolve.** `UBA3` and `UBE2M` match the
ubiquitin prefixes but are neddylation enzymes. `RPS27A` and `UBA52` match the
ribosomal prefixes but are ubiquitin fusion proteins.

**An example of the failure.** A dictionary with `COPS`, `UBA` and `UBE` in the
proteasome bucket scores `GOBP_PROTEIN_NEDDYLATION`'s 13-gene core at
COPS3/4/5/7A/8 + UBA3 + UBE2M = 7/13 = 54%, over the 50% floor, and labels the
COP9-signalosome cluster "Proteasome / UPS" on every run. A strict `PSM*`
definition scores 0%. Adjacent machineries are where prefix schemes fail.

**Why two rows have no prefixes.** Functionally-defined families share no naming
pattern. A prefix-only scheme is silently biased toward systematically-named
machinery and would never label them — which matters here, because NRF2/KEAP1 is
the top up cluster and p62/NCOA4 recur across the other notebooks.

**Scoring:** argmax over buckets of the fraction of core (or consensus) genes
matched, accepted only above 50%. Not first-match — a core can match two buckets
and order should not decide.

The dictionary lives in `src/gsea_families.yaml`, not in code. It carries real
biological judgement and will need extending per experiment.

---

## 6. Files

```
{output_root}/{contrast}/
    {date}_raw.csv          all tested terms, leading edges
    {date}_terms_full.csv   significant terms + clusters + full LEs (supplementary)
    review.csv              persistent decision file, no date; shipped with the code
    {date}_dotplot.svg/png  plot_contrast(save=True)
    params.json             thresholds, ranking, gene-universe size, cache, n per group
{output_root}/{date}_{save_stem}.svg/png and _grid.csv   plot_combined(save=True);
                            the grid is the row table

src/gsea_pipeline.py
src/gsea_families.yaml              family dictionary (review aid)
src/gsea_label_genes.yaml           gene list of each row label (purity selection)
src/gsea_genesets_cache.parquet     Reactome + GO biological process
src/gsea_genesets_gocc_cache.parquet  GO cellular component
src/GSEA_PIPELINE.md                this document
GSEA_REVIEW_PROMPT.md               the review task and rules
```

---

## 7. Review

The decision file is sorted by **term** NES within direction — up block, then down
— with a cluster tag on each row. Only HIGH-priority rows (the top 20 clusters per
direction) are handed to the reviewer; the rest stay in the file, inert and without
gene lists, so the file stays small enough to paste.

Two editable columns: `display_choice` (`rep` / blank / `exclude`) and
`display_label`.

Rules, stated in `GSEA_REVIEW_PROMPT.md`:

- **Decide from `core_genes`, never from the term name.** The Ornithine
  Decarboxylase and Eph-ephrin cases both hinge on this.
- **Wrong-lineage exclusions are a system fact, not a pathway blacklist.** For
  HEK293: no synapses, no immune cells, no vasculature. This *inverts* on brain
  tissue, so it is a per-experiment parameter.
- **Distinguish "duplicate" from "wrong label".** Clathrin-mediated endocytosis is
  real biology under a misleading Reactome name → relabel and keep. Synapse terms
  in HEK293 are wrong-lineage → exclude. Confusing the two loses real biology.
- **Reuse an existing label** wherever it fits; mint a new one only for new
  biology. Row identity in the multi-contrast figure is the label, so
  inconsistent labels silently split one story into two rows.

Decisions are keyed by **term name**, never by cluster number (ids renumber).
`review_status` is `ok` / `NEW` / `MOVED`, computed from a hash of cluster
membership, so a term whose cluster changed is re-flagged rather than inheriting a
decision that may no longer mean what was intended.

**After any decision write, re-run the Pass 2 replay assertion** (§10).

---

## 8. Plot — single contrast

- Rows = clusters, labelled `Name (N terms collapsed)`. Horizontal NES dot plot.
- **Row NES = the maximum across cluster members**, never an average. Members are
  the same gene module in progressively vaguer annotation wrappers, so NES falls
  monotonically with detected size — the proteasome cluster runs from +2.63 at 66
  detected genes to +1.59 at 482. Averaging would blend the cleanest measurement
  with deliberately diluted ones and understate the effect in proportion to how
  many parent terms the database happens to define, which is an artefact of
  database structure rather than biology. `nes_min` / `nes_max` are carried in the
  table so the selection-maximum is transparent.
- **Selection.** Candidates are already FDR < 0.05 on the primary null, so all rows
  are significant by construction. Within that pool take the **top 10 per direction
  by NES**. Phenotype significance never selects or gates.
- **Uniform dot styling** — colour = NES, size = detected set size. **No
  significance legend**, because every row is significant by construction; the
  caption states *all terms shown are FDR < 0.05 by gene permutation*.
  This is **specific to the single-contrast figure**. In §9 the row-admission rule
  and the cell-drawing rule differ, so a legend is required there. The asymmetry
  between the two figures is intentional, not an oversight.
- Warn, never silently shrink: unreviewed `NEW`/`MOVED` rows, clusters with no
  representative or more than one, fewer than 10 rows in a direction after
  exclusions.

---

## 9. Plot — multi-contrast

- **Row identity = the reviewed label**, not the cluster. A *term* has a well-defined score in every contrast; a *cluster* does
  not, because leading edges are ranking-dependent and clusters differ per contrast.
- **Anchor term per label** (`anchor_rule="coverage"`, the default):
  `significant in the most contrasts → max |NES| → FDR asc → term name`. It names
  the row in the row table and, under `cell_terms="canonical"`, fills every cell.
  `anchor_rule="peak"` takes `max |NES|` first, so the contrast with the widest
  separation chooses the gene set every other column is measured with.
- **One gene set per row, chosen for how well it measures the module**
  (`cell_terms="purity"`, the default). Every term tested is scored against the
  label's own gene list in `src/gsea_label_genes.yaml`:

      purity   = |leading edge ∩ list| / |leading edge|
      coverage = |leading edge ∩ list| / |list|

  Terms clearing `purity_floor` (0.45) and `coverage_floor` (0.10) in a contrast
  are candidates. A candidate that is significant in some contrast qualifies only
  if its leading edge reaches `cell_purity_floor` (0.30) in every contrast where
  it is significant; among the qualifying terms the highest purity wins, then
  |NES|. The chosen term is read in every column, so the columns compare one gene
  set. If no candidate is significant anywhere, the purest candidate is used and
  its cells are drawn hollow. A label with no candidate in any contrast is not a
  row, and is named in the run output. A cell whose term was not tested in that
  contrast is drawn as a grey × meaning **not measured**.
  `cell_terms="canonical"` reads the anchor term in every column.
- **Collapse counts are dropped** — they are per-contrast and meaningless on a
  shared row.
- **Trim:** keep labels reaching FDR < 0.05 in ≥ 1 of k contrasts. Each contrast
  is then guaranteed `min_rows_per_contrast` rows per direction, ranked by ITS OWN
  |NES| among the rows significant in it; the remaining slots go by max |NES|
  across contrasts.

### Why purity is scored against a written gene list

Scoring a term against its own cluster's core genes is circular: the core is the
intersection of the members' leading edges, so a member scores well against a core
it helped define. It also cannot reach a better term in another cluster, and in
practice that is where the better term sits — `REACTOME_DARPP_32_EVENTS`, the TRiC
term and `REACTOME_HSF1_ACTIVATION` were each in a cluster that received no
representative at all. The list in `gsea_label_genes.yaml` is written from the
label by name, before any leading edge is read, which is what lets it test the
labelling rather than confirm it.

A term's NAME carries no weight in this. `REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS`
anchors the proteasome row because its leading edge here is 100 % proteasome
subunits in one contrast and 87 % in another — higher than `PROTEASOME_ASSEMBLY`
itself. The row label says what the row is; the row table records what measured it.

### Why the quota exists

Ranking rows by |NES| pooled across contrasts hands the row budget to whichever
contrast separates most widely. In the two-contrast iPSC figure the down block cut
at |NES| 2.21 while no module specific to the shallower contrast reached 2.09, so
that column could not contribute one row of its own however well it was reviewed.
`min_rows_per_contrast` (default 3) fixes a floor per contrast and direction.
`k × min_rows_per_contrast` above the row budget is reported, not silently
truncated.

**Settings of the figures.** Autophagy and CrossSystem: `cell_terms="canonical"`,
`anchor_rule="peak"`, `min_rows_per_contrast=0`. iPSCs: the defaults (`"purity"`,
`"coverage"`, 3 rows per contrast). Both settings read the same decision files and
raw tables.

### Per-cell significance is required here

A row qualifies on "significant in at least one contrast", but each cell is drawn
from its own contrast — so a row significant in contrast 1 may be non-significant
in contrast 2. Three states:

| state | meaning | mark |
|---|---|---|
| filled | FDR < 0.05 in this contrast | coloured dot |
| hollow | tested, FDR ≥ 0.05 | open circle |
| absent | outside this contrast's detected-size window, never tested | light grey `×` |

The third state is necessary because the detected gene universe is per-contrast
(different samples → different genes clear the ≥ 3-per-group floor → different
detected set sizes), so a term near the 20 or 500 boundary can be tested in one
contrast and absent from another. That is categorically different from "tested and
not significant" and must not render as hollow. The count of `×` cells is reported;
more than a handful means the size window is doing something unintended.

**The legend returns on this figure only**, because row admission and cell drawing
follow different rules. On the single-contrast figure they coincide.

### Methods sentences required

- FDR is controlled **within** each contrast. "Significant in ≥ 1 of k contrasts"
  therefore carries no family-wise guarantee across contrasts; it is a
  display-selection rule, not an inferential claim. State k.
- Each row shows one gene set in every column: under `cell_terms="purity"` the set
  chosen by purity against the label's gene list (state the floors), under
  `cell_terms="canonical"` the anchor term, which may understate a contrast whose
  own member of the module scores higher. State which was used.
- Row admission uses a per-contrast quota (`min_rows_per_contrast`), so the row
  set is not a pure ranking by effect size. State the floor.

With k = 1 this path degenerates exactly to §8. One pipeline, not two.

---

## 10. Guardrails

Each of these stops or reports a fault that a run would otherwise complete without
a warning.

- **Pass 2 replay assertion.** After any decision write, replay Pass 2's exact
  selection logic and fail if any chosen representative would be dropped.
- **`cluster_id` global-uniqueness assertion.** Pass 2 keys exclusions and
  de-duplication on `cluster_id` alone, so an id shared by the two directions
  would delete rows.
- **Empty-leading-edge assertion.** A term with an empty leading edge overlaps no
  other term, so it would stand as a cluster of its own and become its
  representative; such terms are dropped and counted.
- **Null-safe leading-edge parsing** (§4).
- **No silent truncation.** Anything dropped is counted and logged.

---

## 11. Known judgement calls

Not bugs — decisions a reader should be able to see and disagree with.

- **The clustering threshold `t = 0.5`** directly controls how aggressively
  redundant families collapse. Measured on the PFF-vs-PBS contrast, up direction,
  complete linkage:

  | t | clusters | largest cluster |
  |---|---|---|
  | 0.5 | 107 | 28 |
  | 0.6 | 126 | 14 |
  | 0.7 | 148 | 13 |

  At `t = 0.7` the proteasome family splits into several rows instead of one. The
  value belongs in the methods, with the sensitivity acknowledged.

- **The artefact list is experiment-specific** and inverts between cell lines and
  tissue. Synapse and vasculature terms are artefacts in HEK293 and real biology in
  brain.
- **The family dictionary encodes biological boundaries**, particularly 26S
  proteasome vs cullin-RING ligase / neddylation, and structural ribosome vs
  biogenesis factors.
- **Representatives and labels are chosen by a human or an LLM** reading
  `core_genes`, and recorded in `review.csv`.
- **Row NES is a maximum over correlated members**, which inflates slightly; the
  range is reported alongside.
- **The two nulls run on different gene universes** (§1), stated rather than
  reconciled.
