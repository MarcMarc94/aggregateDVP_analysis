# GSEA cluster review — prompt

Paste this to Claude together with the HIGH-priority rows of
`{output_root}/{contrast}/review.csv`, and the current
`{output_root}/label_vocabulary.txt`.

Method background is in `src/GSEA_PIPELINE.md`. This file is only the review task.

---

## Task

Gene Set Enrichment Analysis produced hundreds of significant terms. Terms that
enriched through the same genes have been grouped into **clusters**. For each
cluster in the review file, decide:

1. **Which member term represents it** in the figure, and
2. **What the row should be called**, or
3. Whether the whole cluster should be **excluded**.

Edit exactly two columns and change nothing else:

| column | values | meaning |
|---|---|---|
| `display_choice` | `rep` | plot this row as the cluster's representative |
| | *(blank)* | this row is collapsed away |
| | `exclude` | drop the **entire cluster** from the figure |
| `display_label` | free text | the row's displayed name (set on the `rep` row) |

Exactly **one** `rep` per cluster you keep. `exclude` on any single row kills the
whole cluster. A `rep` is pre-filled with an automatic suggestion — leave it if you
agree.

---

## Rules

### 1. Decide from `core_genes`, never from the term name

`core_genes` lists the genes shared by every member of the cluster, ordered by how
strongly they moved. That is what the cluster actually *is*. Database names are
frequently misleading.

Two real examples from this dataset:

- A cluster whose core was 24 proteasome subunits (PSMA/PSMB/PSMC/PSMD) had
  `REACTOME_NUCLEAR_EVENTS_MEDIATED_BY_NFE2L2` as its strongest member. Displaying
  that name would have implied NFE2L2 signalling. The correct label was
  **Proteasome / UPS induction**.
- `REACTOME_REGULATION_OF_ORNITHINE_DECARBOXYLASE_ODC` sat in the same
  proteasome cluster. Ornithine decarboxylase is a proteasome substrate, so the
  term enriches whenever the proteasome does. It must never stand in for the row.

If `core_size` is small (< 5), read `consensus_genes` instead — those are genes in
at least half the members, which is more reliable for a large cluster.

### 2. Distinguish "duplicate" from "wrong label"

These need different actions, and confusing them loses real biology.

- **Duplicate** — several terms describing the same thing. Keep one, blank the
  others. *No information lost.*
- **Wrong label** — the biology is real but the database name is misleading.
  **Relabel and keep.** Example: a cluster with core AP2A1/AP2A2/AP2B1/AP2M1 +
  CLTA/CLTB/CLTC was headed `REACTOME_EPH_EPHRIN_MEDIATED_REPULSION_OF_CELLS`.
  That is clathrin machinery. Correct action: set `rep` on
  `GOBP_CLATHRIN_DEPENDENT_ENDOCYTOSIS` and label it
  **Clathrin-mediated endocytosis**. Excluding it would have discarded a genuine
  finding.
- **Wrong lineage** — the genes cannot mean what the term says in this biological
  system. **Exclude.** See rule 3.

### 3. Wrong-lineage exclusions are a system fact, not a blacklist

Judge against the system the experiment was run in, stated at the top of the
review request. For **HEK293** (kidney epithelial): no synapses, no immune cells,
no vasculature, no skin, no viral infection. A term like
`REACTOME_PROTEIN_PROTEIN_INTERACTIONS_AT_SYNAPSES` driven by FLOT1/FLOT2/DLG1 is
generic membrane scaffolding wearing a neuronal name → exclude.

**This inverts on other systems.** In brain tissue, synapse and vasculature terms
are real biology and must be kept. Always apply the rule for the system named in
the request, never a fixed list.

### 4. Reuse labels from the controlled vocabulary

Row identity in the multi-contrast figure **is** the label. If you call the same
biology "Proteasome / UPS induction" in one contrast and "UPS" in another, one
story silently becomes two rows.

- Reuse an existing vocabulary entry wherever it fits.
- Mint a new label only for genuinely new biology.
- Keep labels short, specific, and mechanism-first: *"Mitochondrial protein import
  (TOM / TIM)"*, not *"transport"*.

### 5. Prefer the strongest member as representative

Default is highest |NES|, which makes the plotted value the one that earned the
row its slot. Override it when the strongest member's *name* is badly misleading
and a nearly-as-strong member is clearer — but prefer relabelling over demoting,
because demoting lowers the plotted NES.

### 6. Review each contrast on its own numbers, in both directions

Work through one contrast at a time and rank its clusters by **its own** NES.
Never label from a ranking pooled across contrasts: NES magnitude tracks how
widely a contrast separates, so a pooled list is ordered by statistical depth and
the shallowest contrast is reviewed last or not at all. Only a labelled cluster
can become a row, so an under-reviewed contrast produces a near-empty column that
reads as "little happens here".

This is not hypothetical. Reviewing the six iPSC contrasts from one pooled list
left 3 of 22 significant up clusters labelled in the shallowest contrast; its
column in the two-contrast figure held one filled cell while 61 of its terms were
significant and up in that direction.

Cover **both directions** of every contrast before reporting. Fixing only the
direction that a figure happened to expose repeats the same defect.

Check the result with:

```python
from gsea_pipeline import report_review_coverage
report_review_coverage(contrasts, cfg, min_labelled=3)
```

It prints, per contrast and direction, how many significant HIGH clusters carry a
label, and flags any block below the floor.

### 7. Every label needs a gene list, and merging is your judgement

A label that will be drawn must have an entry in `src/gsea_label_genes.yaml`: the
genes that label's name means, written from the name and before reading any
leading edge. The figure scores every tested term against it and draws the cell
only from a term that measures the module. A new label without a list stops the
plot with a `KeyError` naming it.

After the terms are chosen, `report_row_redundancy` lists row pairs whose drawn
terms share most of a leading edge. **It is advisory.** Two rows that overlap in
genes can still be two findings, and deciding that is the point of a review rather
than a threshold: keep them apart when they behave differently in any column, and
say why in the notes. Ribosome biogenesis and cytoplasmic translation are the
worked example — they share most of their genes, and LSU rRNA maturation is down
in soma (−2.18) where the ribosomal-protein set is not (+1.30), so merging them
would erase a compartment difference.

### 8. Watch for chained clusters

`min_pairwise_overlap` near the threshold together with a small `core_size` on a
large cluster means members may have been linked through intermediates rather than
sharing one core. If the core does not look like one coherent module, say so rather
than forcing a label — it is fine to exclude and flag it.

---

## Useful columns

| column | meaning |
|---|---|
| `core_genes` | genes in **every** member, strongest first — the ground truth |
| `consensus_genes` | genes in ≥ 50% of members — use when `core_size` < 5 |
| `le_genes_top` | this term's own leading edge, strongest first |
| `family_label` | automatic gene-family guess. **A hint only** — verify against `core_genes` |
| `cluster_size` | terms in this cluster |
| `core_size` / `min_pairwise_overlap` | cluster coherence diagnostics |
| `NES`, `nes_min`, `nes_max` | this term's score, and the cluster's spread |
| `detected_size` | genes of this term measured in the experiment |
| `pheno_robust` | survives the stricter phenotype-permutation null. Informational; does **not** gate anything |
| `review_status` | `ok` / `NEW` / `MOVED` — re-check `NEW` and `MOVED` rows |

---

## Output

Edit `review.csv` in place, then **run the replay check**:

```python
import sys; sys.path.insert(0, 'src')
import pandas as pd
from gsea_pipeline import assert_pass2_selection
assert_pass2_selection(pd.read_csv(REVIEW_PATH, keep_default_na=False))
```

This must pass before plotting. It exists because a bug once deleted four chosen
rows from the figure silently — the only symptom was a missing row spotted by eye.

Then **write the report to `{output_root}/REVIEW_NOTES_{name}.md`** as well as
into the chat, where `{name}` identifies the contrast set. A decision file records
*what* was decided; it cannot record *why*, and six months later the why is the
part that is needed. The notes cover:

- the contrasts reviewed, and which figure column each becomes;
- the rows kept, with label and NES per contrast, up and down;
- which clusters were excluded and why (incoherent core / wrong lineage);
- clusters relabelled because the database name was misleading, naming the core
  genes that decided it — these are the decisions a referee will question;
- any cluster you were unsure about, and what would resolve it;
- any **new** labels added to the vocabulary, and which existing ones were reused;
- **any edit beyond labelling** — a representative moved, a decision applied to
  rows outside the HIGH block — with the reason. These change what the figure
  shows and are invisible in the decision file.

Also write `{output_root}/{name}_row_table.csv`: one row per candidate label with
its canonical term, and per contrast the NES, FDR, number of terms collapsed and
number of clusters merged. The multi-contrast figure drops the collapse count
from the axis because it differs per contrast (§9), so this table is the only
place the number survives.

Do **not** edit any column other than `display_choice` and `display_label`.
