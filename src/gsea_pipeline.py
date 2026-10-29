"""
Gene Set Enrichment Analysis pipeline with redundancy collapse and human review.

Full rationale for every choice here lives in ``src/GSEA_PIPELINE.md``. Read that
before changing defaults — several of them exist because the obvious alternative
failed silently on real data.

Vocabulary used consistently: a **term** is one gene set; a **cluster** is a group
of redundant terms; a **row** is what appears in the figure.

Abbreviations: ES = Enrichment Score, NES = Normalized Enrichment Score,
LE = leading edge, FDR = false discovery rate.

Typical use from a notebook::

    import sys; sys.path.insert(0, 'src')
    from gsea_pipeline import GSEAConfig, load_gene_sets, run_contrast, build_review
    from gsea_pipeline import plot_contrast, plot_combined

    cfg = GSEAConfig(output_root=f"{output_dir}/GSEA",
                     artifact_pathways=HEK293_ARTIFACTS)
    gene_sets = load_gene_sets(cfg)
    res = run_contrast(adata, "class_renamed", cond_a, cond_b,
                       "PFF_vs_PBS_minusBafA1", gene_sets, cfg)
    build_review(res, cfg)          # -> review.csv, then review it
    plot_contrast("PFF_vs_PBS_minusBafA1", cfg)
"""

from __future__ import annotations

import hashlib
import json
import time
import warnings
from dataclasses import dataclass
from datetime import date
from itertools import combinations
from math import comb
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

# decoupler imports marsilea, which imports legendkit, and legendkit sets eight of
# matplotlib's legend rcParams when it is imported (font size 10 among them).
# rc_context restores all rcParams on exit, so this import changes no setting.
with plt.rc_context():
    import decoupler as dc

__all__ = [
    "GSEAConfig",
    "FamilyLabeller",
    "load_gene_sets",
    "welch_t_ranking",
    "build_imputed_matrix",
    "run_contrast",
    "run_contrast_from_rank",
    "build_review",
    "assert_pass2_selection",
    "report_undecided",
    "plot_contrast",
    "plot_combined",
    "load_label_genes",
    "report_offlist_cells",
    "report_review_coverage",
    "report_row_redundancy",
]

_MODULE_DIR = Path(__file__).resolve().parent


# ══════════════════════════════════════════════════════════════════════════
# Configuration
# ══════════════════════════════════════════════════════════════════════════

@dataclass
class GSEAConfig:
    """Every knob in one place. Nothing is read from ambient notebook scope."""

    # --- where output goes -------------------------------------------------
    output_root: Path

    # --- experiment-specific ----------------------------------------------
    #: Wrong-lineage terms for THIS system. HEK293 has no synapses, immune cells
    #: or vasculature; on brain tissue the rule inverts, so this is a parameter
    #: rather than a module global.
    artifact_pathways: frozenset = frozenset()

    # --- gene sets ---------------------------------------------------------
    collections: Tuple[str, ...] = ("reactome_pathways", "go_biological_process")
    genesets_cache: Path = _MODULE_DIR / "gsea_genesets_cache.parquet"
    #: Detected-restricted size window. Applied to |annotation ∩ detected|,
    #: not to the raw annotation size.
    min_size: int = 20
    max_size: int = 500

    # --- ranking -----------------------------------------------------------
    #: Real (non-missing) values required per group for the Welch t. No
    #: imputation is used for the primary statistic.
    min_valid_per_group: int = 3
    #: Completeness required for the *phenotype* matrix only, which needs a
    #: rectangular array. Remaining gaps get the gene's overall mean.
    min_valid_frac_pheno: float = 10 / 12

    # --- permutation -------------------------------------------------------
    #: What to do with terms that come back with the NES=1.0 / p=1.0 sentinel.
    #: "raise" halts, "drop" removes them and reports. Never "keep": the value is
    #: not a statistic, and a degenerate term can win a cluster representative.
    sentinel_action: str = "raise"
    n_perm: int = 1000
    seed: int = 42
    threads: int = 4
    run_phenotype: bool = True

    # --- clustering --------------------------------------------------------
    fdr_cut: float = 0.05
    #: distance = 1 - overlap coefficient; complete linkage; cut at 1 - t.
    #: With complete linkage this is exactly "every pairwise overlap in a
    #: cluster exceeds t". See GSEA_PIPELINE.md §11 for sensitivity.
    overlap_thresh: float = 0.5
    n_review_clusters: int = 20
    top_k_genes: int = 20

    # --- family labels -----------------------------------------------------
    families_path: Path = _MODULE_DIR / "gsea_families.yaml"

    # --- plotting ----------------------------------------------------------
    top_n_plot: int = 10
    #: Dot area in points squared, at the smallest and the largest detected set
    #: size. Area, not diameter, so the drawn disc is proportional to the number
    #: the key states. Raise both on a large canvas, where a 200 pt2 dot (16 pt
    #: across) is a speck. Keep hi/lo within about one order of magnitude, or the
    #: smallest sets stop being visible at all.
    dot_area_min: float = 30
    dot_area_max: float = 200
    #: When > 0, the dots are rescaled after layout so that the LARGEST one's
    #: diameter is this fraction of the row pitch, and dot_area_min/max set only
    #: their ratio. Areas are in points squared and do not scale with the canvas,
    #: so a panel made twice as tall otherwise keeps pin-sized dots. 0 disables
    #: it and the areas are taken literally. Single-contrast figure only.
    dot_fraction_of_row: float = 0.0
    cmap_colors: Tuple[str, str, str] = ("#66c2a5", "white", "#fc8d62")
    label_fontsize: int = 11
    tick_fontsize: int = 9
    title_fontsize: int = 13
    #: Line widths and key marker size of the multi-contrast figure
    #: (plot_combined), in points on the canvas. A notebook that draws on a
    #: canvas larger than the print passes its own values, so that the printed
    #: widths stay within the journal's limits.
    edge_lw_sig: float = 0.35       # outline of a significant dot
    edge_lw_ns: float = 0.6         # outline of a tested, non-significant dot
    cross_lw: float = 1.1           # cross of a never-tested cell
    grid_lw: float = 0.5
    key_edge_lw: float = 0.3        # outline of the set-size key dots
    key_marker_size: float = 8      # markers of the significance key
    # Outlines of the two significance-key markers. None leaves them at
    # rcParams["lines.markeredgewidth"].
    key_mew_sig: Optional[float] = None
    key_mew_ns: Optional[float] = None
    #: Line widths of the single-contrast figure (plot_contrast), in points on
    #: the canvas. Its size-key outlines use key_edge_lw above.
    single_edge_lw: float = 0.4     # outline of a dot
    single_zero_lw: float = 0.6     # the NES = 0 line
    #: Layout for narrow panels. keys_position: "right" puts the keys and the
    #: colourbar in a strip right of the plot. "inside" (single-contrast figure)
    #: puts the size key in the upper left of the axes, which is empty because the rows are
    #: sorted by NES: up rows have their dots on the right, down rows on the
    #: left. "corner" (multi-contrast figure) puts both keys, as one key with two
    #: columns, in the empty corner below the row labels and left of the contrast
    #: labels, and writes "NES" under the colourbar instead of beside it.
    keys_position: str = "right"
    #: With keys_position="inside": the corner of the axes the key goes in.
    keys_inside_loc: str = "lower right"
    #: Single-contrast figure: draw the NES colourbar. NES is also the x axis.
    nes_colorbar: bool = True
    #: Single-contrast figure: write the collapse count as "(9)" instead of
    #: "(9 terms collapsed)". The caption then has to say what the number is.
    short_collapse_label: bool = False
    #: Multi-contrast figure: rotation of the contrast labels, in degrees. At 90
    #: they read upwards, centred on their column, and their lines end at the axis.
    contrast_label_rotation: float = 45

    def contrast_dir(self, contrast: str) -> Path:
        d = Path(self.output_root) / contrast
        d.mkdir(parents=True, exist_ok=True)
        return d


def _today() -> str:
    return date.today().strftime("%Y%m%d")


def _cmap(cfg: GSEAConfig):
    return LinearSegmentedColormap.from_list("gsea_nes", list(cfg.cmap_colors))


# ══════════════════════════════════════════════════════════════════════════
# Leading-edge parsing  (null-safe: see GSEA_PIPELINE.md §4)
# ══════════════════════════════════════════════════════════════════════════

def parse_leading_edge(value) -> set:
    """Split a ';'-joined leading edge into a set of gene symbols.

    Null-safe on purpose. ``str(numpy.nan)`` is ``'nan'``, so a naive split
    returns the one-element set ``{'nan'}`` rather than an empty set — and then
    every term with a missing leading edge overlaps every other one at 1.00,
    collapsing them into a single spurious cluster.
    """
    if value is None:
        return set()
    if isinstance(value, float) and np.isnan(value):
        return set()
    return {
        g.strip() for g in str(value).split(";")
        if g.strip() and g.strip().lower() != "nan"
    }


def order_genes_by_rank(
    genes: Iterable[str],
    direction: str,
    rank_lookup: Dict[str, int],
    k: Optional[int] = None,
) -> List[str]:
    """Order genes by strength of contribution *in the term's own direction*.

    Up-terms: highest Welch t first. Down-terms: most negative first. Genes
    absent from the ranking sort last, alphabetically.

    Never order alphabetically — that implies the first genes shown are the
    important ones and biases toward early-alphabet families (ACTB, AP2*, ATP5*).
    """
    genes = list(genes)
    present = [g for g in genes if g in rank_lookup]
    absent = sorted(g for g in genes if g not in rank_lookup)
    if direction == "up":
        present.sort(key=lambda g: rank_lookup[g])       # rank 0 == highest t
    else:
        present.sort(key=lambda g: -rank_lookup[g])      # most negative first
    out = present + absent
    return out[:k] if k else out


def overlap_coefficient(a: set, b: set) -> float:
    """|A∩B| / min(|A|,|B|).

    Not Jaccard: Jaccard scales with set size, so a large leading edge dilutes
    the ratio below threshold and containment escapes. A 68-gene LE sitting ~60
    genes inside a 186-gene LE gives Jaccard 0.31 but overlap 0.88.
    """
    m = min(len(a), len(b))
    return (len(a & b) / m) if m else 0.0


# ══════════════════════════════════════════════════════════════════════════
# Family labelling  (review aid only — never affects stats or plotting)
# ══════════════════════════════════════════════════════════════════════════

class FamilyLabeller:
    """Three-tier gene-family labeller: exact names -> prefixes -> unlabelled.

    Exact names resolve first because prefixes collide: ``UBA3``/``UBE2M`` match
    the ubiquitin prefixes but are neddylation enzymes; ``RPS27A``/``UBA52``
    match the ribosomal prefixes but are ubiquitin fusion proteins.

    Scoring is argmax over buckets of the fraction of genes matched, accepted
    only above ``min_fraction`` — not first-match, since a core can match two
    buckets and ordering should not decide.
    """

    def __init__(self, families: dict, min_fraction: float = 0.5,
                 core_size_fallback: int = 5):
        self.min_fraction = min_fraction
        self.core_size_fallback = core_size_fallback
        self._exact: Dict[str, str] = {}
        self._prefix: List[Tuple[str, str]] = []
        for name, spec in families.items():
            for g in (spec.get("exact") or []):
                self._exact[str(g).upper()] = name
            for p in (spec.get("prefix") or []):
                self._prefix.append((str(p).upper(), name))
        # longest prefix first, so PSMA beats a hypothetical PS
        self._prefix.sort(key=lambda t: -len(t[0]))

    @classmethod
    def from_yaml(cls, path: Path) -> "FamilyLabeller":
        with open(path) as fh:
            doc = yaml.safe_load(fh)
        return cls(
            families=doc.get("families", {}),
            min_fraction=float(doc.get("min_fraction", 0.5)),
            core_size_fallback=int(doc.get("core_size_fallback", 5)),
        )

    def bucket_of(self, gene: str) -> Optional[str]:
        g = str(gene).upper()
        if g in self._exact:                      # tier 1
            return self._exact[g]
        for pre, name in self._prefix:            # tier 2
            if g.startswith(pre):
                return name
        return None                               # tier 3

    def label(self, genes: Sequence[str]) -> str:
        if not genes:
            return ""
        counts: Dict[str, int] = {}
        for g in genes:
            b = self.bucket_of(g)
            if b:
                counts[b] = counts.get(b, 0) + 1
        if not counts:
            return ""
        best = max(counts, key=lambda k: counts[k])
        return best if counts[best] / len(genes) > self.min_fraction else ""

    def label_cluster(self, core: Sequence[str], consensus: Sequence[str]) -> str:
        """Label from the core, falling back to consensus when the core is thin.

        Complete linkage guarantees pairwise overlap but *not* a non-empty
        intersection, and the intersection shrinks as members accumulate. A
        label computed from a 3-gene core is fragile — one atypical gene flips
        the argmax.
        """
        if len(core) < self.core_size_fallback:
            return self.label(consensus)
        return self.label(core)


# ══════════════════════════════════════════════════════════════════════════
# Gene sets
# ══════════════════════════════════════════════════════════════════════════

def load_gene_sets(cfg: GSEAConfig, force_refresh: bool = False) -> pd.DataFrame:
    """Return a tidy gene-set table with columns ``source`` and ``target``.

    Reads the on-disk cache by default. MSigDB is versioned and
    ``decoupler.op.resource`` fetches live, so an uncached pipeline cannot be
    reproduced later — and the fetch is slow. The cache is committed to the repo.
    """
    cache = Path(cfg.genesets_cache)
    if cache.exists() and not force_refresh:
        gs = pd.read_parquet(cache)
        return gs[["source", "target"]]

    msigdb = dc.op.resource("MSigDB", organism="human")
    frames = [msigdb[msigdb["collection"] == c] for c in cfg.collections]
    gs = (pd.concat(frames, ignore_index=True)
            .rename(columns={"geneset": "source", "genesymbol": "target"})
            .drop_duplicates(subset=["source", "target"])
            .loc[:, ["source", "target"]]
            .reset_index(drop=True))
    cache.parent.mkdir(parents=True, exist_ok=True)
    gs.to_parquet(cache, index=False)
    print(f"Gene sets fetched and cached -> {cache} "
          f"({gs['source'].nunique()} sets, {len(gs)} rows)")
    return gs


def _gene_sets_dict(gs: pd.DataFrame) -> Dict[str, List[str]]:
    return gs.groupby("source")["target"].apply(list).to_dict()


# ══════════════════════════════════════════════════════════════════════════
# Ranking
# ══════════════════════════════════════════════════════════════════════════

def _collapse_symbols(df: pd.DataFrame) -> pd.DataFrame:
    """One row per gene symbol: fewest missing, tie-broken by highest mean.

    Deliberately not max-variance or max-|log2 fold change| — both are
    noise-biased, since the most variable row for a symbol is usually the least
    reliable one.
    """
    n_missing = df.isna().sum(axis=1).to_numpy()
    mean_abund = df.mean(axis=1).to_numpy()
    order = np.lexsort((-mean_abund, n_missing))   # primary key = fewest missing
    ordered = df.iloc[order]
    return ordered[~ordered.index.duplicated(keep="first")]


def welch_t_ranking(
    adata,
    group_col: str,
    cond_a: str,
    cond_b: str,
    cfg: GSEAConfig,
    symbol_col: str = "single_gene_name",
) -> Tuple[pd.Series, pd.DataFrame, dict]:
    """Welch t per gene, computed ignoring missing values. No imputation.

    ``prerank`` takes a vector, not a matrix, so the primary statistic never
    needs a complete table. Requiring only ``min_valid_per_group`` real values
    per group recovers genes that a completeness filter would drop, and keeps
    invented numbers out of the statistic that selects, ranks and colours
    everything.

    Returns ``(welch_t, collapsed_matrix, diagnostics)``. The matrix is returned
    because the phenotype null needs it.
    """
    mask = adata.obs[group_col].isin([cond_a, cond_b])
    sub = adata[mask].copy()
    if sub.n_obs == 0:
        raise ValueError(f"No samples matched {cond_a!r} / {cond_b!r} in {group_col!r}")

    expr = sub.to_df().T
    expr.index = pd.Index(sub.var[symbol_col].values, name="gene")
    expr = expr[expr.index.notna()]
    n_rows_before = expr.shape[0]
    expr = _collapse_symbols(expr)

    classes = np.array([cond_a if c == cond_a else cond_b
                        for c in sub.obs[group_col].values])
    cols_a = expr.columns[classes == cond_a]
    cols_b = expr.columns[classes == cond_b]
    n_a, n_b = len(cols_a), len(cols_b)
    if min(n_a, n_b) < cfg.min_valid_per_group:
        raise ValueError(f"Need >= {cfg.min_valid_per_group} samples per group, "
                         f"got {n_a} vs {n_b}")

    A, B = expr[cols_a], expr[cols_b]
    valid_a, valid_b = A.notna().sum(axis=1), B.notna().sum(axis=1)
    keep = (valid_a >= cfg.min_valid_per_group) & (valid_b >= cfg.min_valid_per_group)

    # pandas mean/std are NaN-aware; std is ddof=1, which is what Welch wants
    num = A[keep].mean(axis=1) - B[keep].mean(axis=1)
    den = np.sqrt(A[keep].std(axis=1) ** 2 / valid_a[keep]
                  + B[keep].std(axis=1) ** 2 / valid_b[keep])
    with np.errstate(divide="ignore", invalid="ignore"):
        t = num / den
    t = t.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)

    diag = {
        "n_protein_rows": int(n_rows_before),
        "n_symbols": int(expr.shape[0]),
        "n_samples_a": int(n_a),
        "n_samples_b": int(n_b),
        "n_genes_primary": int(len(t)),
        "pct_missing": float(100 * expr.isna().sum().sum() / max(expr.size, 1)),
    }
    print(f"Symbol collapse: {n_rows_before} protein rows -> {expr.shape[0]} symbols")
    print(f"Welch t ranking: {len(t)} genes "
          f"(>= {cfg.min_valid_per_group} real values per group, no imputation) "
          f"[{t.iloc[0]:+.2f} … {t.iloc[-1]:+.2f}]")
    return t, expr, diag


def build_imputed_matrix(
    expr: pd.DataFrame, cfg: GSEAConfig
) -> Tuple[pd.DataFrame, dict]:
    """Complete matrix for the phenotype null only.

    NOT strictly conservative: filling with the gene's overall mean shrinks the
    between-class difference (deflating the Welch t numerator) but also shrinks
    the within-class standard deviation (inflating it). The effects oppose and
    the net direction is not guaranteed. Immaterial at <1% of cells, but it must
    not be claimed as a one-way guarantee — which is why the primary statistic
    avoids imputation entirely.
    """
    n_samples = expr.shape[1]
    min_valid = int(np.ceil(cfg.min_valid_frac_pheno * n_samples))
    valid = expr.notna().sum(axis=1)
    kept = expr.loc[valid >= min_valid]
    n_cells = int(kept.isna().sum().sum())
    imputed = kept.apply(lambda r: r.fillna(r.mean()), axis=1)
    diag = {
        "n_genes_phenotype": int(imputed.shape[0]),
        "min_valid_phenotype": min_valid,
        "n_imputed_cells": n_cells,
        "pct_imputed_cells": float(100 * n_cells / max(imputed.size, 1)),
    }
    print(f"Phenotype matrix: {imputed.shape[0]} genes (>= {min_valid}/{n_samples} "
          f"observed), {n_cells} cells imputed "
          f"({diag['pct_imputed_cells']:.2f}% of the matrix)")
    return imputed, diag


# ══════════════════════════════════════════════════════════════════════════
# Running the two nulls
# ══════════════════════════════════════════════════════════════════════════

def _num(s) -> np.ndarray:
    return pd.to_numeric(s, errors="coerce").values


def run_contrast(
    adata,
    group_col: str,
    cond_a: str,
    cond_b: str,
    contrast: str,
    gene_sets: pd.DataFrame,
    cfg: GSEAConfig,
    symbol_col: str = "single_gene_name",
) -> pd.DataFrame:
    """Run both nulls for one pairwise contrast; write raw.csv and params.json.

    Primary   : gene permutation (``gseapy.prerank``) on the Welch t ranking.
                Powered in both directions; selects, ranks and colours everything.
    Robustness: phenotype permutation (``gseapy.gsea``). Stored, never plotted —
                at small n it cannot resolve the down direction, and drawing it
                would imply the down programmes are weak rather than untestable.
    """
    welch_t, expr, diag = welch_t_ranking(
        adata, group_col, cond_a, cond_b, cfg, symbol_col=symbol_col)

    detected = set(welch_t.index)
    gsd_all = _gene_sets_dict(gene_sets)
    detected_size = {pw: len(set(g) & detected) for pw, g in gsd_all.items()}
    gsd = {pw: gsd_all[pw] for pw, n in detected_size.items()
           if cfg.min_size <= n <= cfg.max_size}
    print(f"Detected-size filter [{cfg.min_size}..{cfg.max_size}]: "
          f"{len(gsd)} of {len(detected_size)} annotated sets")

    # ---- primary: gene permutation ------------------------------------
    print(f"[1/2] Gene permutation (prerank) — PRIMARY ...")
    t0 = time.time()
    pre = gp.prerank(rnk=welch_t, gene_sets=gsd, permutation_num=cfg.n_perm,
                     min_size=cfg.min_size, max_size=cfg.max_size,
                     outdir=None, seed=cfg.seed, verbose=False,
                     threads=cfg.threads)
    print(f"      {time.time() - t0:.1f} s")
    p2 = pre.res2d.copy()
    res = pd.DataFrame({
        "contrast": contrast,
        "term": p2["Term"].values,
        "NES": _num(p2["NES"]),
        "fdr": _num(p2["FDR q-val"]),
        "pval": _num(p2["NOM p-val"]),
        "leading_edge": p2["Lead_genes"].values,
    })

    # ---- robustness: phenotype permutation ----------------------------
    n_a, n_b = diag["n_samples_a"], diag["n_samples_b"]
    n_unique_splits = comb(n_a + n_b, n_a) // 2
    pheno_diag = {}
    if cfg.run_phenotype:
        imputed, pheno_diag = build_imputed_matrix(expr, cfg)
        cls = ["A" if c == cond_a else "B"
               for c in adata.obs.loc[imputed.columns, group_col]]
        print(f"[2/2] Phenotype permutation — ROBUSTNESS flag ...")
        print(f"      {cfg.n_perm} permutations drawn WITH REPLACEMENT from the "
              f"{n_unique_splits} distinct {n_a}-vs-{n_b} relabellings that exist; "
              f"the null is saturated and the resolution ceiling is "
              f"{1 / n_unique_splits:.2e}")
        t0 = time.time()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ph = gp.gsea(data=imputed, gene_sets=gsd, cls=cls,
                         permutation_type="phenotype", permutation_num=cfg.n_perm,
                         method="t_test", min_size=cfg.min_size,
                         max_size=cfg.max_size, outdir=None, seed=cfg.seed,
                         verbose=False, threads=cfg.threads)
        print(f"      {time.time() - t0:.1f} s")
        h2 = ph.res2d.copy()
        res = res.merge(pd.DataFrame({
            "term": h2["Term"].values,
            "NES_pheno": _num(h2["NES"]),
            "fdr_pheno": _num(h2["FDR q-val"]),
        }), on="term", how="left")
    else:
        res["NES_pheno"] = np.nan
        res["fdr_pheno"] = np.nan

    res["pheno_robust"] = res["fdr_pheno"] < cfg.fdr_cut
    res["direction"] = np.where(res["NES"] > 0, "up", "down")
    res["detected_size"] = [detected_size[t] for t in res["term"]]
    res = _reject_sentinel_terms(res, contrast, cfg)

    # ---- report --------------------------------------------------------
    n_sig = int((res["fdr"] < cfg.fdr_cut).sum())
    n_rob = int(res["pheno_robust"].sum())
    print(f"\nPRIMARY (gene permutation) FDR<{cfg.fdr_cut}: {n_sig} terms")
    print(f"ROBUST  (phenotype)        FDR<{cfg.fdr_cut}: {n_rob} terms "
          f"({100 * n_rob / max(n_sig, 1):.0f}% of primary hits)")
    for d in ("up", "down"):
        sub = res[res["direction"] == d]
        print(f"  {d:<4}: {len(sub):4d} terms | primary sig {int((sub['fdr'] < cfg.fdr_cut).sum()):4d}"
              f" | nominal p<0.05 {int((sub['pval'] < 0.05).sum()):4d}"
              f" | pheno-robust {int(sub['pheno_robust'].sum()):4d}")

    # ---- write ---------------------------------------------------------
    out = cfg.contrast_dir(contrast)
    raw_path = out / f"{_today()}_raw.csv"
    res.to_csv(raw_path, index=False)

    params = {
        "contrast": contrast, "cond_a": cond_a, "cond_b": cond_b,
        "group_col": group_col, "symbol_col": symbol_col,
        "date": _today(),
        "ranking": "Welch t (unequal variance), NaN-aware, no imputation",
        "primary_null": "gene permutation (gseapy.prerank)",
        "secondary_null": "phenotype permutation (gseapy.gsea)" if cfg.run_phenotype else None,
        "n_perm": cfg.n_perm, "seed": cfg.seed,
        "n_unique_relabellings": n_unique_splits,
        "phenotype_p_floor": 1 / n_unique_splits,
        "min_size": cfg.min_size, "max_size": cfg.max_size,
        "min_valid_per_group": cfg.min_valid_per_group,
        "overlap_thresh": cfg.overlap_thresh, "fdr_cut": cfg.fdr_cut,
        "collections": list(cfg.collections),
        "genesets_cache": str(cfg.genesets_cache),
        "n_gene_sets_tested": len(gsd),
        "n_artifact_pathways": len(cfg.artifact_pathways),
        **diag, **pheno_diag,
    }
    with open(out / "params.json", "w") as fh:
        json.dump(params, fh, indent=2)
    print(f"\nWrote {raw_path}\n      {out / 'params.json'}")

    res.attrs["welch_t"] = welch_t
    res.attrs["contrast"] = contrast
    return res


def run_contrast_from_rank(
    ranking: pd.Series,
    contrast: str,
    gene_sets: pd.DataFrame,
    cfg: GSEAConfig,
    ranking_description: str,
    extra_params: Optional[dict] = None,
) -> pd.DataFrame:
    """Gene-permutation GSEA on a ranking computed elsewhere.

    ``run_contrast`` builds its own ranking with :func:`welch_t_ranking`, a
    two-group Welch t over samples. That is wrong whenever the design is not two
    independent groups of samples: for a paired or donor-level contrast it pools
    replicate dissections from the same donor into one group and treats them as
    independent, which inflates significance. This entry point takes the ranking
    from the design that fits the data and runs only the prerank half.

    ``ranking`` must be indexed by gene symbol, unique, finite, and on a
    statistic scale rather than a p-value scale: ``prerank`` weights the
    enrichment score by |metric|, so -log10(p) would let a handful of genes
    dominate every set. A signed t or z serves; a signed -log10(p) does not.

    ``ranking_description`` is written verbatim into params.json. It is required
    because raw.csv on its own cannot say what was ranked, and a stored result
    whose statistic is unrecoverable cannot be reported.

    The phenotype null is not available here: it permutes sample labels and
    needs the expression matrix and a class vector, neither of which a ranking
    carries. ``cfg.run_phenotype`` must therefore be False, and params.json
    records the absence rather than leaving the field to be misread as untested.

    Returns the same columns and the same ``attrs`` as :func:`run_contrast`, so
    :func:`build_review` and Pass 2 take this result unchanged.
    """
    if cfg.run_phenotype:
        raise ValueError(
            "run_contrast_from_rank cannot run the phenotype null: it permutes "
            "sample labels and needs the expression matrix. Set "
            "GSEAConfig(run_phenotype=False) for this contrast."
        )
    if not isinstance(ranking, pd.Series):
        raise TypeError(f"ranking must be a Series, got {type(ranking).__name__}")
    if not ranking_description:
        raise ValueError("ranking_description is required and must say what was ranked")
    if ranking.index.has_duplicates:
        dup = ranking.index[ranking.index.duplicated()].unique()[:5].tolist()
        raise ValueError(f"ranking index has duplicate genes, e.g. {dup}. Collapse "
                         f"protein groups onto one representative per symbol first.")
    if not np.isfinite(ranking.values).all():
        raise ValueError(f"ranking has {int((~np.isfinite(ranking.values)).sum())} "
                         f"non-finite values")
    ranking = ranking.sort_values(ascending=False)

    detected = set(ranking.index)
    gsd_all = _gene_sets_dict(gene_sets)
    detected_size = {pw: len(set(g) & detected) for pw, g in gsd_all.items()}
    gsd = {pw: gsd_all[pw] for pw, n in detected_size.items()
           if cfg.min_size <= n <= cfg.max_size}
    print(f"Ranking: {len(ranking)} genes [{ranking.iloc[0]:+.2f} … "
          f"{ranking.iloc[-1]:+.2f}]  ({ranking_description})")
    print(f"Detected-size filter [{cfg.min_size}..{cfg.max_size}]: "
          f"{len(gsd)} of {len(detected_size)} annotated sets")
    if not gsd:
        raise ValueError("No gene set survives the detected-size filter. Lower "
                         "cfg.min_size or check that the ranking index holds gene "
                         "symbols matching the gene-set table.")

    print("[1/1] Gene permutation (prerank) — PRIMARY ...")
    t0 = time.time()
    pre = gp.prerank(rnk=ranking, gene_sets=gsd, permutation_num=cfg.n_perm,
                     min_size=cfg.min_size, max_size=cfg.max_size,
                     outdir=None, seed=cfg.seed, verbose=False,
                     threads=cfg.threads)
    print(f"      {time.time() - t0:.1f} s")

    p2 = pre.res2d.copy()
    res = pd.DataFrame({
        "contrast": contrast,
        "term": p2["Term"].values,
        "NES": _num(p2["NES"]),
        "fdr": _num(p2["FDR q-val"]),
        "pval": _num(p2["NOM p-val"]),
        "leading_edge": p2["Lead_genes"].values,
    })
    res["NES_pheno"] = np.nan
    res["fdr_pheno"] = np.nan
    res["pheno_robust"] = res["fdr_pheno"] < cfg.fdr_cut
    res["direction"] = np.where(res["NES"] > 0, "up", "down")
    res["detected_size"] = [detected_size[t] for t in res["term"]]
    res = _reject_sentinel_terms(res, contrast, cfg)

    n_sig = int((res["fdr"] < cfg.fdr_cut).sum())
    print(f"\nPRIMARY (gene permutation) FDR<{cfg.fdr_cut}: {n_sig} terms")
    for d in ("up", "down"):
        sub = res[res["direction"] == d]
        print(f"  {d:<4}: {len(sub):4d} terms | primary sig "
              f"{int((sub['fdr'] < cfg.fdr_cut).sum()):4d}"
              f" | nominal p<0.05 {int((sub['pval'] < 0.05).sum()):4d}")

    out = cfg.contrast_dir(contrast)
    raw_path = out / f"{_today()}_raw.csv"
    res.to_csv(raw_path, index=False)

    params = {
        "contrast": contrast, "date": _today(),
        "ranking": ranking_description,
        "ranking_source": "supplied to run_contrast_from_rank",
        "primary_null": "gene permutation (gseapy.prerank)",
        "secondary_null": None,
        "secondary_null_absent_because":
            "the phenotype null needs the expression matrix and a class vector; "
            "a ranking carries neither",
        "n_perm": cfg.n_perm, "seed": cfg.seed,
        "min_size": cfg.min_size, "max_size": cfg.max_size,
        "overlap_thresh": cfg.overlap_thresh, "fdr_cut": cfg.fdr_cut,
        "collections": list(cfg.collections),
        "genesets_cache": str(cfg.genesets_cache),
        "n_genes_primary": int(len(ranking)),
        "n_gene_sets_tested": len(gsd),
        "n_artifact_pathways": len(cfg.artifact_pathways),
        **(extra_params or {}),
    }
    with open(out / "params.json", "w") as fh:
        json.dump(params, fh, indent=2)
    print(f"\nWrote {raw_path}\n      {out / 'params.json'}")

    res.attrs["welch_t"] = ranking
    res.attrs["contrast"] = contrast
    return res


def _reject_sentinel_terms(res: pd.DataFrame, contrast: str,
                           cfg: GSEAConfig) -> pd.DataFrame:
    """Remove or refuse terms that returned the NES=1.0 / p=1.0 sentinel.

    NES is ``es / mean(esnull on the sign of es)`` and the nominal p is
    ``#(esnull >= es) / #(esnull >= 0)``. Both can equal exactly 1.0 only if the
    permutation null is CONSTANT and equal to the observed score, i.e. the
    permutation did nothing for that term. The value is not a statistic, and
    because every affected term shares one NES it also shares one FDR.

    This runs before anything is written, cached, clustered or reviewed: a
    degenerate term can be picked as a cluster representative, and the Pass 2
    replay assertion would then faithfully reproduce a decision made on an
    invalid value.

    Terms already on ``cfg.artifact_pathways`` are reported but not acted on,
    because ``build_review`` drops them before clustering, so they cannot reach
    a figure.
    """
    nes = pd.to_numeric(res["NES"], errors="coerce")
    pv = pd.to_numeric(res["pval"], errors="coerce")
    flagged = (np.isclose(pv, 1.0, rtol=0, atol=1e-12)
               & np.isclose(nes, 1.0, rtol=0, atol=1e-9))
    excluded = res["term"].isin(cfg.artifact_pathways).to_numpy()
    n_excl = int((flagged & excluded).sum())
    if n_excl:
        names = sorted(res.loc[flagged & excluded, "term"])
        print(f"Sentinel check: {n_excl} affected term(s) are on the artifact list "
              f"and are dropped before clustering anyway: {names[:5]}")
    sentinel = flagged & ~excluded
    n = int(sentinel.sum())
    if n == 0:
        print(f"Sentinel check: no usable term of {len(res)} returned "
              f"NES=1.0 / p=1.0")
        return res

    bad = res.loc[sentinel]
    fdrs = sorted(pd.to_numeric(bad["fdr"], errors="coerce").round(6).unique())
    sizes = pd.to_numeric(bad["detected_size"], errors="coerce")
    allsz = pd.to_numeric(res["detected_size"], errors="coerce")
    msg = (f"{n} of {len(res)} terms returned the NES=1.0 / p=1.0 sentinel in "
           f"contrast {contrast!r}. The permutation null was constant for these "
           f"terms, so NES, p and FDR are all invalid.\n"
           f"  shared FDR value(s): {fdrs}\n"
           f"  detected size range: {int(sizes.min())}-{int(sizes.max())} "
           f"(all terms: {int(allsz.min())}-{int(allsz.max())})\n"
           f"  first five: {sorted(bad['term'])[:5]}")
    if cfg.sentinel_action == "raise":
        raise ValueError(msg + "\n  Set GSEAConfig(sentinel_action='drop') to "
                               "continue with these terms removed. Pass 1 must "
                               "then be regenerated and every MOVED row "
                               "re-reviewed, because clustering ran on these "
                               "values.")
    if cfg.sentinel_action != "drop":
        raise ValueError(f"sentinel_action must be 'raise' or 'drop', got "
                         f"{cfg.sentinel_action!r}")
    print("WARNING: " + msg)
    print(f"  Dropping them. {len(res) - n} terms remain. Regenerate Pass 1 and "
          f"re-review any MOVED row.")
    return res.loc[~sentinel].reset_index(drop=True)


# ══════════════════════════════════════════════════════════════════════════
# Pass 1 — clustering and the decision file
# ══════════════════════════════════════════════════════════════════════════

def _cluster_one_direction(les: List[set], thresh: float) -> np.ndarray:
    """Complete-linkage clustering on 1 - overlap coefficient.

    Complete linkage, not single: single chains catastrophically here (377 of
    402 up-terms merged into one core-less cluster, and raising the threshold
    does not help because containment already scores 1.0). With complete
    linkage, cutting at ``1 - t`` is exactly equivalent to requiring every
    pairwise overlap within a cluster to exceed ``t``.
    """
    n = len(les)
    if n == 0:
        return np.array([], dtype=int)
    if n == 1:
        return np.array([1])
    sim = np.eye(n)
    for i, j in combinations(range(n), 2):
        sim[i, j] = sim[j, i] = overlap_coefficient(les[i], les[j])
    Z = linkage(squareform(1 - sim, checks=False), method="complete")
    return fcluster(Z, t=1 - thresh, criterion="distance")


def build_review(
    res: pd.DataFrame,
    cfg: GSEAConfig,
    welch_t: Optional[pd.Series] = None,
    contrast: Optional[str] = None,
) -> pd.DataFrame:
    """Pass 1: cluster redundant terms and write the decision file.

    Drops nothing from the analysis — it only *proposes* collapses for review.
    Writes ``review.csv`` (persistent, holds decisions) and a dated
    ``_terms_full.csv`` (supplementary, full leading edges).
    """
    welch_t = welch_t if welch_t is not None else res.attrs.get("welch_t")
    contrast = contrast or res.attrs.get("contrast")
    if welch_t is None or contrast is None:
        raise ValueError("Need welch_t and contrast (pass explicitly, or use the "
                         "DataFrame returned by run_contrast)")

    rank_lookup = {g: i for i, g in enumerate(welch_t.index)}
    labeller = FamilyLabeller.from_yaml(cfg.families_path)

    cand = res[res["fdr"] < cfg.fdr_cut].copy()
    cand = cand[~cand["term"].isin(cfg.artifact_pathways)]
    cand["le"] = cand["leading_edge"].map(parse_leading_edge)

    # A term with an empty leading edge overlaps no other term, so it would stand
    # as a cluster of its own and become its representative.
    n_empty = int((cand["le"].map(len) == 0).sum())
    if n_empty:
        print(f"WARNING: dropping {n_empty} candidate(s) with an empty leading edge:")
        for t in cand.loc[cand["le"].map(len) == 0, "term"]:
            print(f"    {t}")
    cand = cand[cand["le"].map(len) > 0].reset_index(drop=True)
    assert int((cand["le"].map(len) == 0).sum()) == 0, "empty leading edge survived"
    print(f"Candidates at FDR<{cfg.fdr_cut}: {len(cand)} terms "
          f"({(cand['direction'] == 'up').sum()} up, "
          f"{(cand['direction'] == 'down').sum()} down)")

    # ---- cluster within direction, with GLOBALLY unique ids -------------
    cand["cluster_key"] = ""
    for direction in ("up", "down"):
        m = cand["direction"] == direction
        idx = cand.index[m].tolist()
        if not idx:
            continue
        lab = _cluster_one_direction(cand.loc[m, "le"].tolist(), cfg.overlap_thresh)
        cand.loc[idx, "cluster_key"] = [f"{direction}_{v}" for v in lab]
        big = pd.Series(lab).value_counts().iloc[0]
        print(f"  {direction}: {len(set(lab))} clusters, "
              f"largest {big} ({100 * big / len(idx):.1f}%)")

    # ---- per-cluster statistics -----------------------------------------
    info: Dict[str, dict] = {}
    for key, grp in cand.groupby("cluster_key"):
        members = grp["term"].tolist()
        direction = grp.iloc[0]["direction"]
        les = grp["le"].tolist()
        core = set.intersection(*les)
        counts: Dict[str, int] = {}
        for s in les:
            for g in s:
                counts[g] = counts.get(g, 0) + 1
        consensus = {g for g, c in counts.items()
                     if c >= np.ceil(len(members) / 2)}
        min_ov = (min(overlap_coefficient(a, b) for a, b in combinations(les, 2))
                  if len(les) > 1 else np.nan)
        core_o = order_genes_by_rank(core, direction, rank_lookup)
        cons_o = order_genes_by_rank(consensus, direction, rank_lookup)

        # Representative: |NES| desc -> FDR asc -> smaller LE -> name.
        # |NES| leads so the plotted dot IS the value that earned the row its
        # slot; best-FDR-first could plot a weak member of a strong cluster.
        g2 = grp.assign(_abs=grp["NES"].abs(), _le=grp["le"].map(len))
        g2 = g2.sort_values(["_abs", "fdr", "_le", "term"],
                            ascending=[False, True, True, True])
        info[key] = dict(
            members=members, direction=direction, rep=g2.iloc[0]["term"],
            core=core_o, consensus=cons_o, core_size=len(core), min_ov=min_ov,
            nes_min=float(grp["NES"].min()), nes_max=float(grp["NES"].max()),
            signed_max=float(grp["NES"].max() if direction == "up" else grp["NES"].min()),
            family=labeller.label_cluster(core_o, cons_o),
            hash=hashlib.md5(";".join(sorted(members)).encode()).hexdigest()[:8],
        )

    for direction in ("up", "down"):
        keys = [k for k, d in info.items() if d["direction"] == direction]
        keys.sort(key=lambda k: -abs(info[k]["signed_max"]))
        for rank, k in enumerate(keys, start=1):
            info[k]["plot_rank"] = rank

    # Globally unique integer ids: Pass 2 keys exclusions and de-duplication on
    # cluster_id alone, so an id shared by an up and a down cluster deletes rows.
    cl_num = {k: i + 1 for i, k in enumerate(sorted(info))}
    term_to_key = {m: k for k, d in info.items() for m in d["members"]}

    rows = []
    for _, r in cand.iterrows():
        key = term_to_key[r["term"]]
        d = info[key]
        high = d["plot_rank"] <= cfg.n_review_clusters
        rows.append({
            "term": r["term"], "direction": r["direction"],
            "NES": round(float(r["NES"]), 3), "fdr": float(r["fdr"]),
            "pval": float(r["pval"]),
            "fdr_pheno": float(r["fdr_pheno"]) if pd.notna(r["fdr_pheno"]) else np.nan,
            "pheno_robust": bool(r["pheno_robust"]),
            "cluster_id": cl_num[key], "cluster_size": len(d["members"]),
            "cluster_plot_rank": d["plot_rank"],
            "review_priority": "HIGH" if high else "low",
            "family_label": d["family"],
            "detected_size": int(r["detected_size"]), "le_size": len(r["le"]),
            "core_size": d["core_size"],
            "min_pairwise_overlap": (round(d["min_ov"], 3)
                                     if pd.notna(d["min_ov"]) else ""),
            "nes_min": round(d["nes_min"], 3), "nes_max": round(d["nes_max"], 3),
            "le_genes_top": ";".join(order_genes_by_rank(
                r["le"], r["direction"], rank_lookup, cfg.top_k_genes)) if high else "",
            "core_genes": ";".join(d["core"][:cfg.top_k_genes]) if high else "",
            "consensus_genes": ";".join(d["consensus"][:cfg.top_k_genes]) if high else "",
            "cluster_members_hash": d["hash"],
            "is_suggested_rep": r["term"] == d["rep"],
            "display_choice": "rep" if r["term"] == d["rep"] else "",
            "display_label": "",
            "_full_le": ";".join(order_genes_by_rank(
                r["le"], r["direction"], rank_lookup)),
        })

    review = pd.DataFrame(rows)
    review = _sort_review(review)
    review["term_rank"] = review.groupby("direction").cumcount() + 1

    out = cfg.contrast_dir(contrast)
    audit = (review.drop(columns=["display_choice", "display_label", "le_genes_top",
                                  "core_genes", "consensus_genes"])
                   .rename(columns={"_full_le": "leading_edge_full"}))
    audit_path = out / f"{_today()}_terms_full.csv"
    audit.to_csv(audit_path, index=False)
    review = review.drop(columns="_full_le")

    review = _merge_prior_decisions(review, out / "review.csv")
    review = _sort_review(review)
    review_path = out / "review.csv"
    review.to_csv(review_path, index=False)

    n_high = int((review["review_priority"] == "HIGH").sum())
    print(f"\nAudit table : {audit_path}")
    print(f"Decision file: {review_path}")
    print(f"  {len(review)} terms, {review['cluster_id'].nunique()} clusters, "
          f"{n_high} HIGH-priority rows (top {cfg.n_review_clusters} clusters "
          f"per direction)")
    print("  -> hand ONLY the HIGH rows to the reviewer; low rows carry no gene lists")
    return review


def _sort_review(df: pd.DataFrame) -> pd.DataFrame:
    """Sort by term strength within direction: up block first, then down."""
    o = np.where(df["direction"] == "up", -df["NES"], df["NES"])
    return df.assign(_o=o).sort_values(["direction", "_o"]).drop(columns="_o")


def _merge_prior_decisions(review: pd.DataFrame, path: Path) -> pd.DataFrame:
    """Carry forward decisions, keyed by TERM NAME — never by cluster number.

    Cluster ids renumber whenever the data changes. A term whose *cluster
    membership* changed is flagged MOVED rather than silently inheriting a
    decision that may no longer mean what was intended.
    """
    if not Path(path).exists():
        review["review_status"] = "NEW"
        print(f"No prior review file — all {len(review)} terms marked NEW")
        return review

    prior = pd.read_csv(path, keep_default_na=False)
    cols = [c for c in ["term", "display_choice", "display_label",
                        "cluster_members_hash"] if c in prior.columns]
    prior = prior[cols].rename(columns={"display_choice": "_pc",
                                        "display_label": "_pl",
                                        "cluster_members_hash": "_ph"})
    review = review.merge(prior, on="term", how="left")
    seen = review["_pc"].notna()
    moved = seen & (review["_ph"] != review["cluster_members_hash"])
    review["review_status"] = np.where(~seen, "NEW", np.where(moved, "MOVED", "ok"))
    take = seen & ~moved
    review.loc[take, "display_choice"] = review.loc[take, "_pc"]
    review.loc[take, "display_label"] = review.loc[take, "_pl"].fillna("")
    review = review.drop(columns=["_pc", "_pl", "_ph"])
    print(f"Merged prior review: {int((review['review_status'] == 'ok').sum())} ok, "
          f"{int((review['review_status'] == 'NEW').sum())} NEW, "
          f"{int((review['review_status'] == 'MOVED').sum())} MOVED")
    return review


# ══════════════════════════════════════════════════════════════════════════
# Guardrail — replay Pass 2's selection and fail if a decision is lost
# ══════════════════════════════════════════════════════════════════════════

def assert_pass2_selection(review: pd.DataFrame) -> pd.DataFrame:
    """Replay exactly what Pass 2 does, and fail if any representative is lost.

    Checks that no ``cluster_id`` is shared by the two directions, and that every
    declared representative outside an excluded cluster survives the selection.
    """
    if review["cluster_id"].dtype == object:
        review = review.copy()
        review["cluster_id"] = pd.to_numeric(review["cluster_id"], errors="coerce")

    bad = review.groupby("cluster_id")["direction"].nunique()
    if not bad.eq(1).all():
        offenders = list(bad[bad > 1].index)[:5]
        raise AssertionError(
            f"cluster_id reused across directions: {offenders}. Pass 2 keys "
            f"exclusions and de-duplication on cluster_id alone, so this "
            f"silently deletes rows.")

    choice = review["display_choice"].astype(str).str.strip().str.lower()
    excluded = set(review.loc[choice == "exclude", "cluster_id"])
    cluster_of = dict(zip(review["term"], review["cluster_id"]))

    # Representatives the reviewer declared, minus those in deliberately
    # excluded clusters — every one of these must survive Pass 2's selection.
    expected = {t for t in review.loc[choice == "rep", "term"]
                if cluster_of[t] not in excluded}

    surviving = review[~review["cluster_id"].isin(excluded)]
    reps = surviving[surviving["display_choice"].astype(str).str.strip().str.lower() == "rep"]
    reps = (reps.assign(_abs=reps["NES"].abs())
                .sort_values("_abs", ascending=False)
                .drop_duplicates("cluster_id", keep="first"))

    lost = expected - set(reps["term"])
    if lost:
        raise AssertionError(
            f"{len(lost)} representative(s) would be dropped by Pass 2: "
            f"{sorted(lost)[:5]}")
    print(f"Pass 2 replay check: all {len(reps)} representatives survive selection")
    return reps


def report_undecided(review: pd.DataFrame, name: str = "") -> pd.DataFrame:
    """Print and return the HIGH-priority clusters that carry no decision.

    ``review_status`` records only changes of cluster membership, and Pass 1
    marks every cluster ``display_choice = "rep"``, so an unreviewed cluster also
    reads "ok". A representative without a ``display_label`` is the mark of a
    cluster nobody decided on; it would be plotted under its raw term name.
    """
    choice = review["display_choice"].astype(str).str.strip().str.lower()
    label = review["display_label"].astype(str).str.strip()
    undecided = review[(review["review_priority"] == "HIGH") & (choice == "rep")
                       & (label == "")]
    if len(undecided):
        prefix = f"{name}: " if name else ""
        print(f"!! {prefix}{undecided['cluster_id'].nunique()} HIGH cluster(s) are "
              f"'rep' with no display_label (not reviewed):")
        for _, r in undecided.iterrows():
            print(f"     cid={r['cluster_id']:<4} NES={float(r['NES']):+.2f}  "
                  f"{r['term'][:60]}")
    return undecided


LABEL_GENES_PATH = _MODULE_DIR / "gsea_label_genes.yaml"


def load_label_genes(path: Optional[Path] = None) -> Dict[str, set]:
    """Read the gene list that defines each row label.

    A label names a module; the list says which genes that name means. Terms are
    scored against it, never against the cluster's own core genes, which are the
    intersection of the same leading edges and therefore cannot tell a
    mislabelled cluster from a correct one.
    """
    p = Path(path or LABEL_GENES_PATH)
    if not p.exists():
        raise FileNotFoundError(
            f"{p} is missing. Every label drawn needs a gene list; see the "
            f"header of that file for how they are written.")
    with open(p) as fh:
        raw = yaml.safe_load(fh)
    return {k: {str(g).strip() for g in v if str(g).strip()} for k, v in raw.items()}


def score_terms(raw: pd.DataFrame, genes: Iterable[str]) -> pd.DataFrame:
    """Purity and coverage of every term's leading edge against a gene list.

    purity   = |leading edge ∩ list| / |leading edge|   — is this term measuring
               the module, or something else that happens to share a few genes
    coverage = |leading edge ∩ list| / |list|           — how much of the module
               the term actually moved

    Purity is the discriminating measure. Coverage is a floor only: a large list
    is rarely covered above 0.5 by any single set.
    """
    genes = set(genes)
    idx = raw.index if raw.index.name == "term" else raw["term"]
    le = raw["leading_edge"] if "leading_edge" in raw.columns else raw["le_genes_top"]
    rows = []
    for term, s in zip(idx, le):
        e = {g.strip() for g in str(s).split(";") if g.strip()}
        n = len(e & genes)
        rows.append((term, len(e), n,
                     n / len(e) if e else np.nan,
                     n / len(genes) if genes else np.nan))
    return pd.DataFrame(rows, columns=["term", "le_size", "n_in_list",
                                       "purity", "coverage"]).set_index("term")


def report_row_redundancy(grid: pd.DataFrame, raws: Dict[str, pd.DataFrame],
                          contrasts: Sequence[str], thresh: float = 0.6) -> pd.DataFrame:
    """Row pairs whose drawn terms share most of a leading edge, per contrast.

    Two rows that share 0.95 of a leading edge are one measurement drawn twice.
    Run this AFTER the terms are chosen — re-choosing a term often removes the
    overlap, because the overlap came from two off-module terms sharing the same
    unrelated genes. Merge only pairs that repeat in EVERY column: a pair that
    agrees in one contrast and not another is carrying a real difference.
    """
    def _le(c, t):
        r = raws[c]
        if not isinstance(t, str) or t not in r.index:
            return set()
        return {g.strip() for g in str(r.loc[t, "leading_edge"]).split(";") if g.strip()}

    hits = []
    for c in contrasts:
        sets = {r["label"]: _le(c, r.get(f"term__{c}")) for _, r in grid.iterrows()}
        sets = {k: v for k, v in sets.items() if len(v) >= 10}
        for a, b in combinations(sorted(sets), 2):
            ov = len(sets[a] & sets[b]) / min(len(sets[a]), len(sets[b]))
            if ov >= thresh:
                hits.append({"contrast": c, "row_a": a, "row_b": b, "overlap": round(ov, 2)})
    out = pd.DataFrame(hits)
    if out.empty:
        print(f"Row redundancy: no pair shares >= {thresh:.0%} of a leading edge")
        return out
    every = (out.groupby(["row_a", "row_b"])["contrast"].nunique() == len(contrasts))
    print(f"Row redundancy: {len(out)} pair-contrast hit(s) at >= {thresh:.0%}")
    for (a, b), n in out.groupby(["row_a", "row_b"]).size().sort_values(ascending=False).items():
        mark = "  << in EVERY column: merge" if every.get((a, b), False) else ""
        print(f"  {n}/{len(contrasts)} columns  {a[:34]:34s} vs {b[:34]:34s}{mark}")
    return out


def _le_purity(raw: pd.DataFrame, term: str, genes: set) -> float:
    """Fraction of one term's leading edge that lies in a label's gene list."""
    if term not in raw.index:
        return 0.0
    col = "leading_edge" if "leading_edge" in raw.columns else "le_genes_top"
    e = {g.strip() for g in str(raw.loc[term, col]).split(";") if g.strip()}
    return len(e & set(genes)) / len(e) if e else 0.0


def report_offlist_cells(grid: pd.DataFrame, raws: Dict[str, pd.DataFrame],
                         contrasts: Sequence[str], label_genes: Dict[str, set],
                         cfg: GSEAConfig, purity_floor: float = 0.45) -> pd.DataFrame:
    """Hollow cells where some OTHER set for that module is significant.

    The row shows one gene set in every column. A cell can therefore read
    "tested, not significant" while a different set covering the same biology is
    significant in that contrast but fails the purity floor against the label's
    gene list. Every one of those is a question about the LIST, not about the
    figure: it is the label's genes that failed to recognise a real enrichment.

    Reported rather than drawn. Three cell states are already the limit of what a
    reader decodes, and a fourth marker would say "the label's gene list is
    imperfect here", which belongs in the notes.
    """
    hits = []
    for _, r in grid.iterrows():
        lab = r["label"]
        genes = label_genes.get(lab)
        if not genes:
            continue
        for c in contrasts:
            fdr = r.get(f"fdr__{c}")
            if pd.isna(fdr) or fdr < cfg.fdr_cut:
                continue                      # absent, or already significant
            raw = raws[c]
            # The caller passes the full raw table, so artifact terms are
            # removed here too; otherwise one can be named as the other set.
            sig = raw[(pd.to_numeric(raw["fdr"], errors="coerce") < cfg.fdr_cut)
                      & ~raw.index.isin(cfg.artifact_pathways)]
            if sig.empty:
                continue
            sc = score_terms(sig, genes)
            best = sc.sort_values("purity", ascending=False).head(1)
            if len(best) and float(best["purity"].iloc[0]) >= purity_floor * 0.6:
                t = best.index[0]
                hits.append({"label": lab, "contrast": c,
                             "row_term_NES": round(float(r.get(f"nes__{c}", np.nan)), 2),
                             "row_term_FDR": float(fdr),
                             "other_term": t,
                             "other_NES": round(float(raw.loc[t, "NES"]), 2),
                             "other_FDR": float(raw.loc[t, "fdr"]),
                             "other_purity": round(float(best["purity"].iloc[0]), 2)})
    out = pd.DataFrame(hits)
    if out.empty:
        print("Off-list cells: none — no hollow cell is covered by a significant "
              "set that the label's gene list failed to recognise")
        return out
    print(f"Off-list cells: {len(out)} hollow cell(s) where another significant "
          f"set covers the module")
    for _, h in out.iterrows():
        print(f"  {h['label'][:34]:34s} {h['contrast'][:22]:24s} "
              f"row {h['row_term_NES']:+.2f} (FDR {h['row_term_FDR']:.2g})  vs  "
              f"{h['other_term'][:40]} {h['other_NES']:+.2f} "
              f"(FDR {h['other_FDR']:.2g}, purity {h['other_purity']:.2f})")
    return out


def report_review_coverage(contrasts: Sequence[str], cfg: GSEAConfig,
                           min_labelled: int = 3, strict: bool = False) -> pd.DataFrame:
    """How much of each contrast's significant biology carries a label.

    Only a labelled cluster can become a row, so a contrast reviewed less than
    the others contributes fewer rows however strong its terms are — and the
    figure then reads as "less happens here", which is a statement about the
    review, not about the experiment.

    Returns one row per contrast and direction. With ``strict`` it raises
    instead of printing when a cell falls below ``min_labelled``.
    """
    rows = []
    for c in contrasts:
        rv = pd.read_csv(cfg.contrast_dir(c) / "review.csv", keep_default_na=False)
        rv["fdr"] = pd.to_numeric(rv["fdr"], errors="coerce")
        choice = rv["display_choice"].astype(str).str.strip().str.lower()
        label = rv["display_label"].astype(str).str.strip()
        hi = rv[rv["review_priority"] == "HIGH"]
        for d in ("up", "down"):
            blk = hi[hi["direction"] == d]
            sig_cl = set(blk.loc[blk["fdr"] < cfg.fdr_cut, "cluster_id"])
            lab_cl = set(rv.loc[(choice == "rep") & (label != ""), "cluster_id"])
            rows.append({"contrast": c, "direction": d,
                         "sig_HIGH_clusters": len(sig_cl),
                         "labelled": len(sig_cl & lab_cl),
                         "unlabelled": len(sig_cl - lab_cl)})
    cov = pd.DataFrame(rows)
    print("Review coverage (HIGH clusters significant at FDR < "
          f"{cfg.fdr_cut}):")
    for _, r in cov.iterrows():
        flag = "  << below min_labelled" if r["labelled"] < min_labelled else ""
        print(f"  {r['contrast']:<28s} {r['direction']:<4s} "
              f"{r['labelled']:3d} labelled of {r['sig_HIGH_clusters']:3d}{flag}")
    short = cov[cov["labelled"] < min_labelled]
    if len(short):
        msg = (f"{len(short)} contrast/direction block(s) carry fewer than "
               f"{min_labelled} labelled significant clusters: "
               f"{[f'{r.contrast} {r.direction}' for r in short.itertuples()]}")
        if strict:
            raise AssertionError(msg)
        print(f"  !! {msg}")
    return cov


# ══════════════════════════════════════════════════════════════════════════
# Plotting
# ══════════════════════════════════════════════════════════════════════════

def _clean_term(s: str) -> str:
    for p in ("REACTOME_", "GOBP_", "GOCC_", "GOMF_"):
        s = s.replace(p, "")
    return s.replace("_", " ").title()


def _latest(folder: Path, suffix: str) -> Path:
    hits = sorted(Path(folder).glob(f"*{suffix}"))
    if not hits:
        raise FileNotFoundError(f"No *{suffix} in {folder}")
    return hits[-1]                      # YYYYMMDD prefix sorts chronologically


def _size_to_area(v, vmin, vmax, lo=40, hi=260):
    return lo + (v - vmin) / max(vmax - vmin, 1) * (hi - lo)


def plot_contrast(contrast: str, cfg: GSEAConfig, save: bool = True,
                  figsize: Optional[Tuple[float, float]] = None,
                  title: Optional[str] = None):
    """Pass 2, single contrast. Rows are clusters; one dot each.

    No significance legend: candidates were filtered at FDR < ``fdr_cut`` before
    clustering, so every row is significant by construction. State that in the
    caption. The multi-contrast figure DOES need a legend, because there row
    admission and cell drawing follow different rules — the asymmetry is
    intentional.
    """
    out = cfg.contrast_dir(contrast)
    rev = pd.read_csv(out / "review.csv", keep_default_na=False)
    for c in ("NES", "fdr", "cluster_id", "cluster_size", "detected_size"):
        rev[c] = pd.to_numeric(rev[c], errors="coerce")

    unrev = rev[(rev["review_status"].isin(["NEW", "MOVED"]))
                & (rev["review_priority"] == "HIGH")]
    if len(unrev):
        print(f"!! {len(unrev)} HIGH-priority term(s) UNREVIEWED — using the "
              f"auto-suggestion:")
        for _, r in unrev.iterrows():
            print(f"   [{r['review_status']:<5}] {r['term'][:70]}")

    reps = assert_pass2_selection(rev)
    choice = rev["display_choice"].astype(str).str.strip().str.lower()
    excluded = set(rev.loc[choice == "exclude", "cluster_id"])
    if excluded:
        print(f"Excluded {len(excluded)} cluster(s) by review")

    missing = set(rev["cluster_id"]) - set(reps["cluster_id"]) - excluded
    if missing:
        print(f"!! {len(missing)} cluster(s) have no representative and are omitted")

    def _label(r):
        base = str(r["display_label"]).strip() or _clean_term(r["term"])
        n = int(r["cluster_size"])
        if n <= 1:
            return base
        return f"{base} ({n})" if cfg.short_collapse_label else f"{base}  ({n} terms collapsed)"

    reps = reps.copy()
    reps["label"] = reps.apply(_label, axis=1)
    sel = pd.concat([
        reps[reps["NES"] > 0].nlargest(cfg.top_n_plot, "NES"),
        reps[reps["NES"] < 0].nsmallest(cfg.top_n_plot, "NES"),
    ]).sort_values("NES").reset_index(drop=True)

    n_up = int((sel["NES"] > 0).sum())
    n_dn = int((sel["NES"] < 0).sum())
    print(f"Plotting {len(sel)} clusters ({n_up} up, {n_dn} down)")
    if n_up < cfg.top_n_plot or n_dn < cfg.top_n_plot:
        print(f"   NOTE: fewer than {cfg.top_n_plot} available in a direction after "
              f"exclusions — raise n_review_clusters for more headroom")

    vmax = np.ceil(sel["NES"].abs().max())
    norm = plt.Normalize(vmin=-vmax, vmax=vmax)
    smin, smax = sel["detected_size"].min(), sel["detected_size"].max()
    cmap = _cmap(cfg)

    # Width is budgeted explicitly: margin + label column + data area + legend
    # strip. tight_layout() does not reserve space for legends anchored outside
    # the axes, so the strip is held back with rect= further down. Every
    # constant scales with cfg's font sizes: fixed inches silently clip the keys
    # and the row labels once the fonts are raised to a notebook's own sizes.
    char_w = cfg.tick_fontsize / 72 * 0.62          # mean advance width, inches
    title_w = cfg.label_fontsize / 72 * 0.62
    lab_w = max(len(str(l)) for l in sel["label"]) * char_w
    LEGEND_IN = max(len("Detected set size") * title_w,
                    max(len(str(int(v))) for v in np.linspace(smin, smax, 3))
                    * char_w + 0.55) + 0.30
    strip = LEGEND_IN if cfg.keys_position == "right" else 0.0
    fig_w = 1.0 + lab_w + 3.0 + strip
    if figsize is None:
        figsize = (fig_w,
                   max(5, len(sel) * max(0.36, cfg.tick_fontsize / 72 * 1.6)))
    else:
        fig_w = figsize[0]
    fig, ax = plt.subplots(figsize=figsize)
    for i, r in sel.iterrows():
        ax.scatter(r["NES"], i, c=[r["NES"]], cmap=cmap, norm=norm,
                   s=_size_to_area(r["detected_size"], smin, smax,
                                   cfg.dot_area_min, cfg.dot_area_max),
                   edgecolors="black", linewidths=cfg.single_edge_lw, zorder=3)
    ax.set_yticks(range(len(sel)))
    ax.set_yticklabels(sel["label"], fontsize=cfg.tick_fontsize)
    ax.axvline(0, color="black", linewidth=cfg.single_zero_lw)
    ax.set_xlabel("NES", fontsize=cfg.label_fontsize)
    # The folder name is an identifier, not a caption. `title` overrides it;
    # pass "" for no title at all, which is what a composed figure wants.
    ax.set_title(f"GSEA — {contrast}" if title is None else title,
                 fontsize=cfg.title_fontsize, fontweight="bold")
    ax.grid(axis="x", alpha=0.2)
    ax.set_ylim(-0.7, len(sel) - 0.3)

    if cfg.nes_colorbar:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.02)
        cb.set_label("NES", fontsize=cfg.label_fontsize)
        cb.ax.tick_params(labelsize=cfg.tick_fontsize)

    fig.tight_layout(rect=[0, 0, 1 - strip / fig_w, 1])

    # Dots are drawn in points squared, which is independent of the canvas. On a
    # tall panel that leaves them as specks. Rescaling here rather than at
    # scatter() time is deliberate: the row pitch is only known once the layout
    # engine has run. The key is rescaled by the same factor, so it keeps
    # meaning what the dots mean.
    factor = 1.0
    if cfg.dot_fraction_of_row > 0 and len(sel) > 0:
        fig.canvas.draw()
        pitch_pt = ax.get_window_extent().height / fig.dpi * 72 / len(sel)
        factor = cfg.dot_fraction_of_row * pitch_pt / np.sqrt(cfg.dot_area_max)
        for coll in ax.collections:
            coll.set_sizes(np.asarray(coll.get_sizes()) * factor ** 2)
        print(f"Dots scaled x{factor:.2f}: largest {np.sqrt(cfg.dot_area_max) * factor:.0f} pt "
              f"across against a {pitch_pt:.0f} pt row pitch")

    # The key is built from the FINAL areas, after the rescale, so its markers
    # are the drawn ones. handlelength and handleheight are in font units, so a
    # marker larger than the font needs them raised or it prints over its label.
    legend_sizes = np.linspace(smin, smax, 3).astype(int)
    _ms = [np.sqrt(_size_to_area(v, smin, smax, cfg.dot_area_min, cfg.dot_area_max))
           * factor for v in legend_sizes]
    # Inside the axes the key has to be compact: handles just larger than the
    # largest marker, and tight spacing between the entries.
    inside = cfg.keys_position == "inside"
    if inside:
        _units = max(1.0, 1.05 * max(_ms) / cfg.tick_fontsize)
    else:
        _units = max(2.0, 1.25 * max(_ms) / cfg.tick_fontsize)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor="grey",
               markeredgecolor="black", markeredgewidth=cfg.key_edge_lw,
               markersize=m, label=str(v))
        for v, m in zip(legend_sizes, _ms)],
        title="Detected set size",
        loc=cfg.keys_inside_loc if inside else "upper left",
        bbox_to_anchor=None if inside else (1.02, 1.0),
        frameon=False, fontsize=cfg.tick_fontsize,
        title_fontsize=cfg.label_fontsize,
        handlelength=_units, handleheight=_units,
        labelspacing=0.3 if inside else 0.8)

    if save:
        stem = out / f"{_today()}_dotplot"
        fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
        fig.savefig(f"{stem}.svg", bbox_inches="tight")
        print(f"Saved {stem}.png / .svg")
    return fig, ax


def plot_combined(contrasts: Sequence[str], cfg: GSEAConfig,
                  contrast_labels: Optional[Sequence[str]] = None,
                  max_rows: int = 25, save: bool = True,
                  rows_up: Optional[int] = None,
                  rows_down: Optional[int] = None,
                  figsize: Optional[Tuple[float, float]] = None,
                  title: Optional[str] = None,
                  save_stem: str = "combined_dotplot",
                  cell_terms: str = "purity",
                  anchor_rule: str = "coverage",
                  min_rows_per_contrast: int = 3,
                  label_genes: Optional[Dict[str, set]] = None,
                  purity_floor: float = 0.45,
                  coverage_floor: float = 0.10,
                  cell_purity_floor: float = 0.30):
    """Multi-contrast figure. Rows are LABELS, not clusters.

    A term has a well-defined score in every contrast; a cluster does not,
    because leading edges are ranking-dependent. So each row is keyed by its
    reviewed label.

    ``cell_terms`` decides which gene set fills each cell.

    ``"purity"`` (default) picks one gene set per row from every term tested, by
    the purity of its leading edge against the label's gene list
    (gsea_label_genes.yaml), and reads that set in every column; see the
    comments of the purity selection below.

    ``"canonical"`` reads the row's anchor term in every column (Autophagy and
    CrossSystem). The columns then compare one gene set, but a cell can read
    non-significant where another set of the same module is significant.

    The anchor term is picked for COVERAGE — significant in the most contrasts,
    ties broken by max |NES| — not for peak strength, so the row table names a
    term the module supports in more than one column. ``anchor_rule="peak"``
    picks the term of largest |NES| instead.

    Row admission is "significant in >= 1 contrast", then a quota: each contrast
    is guaranteed ``min_rows_per_contrast`` rows per direction, ranked by ITS OWN
    NES, before the remaining slots are filled by max |NES| across contrasts.
    Without the quota the contrast with the widest separation buys nearly every
    row and a shallower contrast's biology cannot appear at all, however well it
    was reviewed.

    ``save_stem`` names the output files, so a second figure over a SUBSET of the
    contrasts does not overwrite the first. Both the label set and the terms are
    computed over the contrasts passed in, so a subset figure can legitimately
    show different values for the contrasts it keeps.
    """
    if cell_terms not in ("purity", "canonical"):
        raise ValueError(f"cell_terms must be 'purity' or 'canonical', "
                         f"got {cell_terms!r}")
    if anchor_rule not in ("coverage", "peak"):
        raise ValueError(f"anchor_rule must be 'coverage' or 'peak', "
                         f"got {anchor_rule!r}")
    contrast_labels = list(contrast_labels or contrasts)
    k = len(contrasts)

    raws, reviews = {}, {}
    for c in contrasts:
        d = cfg.contrast_dir(c)
        raws[c] = pd.read_csv(_latest(d, "_raw.csv"), keep_default_na=False)
        # Artifact terms are removed before clustering (build_review), but the
        # raw table still holds every tested term, and the purity selection
        # below searches all of them, so an artifact-listed term could otherwise
        # fill a row.
        raws[c] = raws[c][~raws[c]["term"].isin(cfg.artifact_pathways)].copy()
        raws[c]["NES"] = pd.to_numeric(raws[c]["NES"], errors="coerce")
        raws[c]["fdr"] = pd.to_numeric(raws[c]["fdr"], errors="coerce")
        rv = pd.read_csv(d / "review.csv", keep_default_na=False)
        rv["NES"] = pd.to_numeric(rv["NES"], errors="coerce")
        rv["fdr"] = pd.to_numeric(rv["fdr"], errors="coerce")
        rv["cluster_id"] = pd.to_numeric(rv["cluster_id"], errors="coerce")
        reviews[c] = rv

    # ---- label -> candidate terms, across all contrasts ------------------
    # cand[label] is the union across contrasts of the representatives that
    # carry the label; the anchor is drawn from it.
    cand: Dict[str, set] = {}
    _multi: Dict[Tuple[str, str], list] = {}
    for c in contrasts:
        rv = reviews[c]
        choice = rv["display_choice"].astype(str).str.strip().str.lower()
        excl = set(rv.loc[choice == "exclude", "cluster_id"])
        reps = rv[(choice == "rep") & (~rv["cluster_id"].isin(excl))]
        for _, r in reps.iterrows():
            lab = str(r["display_label"]).strip() or _clean_term(r["term"])
            cand.setdefault(lab, set()).add(r["term"])
            _multi.setdefault((c, lab), []).append((abs(float(r["NES"])), r["term"]))
    # One label may be carried by several clusters in the same contrast, which
    # is a legitimate reviewer decision (one module split across clusters). The
    # note names their |NES| range: a large spread means the label merged
    # clusters that do not belong together.
    for (c, lab), hits in _multi.items():
        hits.sort(reverse=True)
        if len(hits) > 1:
            print(f"  note: {c} has {len(hits)} representatives labelled "
                  f"{lab!r} (|NES| {hits[0][0]:.2f}..{hits[-1][0]:.2f})")
    if not cand:
        print("No reviewed representatives found in any contrast.")
        return None, None

    lookup = {c: raws[c].set_index("term") for c in contrasts}

    # ---- anchor term per label: COVERAGE first --------------------------
    # Significant in the most contrasts, then max |NES|, then lowest FDR, then
    # name. Ranking by peak |NES| alone hands the choice to whichever contrast
    # separates most widely, and that term is then read in every other column.
    def _stats(term):
        n_sig, best_abs, best_fdr = 0, -np.inf, np.inf
        for c in contrasts:
            if term in lookup[c].index:
                row = lookup[c].loc[term]
                if row["fdr"] < cfg.fdr_cut:
                    n_sig += 1
                if abs(row["NES"]) > best_abs:
                    best_abs, best_fdr = abs(row["NES"]), row["fdr"]
        return n_sig, best_abs, best_fdr

    anchor = {}
    for lab, terms in cand.items():
        scored = [(*_stats(t), t) for t in terms]
        if anchor_rule == "coverage":
            scored.sort(key=lambda x: (-x[0], -x[1], x[2], x[3]))
        else:                       # "peak": max |NES| only
            scored.sort(key=lambda x: (-x[1], x[2], x[3]))
        anchor[lab] = scored[0][3]

    # ---- purity selection ------------------------------------------------
    # Under cell_terms="purity" the term of a row is chosen from EVERY term
    # tested, not only the members of the cluster that received the label, and
    # it must clear the floors against the label's own gene list. Both matter:
    # the better term is often in a cluster that never received a representative
    # (DARPP-32, the TRiC term and HSF1 activation all were), and scoring against
    # the cluster's core cannot detect a mislabelled cluster.
    # Among the terms that pass the floors, the highest purity wins, then max
    # |NES|. Taking the largest |NES| first lands on sets that clear the floor by
    # a hair: for the proteasome row it chose a cell-cycle set at purity 0.36 over
    # REACTOME_PROTEASOME_ASSEMBLY at 0.80, and for the chaperonin row a
    # DNA-biosynthesis set at 0.37 over the TRiC term at 0.54.
    picked: Dict[Tuple[str, str], str] = {}
    dropped_labels: List[str] = []
    n_rowwide = 0
    if cell_terms == "purity":
        if label_genes is None:
            label_genes = load_label_genes()
        missing = [l for l in anchor if l not in label_genes]
        if missing:
            raise KeyError(f"{len(missing)} label(s) have no gene list: "
                           f"{sorted(missing)[:5]} — add them to "
                           f"{LABEL_GENES_PATH.name} before plotting")
        scores = {c: {} for c in contrasts}
        purity = {c: {} for c in contrasts}
        for c in contrasts:
            for lab in anchor:
                s = score_terms(raws[c], label_genes[lab])
                ok = s[(s["purity"] >= purity_floor) & (s["coverage"] >= coverage_floor)]
                scores[c][lab] = set(ok.index)
                purity[c][lab] = ok["purity"].to_dict()

        # A term is chosen ONCE per row, in the contrast where the module is
        # enriched, and then read in every column. Purity is a property of an
        # enrichment, not of a gene set: where a module does not move, the
        # leading edge is the top of a flat ranking and no term clears the
        # floor. Choosing per cell therefore left a hole in exactly the columns
        # where the answer is "tested, not significant", and drew it as though
        # nothing had been measured.
        for lab in anchor:
            # The row's term is fixed, its leading edge is not: a set can be
            # significant in one contrast through the module and in another
            # through different genes. A candidate is therefore DISQUALIFIED if
            # any contrast where it is significant — and so drawn filled — has a
            # leading edge below cell_purity_floor against the label's list.
            # That is a filter, not a ranking: ranking on the weakest purity
            # instead favours sets significant in a single column and lowers the
            # purity of the figure as a whole. Without the filter the actomyosin
            # row was carried by a cochlear hair-cell set scoring 0.00 and 0.10
            # in two of the three columns where it was significant.
            sig_cand = {}
            for t in set().union(*(scores[c][lab] for c in contrasts)):
                pur = [_le_purity(lookup[c], t, label_genes[lab])
                       for c in contrasts
                       if t in lookup[c].index
                       and float(lookup[c].loc[t, "fdr"]) < cfg.fdr_cut]
                if pur and min(pur) >= cell_purity_floor:
                    sig_cand[t] = (max(pur),
                                   max(abs(float(lookup[c].loc[t, "NES"]))
                                       for c in contrasts if t in lookup[c].index))
            if sig_cand:
                # sorted(): two Reactome sets with identical detected genes tie
                # exactly, and set order changes between processes; the
                # alphabetically first of a tie wins.
                best = max(sorted(sig_cand), key=lambda t: sig_cand[t])
                for c in contrasts:
                    picked[(c, lab)] = best
                n_rowwide += 1
                continue
            # Nothing significant anywhere clears the floors. Keep the purest
            # term that does, so the row can still be drawn hollow, and let the
            # admission rule below decide whether it earns a place at all.
            passing = [scores[c][lab] for c in contrasts]
            shared = set.union(*passing) if passing else set()
            if shared:
                def _pur(t):
                    return max(purity[c][lab].get(t, 0.0) for c in contrasts)
                best = max(sorted(shared), key=lambda t: (
                    _pur(t),
                    max(abs(float(lookup[c].loc[t, "NES"])) for c in contrasts
                        if t in lookup[c].index)))
                for c in contrasts:
                    picked[(c, lab)] = best
                continue
            dropped_labels.append(lab)
        print(f"  purity selection: floors {purity_floor:.2f} purity / "
              f"{coverage_floor:.2f} coverage; one gene set per row, chosen where "
              f"the module is enriched and read in every column; {n_rowwide} of "
              f"{len(anchor)} label(s) anchored on a significant term; "
              f"{len(dropped_labels)} label(s) have no qualifying term at all")
        for lab in sorted(dropped_labels):
            print(f"     dropped: {lab}")

    # ---- build the grid --------------------------------------------------
    # One term per row, read in every column: the purity pick, or the anchor
    # under "canonical". The term used is recorded per cell so the row table
    # states what was measured where.
    recs, n_absent_pool = [], 0
    for lab, a_term in anchor.items():
        if lab in dropped_labels:
            continue
        row = {"label": lab, "term": a_term}
        for c in contrasts:
            if cell_terms == "canonical":
                term = a_term
            else:
                term = picked.get((c, lab))
            row[f"term__{c}"] = term
            if term in lookup[c].index:
                e = lookup[c].loc[term]
                row[f"nes__{c}"] = float(e["NES"])
                row[f"fdr__{c}"] = float(e["fdr"])
                row[f"size__{c}"] = float(e["detected_size"])
            else:
                row[f"nes__{c}"] = np.nan
                row[f"fdr__{c}"] = np.nan
                row[f"size__{c}"] = np.nan
                n_absent_pool += 1
        recs.append(row)
    grid = pd.DataFrame(recs)
    grid["cell_terms"] = cell_terms

    fdr_cols = [f"fdr__{c}" for c in contrasts]
    nes_cols = [f"nes__{c}" for c in contrasts]
    grid["any_sig"] = (grid[fdr_cols] < cfg.fdr_cut).any(axis=1)
    grid["max_abs_nes"] = grid[nes_cols].abs().max(axis=1)
    # Admit on "significant in >=1 of k", trim by strength, then order the rows
    # by mean NES so down-programmes sit at the bottom and up at the top.
    grid = grid[grid["any_sig"]]
    if rows_up is None and rows_down is None:
        grid = grid.nlargest(max_rows, "max_abs_nes")
    else:
        # A single global trim by |NES| is direction-blind, and the two
        # directions do not compete on equal terms: up-NES runs higher, so the
        # down programme is squeezed out entirely long before it is exhausted.
        # Row direction = sign of the STRONGEST significant cell, the same rule
        # that picks the canonical term (max |NES| across contrasts), so a row
        # is filed under the direction of the measurement it is anchored on. A
        # non-significant cell does not vote. Summing the significant cells
        # instead lets several weak opposite cells outvote one strong one:
        # ER to Golgi (COPII) is +3.77 in HEK and mildly negative in three
        # shallower contrasts, and would be filed as a down row.
        signed = grid[nes_cols].where(grid[fdr_cols].lt(cfg.fdr_cut).to_numpy())
        strongest = signed.abs().idxmax(axis=1)
        is_up = pd.Series([signed.at[i, strongest[i]] >= 0 for i in grid.index],
                          index=grid.index)
        n_up = len(grid) if rows_up is None else rows_up
        n_dn = len(grid) if rows_down is None else rows_down

        def _pick(block, n, side):
            """Quota first, then the strongest of what is left.

            Each contrast is guaranteed up to ``min_rows_per_contrast`` rows,
            ranked by ITS OWN |NES| among the rows significant in it. A global
            ranking alone is decided by the contrast with the widest separation:
            in one two-contrast figure the down block cut at |NES| 2.21 while no
            module specific to the shallower contrast reached 2.09, so that
            column could not contribute a single row of its own.
            """
            if n >= len(block):
                return block
            quota, by_contrast = [], {}
            for c in contrasts:
                sig = block[block[f"fdr__{c}"] < cfg.fdr_cut]
                order = sig[f"nes__{c}"].abs().sort_values(ascending=False).index
                picked = [i for i in order[:min_rows_per_contrast]]
                by_contrast[c] = picked
                quota += [i for i in picked if i not in quota]
            rest = (block.drop(index=quota)
                         .nlargest(max(0, n - len(quota)), "max_abs_nes"))
            out = pd.concat([block.loc[quota], rest])
            print(f"  {side}: {len(block)} available, {len(quota)} held by the "
                  f"per-contrast quota ({min_rows_per_contrast}/contrast), "
                  f"{len(rest)} filled by max |NES|")
            for c, ids in by_contrast.items():
                own = [i for i in ids if i in out.index]
                if len(own) < min_rows_per_contrast:
                    print(f"    note: {c} contributes {len(own)} {side} row(s); "
                          f"it has fewer significant labels than the quota")
            if len(out) > n:
                print(f"    !! {len(out)} {side} rows exceeds the {n} requested: "
                      f"{len(contrasts)} contrasts x {min_rows_per_contrast} "
                      f"guaranteed. Lower min_rows_per_contrast or raise rows_"
                      f"{'up' if side == 'up' else 'down'}.")
            return out

        grid = pd.concat([_pick(grid[is_up], n_up, "up"),
                          _pick(grid[~is_up], n_dn, "down")])
        print(f"  per-direction trim: {int(is_up.sum())} up / "
              f"{int((~is_up).sum())} down available")
    grid = (grid.assign(_mean_nes=grid[nes_cols].mean(axis=1))
                .sort_values("_mean_nes")
                .drop(columns="_mean_nes")
                .reset_index(drop=True))

    print(f"Combined figure: {len(grid)} rows from {k} contrasts")
    print(f"  rows admitted on FDR<{cfg.fdr_cut} in >=1 of k={k} contrasts "
          f"(display rule, not an inferential claim — FDR is controlled WITHIN "
          f"each contrast)")
    # Counted on the final grid: only drawn cells are drawn as grey x.
    n_absent = int(grid[nes_cols].isna().to_numpy().sum())
    n_cells = len(grid) * k
    if n_absent:
        why = ("no term cleared the purity floor for that module in that contrast"
               if cell_terms == "purity"
               else "outside a contrast's detected-size window (never tested)")
        print(f"  {n_absent} of {n_cells} drawn cell(s) ({100 * n_absent / n_cells:.1f}%) "
              f"-> drawn as grey x: {why}")
        if n_absent > 0.1 * n_cells:
            print("  !! more than 10% of DRAWN cells; check the size window")
    if n_absent_pool:
        print(f"  (across all {len(anchor)} candidate labels before admission "
              f"and trimming: {n_absent_pool} untested cells)")

    # ---- draw -------------------------------------------------------------
    vmax = np.ceil(np.nanmax(np.abs(grid[nes_cols].values)))
    norm = plt.Normalize(vmin=-vmax, vmax=vmax)
    sizes = grid[[f"size__{c}" for c in contrasts]].values
    smin, smax = np.nanmin(sizes), np.nanmax(sizes)
    cmap = _cmap(cfg)

    # Width scales with the longest row label, else long labels squeeze the
    # axes and push the legends off the canvas. Width is budgeted explicitly:
    # margin + label column + data area + legend strip, as in plot_contrast.
    # tight_layout() does not reserve space for legends anchored outside the
    # axes, so the strip is held back with rect= further down. The allowances
    # scale with the font size; pass `figsize` to override, e.g. to fit a
    # composed-figure slot.
    # The key strings are defined once and the strip is measured from them, so
    # editing a label cannot leave the reservation stale.
    # The greater-or-equal sign is typeset, not ">=" — this is figure text.
    # Under purity selection the third state means "no gene set measured this
    # module in this contrast", which is not the same statement as "the term was
    # outside the size window", so the key says so.
    SIG_LABELS = (f"FDR < {cfg.fdr_cut}", f"FDR ≥ {cfg.fdr_cut}",
                  "not measured" if cell_terms == "purity" else "not tested")
    # A key entry for a mark the figure does not contain is noise, so the third
    # state is keyed only when at least one drawn cell is untested.
    if not n_absent:
        SIG_LABELS = SIG_LABELS[:2]
    KEY_TITLES = ("Detected set size", "Significance")
    char_w = cfg.tick_fontsize / 72 * 0.62           # mean advance width, inches
    title_w = cfg.label_fontsize / 72 * 0.62
    lab_w = max(len(str(l)) for l in grid["label"]) * char_w
    key_text = max(max(len(t) for t in KEY_TITLES) * title_w,
                   max(len(e) for e in SIG_LABELS) * char_w + 0.55)
    LEGEND_IN = key_text + 0.30                      # + gap to the colourbar
    strip = LEGEND_IN if cfg.keys_position != "corner" else 0.0
    fig_w = 1.0 + lab_w + k * 1.0 + strip
    if figsize is None:
        figsize = (fig_w,
                   max(5, len(grid) * max(0.36, cfg.tick_fontsize / 72 * 1.6)))
    else:
        fig_w = figsize[0]
    fig, ax = plt.subplots(figsize=figsize)
    for j, c in enumerate(contrasts):
        for i, r in grid.iterrows():
            nes, fdr, sz = r[f"nes__{c}"], r[f"fdr__{c}"], r[f"size__{c}"]
            if np.isnan(nes):                                   # never tested
                # The cross tracks the smallest dot, so raising the dot range
                # cannot leave the never-tested cells as invisible specks.
                ax.scatter(j, i, marker="x", s=cfg.dot_area_min, c="lightgrey",
                           linewidths=cfg.cross_lw, zorder=3)
            elif fdr < cfg.fdr_cut:                             # significant
                ax.scatter(j, i, c=[nes], cmap=cmap, norm=norm,
                           s=_size_to_area(sz, smin, smax,
                                           cfg.dot_area_min, cfg.dot_area_max),
                           edgecolors="black", linewidths=cfg.edge_lw_sig, zorder=3)
            else:                                               # tested, n.s.
                ax.scatter(j, i, facecolors="none", edgecolors="grey",
                           s=_size_to_area(sz, smin, smax,
                                           cfg.dot_area_min, cfg.dot_area_max),
                           linewidths=cfg.edge_lw_ns, zorder=3)

    ax.set_xticks(range(k))
    rot = cfg.contrast_label_rotation
    tick_kw = dict(rotation=rot, ha="right" if 0 < rot < 90 else "center")
    if rot == 90:
        tick_kw["multialignment"] = "right"        # every line ends at the axis
    ax.set_xticklabels(contrast_labels, fontsize=cfg.tick_fontsize, **tick_kw)
    ax.set_yticks(range(len(grid)))
    ax.set_yticklabels(grid["label"], fontsize=cfg.tick_fontsize)
    ax.set_xlim(-0.6, k - 0.4)
    ax.set_ylim(-0.7, len(grid) - 0.3)
    ax.grid(alpha=0.2, linewidth=cfg.grid_lw)
    # No title by default: panel titles belong to the composed figure, not to
    # the panel. Pass `title` to set one.
    if title:
        ax.set_title(title, fontsize=cfg.title_fontsize, fontweight="bold")

    corner = cfg.keys_position == "corner"
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.02)
    if corner:
        # Under the bar, horizontal: beside it, the label takes a column of width.
        cb.ax.set_xlabel("NES", fontsize=cfg.label_fontsize)
    else:
        cb.set_label("NES", fontsize=cfg.label_fontsize)
    cb.ax.tick_params(labelsize=cfg.tick_fontsize)

    sig_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="lightgrey",
               markeredgecolor="black", markersize=cfg.key_marker_size,
               markeredgewidth=cfg.key_mew_sig, label=SIG_LABELS[0]),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
               markeredgecolor="grey", markersize=cfg.key_marker_size,
               markeredgewidth=cfg.key_mew_ns, label=SIG_LABELS[1]),
    ]
    if len(SIG_LABELS) == 3:
        sig_handles.append(Line2D([0], [0], marker="x", color="lightgrey", lw=0,
                                  markersize=np.sqrt(cfg.dot_area_min),
                                  label=SIG_LABELS[2]))
    legend_sizes = np.linspace(smin, smax, 3).astype(int)
    size_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="grey",
               markeredgecolor="black", markeredgewidth=cfg.key_edge_lw,
               markersize=np.sqrt(_size_to_area(s, smin, smax, cfg.dot_area_min,
                                                cfg.dot_area_max)), label=str(s))
        for s in legend_sizes]
    key_style = dict(frameon=False, fontsize=cfg.tick_fontsize,
                     title_fontsize=cfg.label_fontsize)
    if corner:
        # One key, two columns, in the bottom-left corner: under the row labels
        # and left of the rotated contrast labels, which is empty. Its upper
        # right corner is at the lower left corner of the axes. Column 1 holds
        # the set sizes, column 2 the significance states; "genes" in the size
        # entries replaces the two key titles.
        for h, v in zip(size_handles, legend_sizes):
            h.set_label(f"{v} genes")
        ax.legend(handles=size_handles + sig_handles, ncol=2, loc="upper right",
                  bbox_to_anchor=(0, 0), **key_style)
        fig.tight_layout()
    else:
        sig_legend = ax.legend(handles=sig_handles, title="Significance",
                               loc="lower left", bbox_to_anchor=(1.02, 0.0),
                               **key_style)
        ax.add_artist(sig_legend)
        ax.legend(handles=size_handles, title="Detected set size",
                  loc="upper left", bbox_to_anchor=(1.02, 1.0), **key_style)
        fig.tight_layout(rect=[0, 0, 1 - LEGEND_IN / fig_w, 1])
    if save:
        stem = Path(cfg.output_root) / f"{_today()}_{save_stem}"
        fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
        fig.savefig(f"{stem}.svg", bbox_inches="tight")
        grid_path = Path(cfg.output_root) / f"{_today()}_{save_stem}_grid.csv"
        grid.to_csv(grid_path, index=False)
        print(f"Saved {stem}.png / .svg and {grid_path.name}")
    return fig, ax
