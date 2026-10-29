"""Dot plots of log2 fold change and FDR, genes x comparisons, shared by the aggregate DVP notebooks.

ttest_dotplot() draws one row per gene and one column per comparison. Colour is
the log2 fold change, dot area is -log10 FDR, a black outline marks FDR < 0.05,
and a grey cross marks a gene symbol that has no row in that comparison's table.
The genes come in named groups, whose names are drawn rotated left of the gene
names. Sizes follow figure_style: printed points and panel units.

tested_dotplot() draws the same encoding from one long table that records per
gene and comparison whether the gene was tested. The colour scale and the dot
area are capped at fixed values, the grey cross marks "not tested", and each
group of genes has a bracket next to its rotated name.

cross_system_dotplot() draws the same encoding from two gene x column matrices,
log2 fold change and FDR, which the caller has already reduced to one value per
gene. The colour scale and the dot area run to the largest value in the panel;
a solid line separates sections of groups and a dashed line the groups inside a
section; the group names stand in staggered columns left of the gene names.
"""
import textwrap

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from alphapepttools.pl.figure import create_figure
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation, blended_transform_factory

import figure_style as fs


def ttest_dotplot(groups, tables, gene_map, comparisons, labels, letter, name, folder,
                  size=(1, 3), separator_before=None, cbar_label="log2 FC",
                  dot_min_pt=1.0, dot_max_pt=4.0, group_gap_pt=0.0, stagger_gap_pt=0.0,
                  key_pad=1):
    """Draw the dot plot and save it with figure_style.save_panel.

    groups           list of (group name, [gene symbols]), in row order
    tables           dict: comparison -> t-test table with the columns protein,
                     log2fc and fdr
    gene_map         dict: protein ID -> gene symbol. A gene is matched to the
                     first table row whose protein maps to it.
    comparisons      keys of `tables`, in column order
    labels           column names, one per comparison
    letter, name, folder   passed to save_panel
    size             panel size in units (columns, rows)
    separator_before k: a dashed line before column k
    dot_min_pt, dot_max_pt  printed dot diameters at FDR = 1 and at the smallest
                     FDR in the panel; the area is linear in -log10 FDR
    group_gap_pt     printed gap from the gene names to the group labels
    stagger_gap_pt   printed gap between the two columns of group labels
    key_pad          distance of the key from the lower left corner of the axes,
                     in font sizes, to the left and down
    """
    dot_min, dot_max = fs.native_pt(dot_min_pt), fs.native_pt(dot_max_pt)
    group_gap, stagger_gap = fs.native_pt(group_gap_pt), fs.native_pt(stagger_gap_pt)
    tables = {comp: tables[comp].assign(gene=tables[comp]["protein"].map(gene_map))
              for comp in comparisons}

    rows, protein_order, missing = [], [], []
    for group_name, genes in groups:
        for gene in genes:
            protein_order.append(gene)
            found = False
            for comp, comp_label in zip(comparisons, labels):
                match = tables[comp][tables[comp]["gene"] == gene]
                found = found or not match.empty
                rows.append({"protein": gene, "comparison": comp_label,
                             "log2fc": match.iloc[0]["log2fc"] if not match.empty else np.nan,
                             "fdr":    match.iloc[0]["fdr"]    if not match.empty else np.nan})
            if not found:
                missing.append((group_name, gene))
    if missing:
        print(f"  {len(missing)} gene(s) not quantified in any comparison, drawn as a grey cross:")
        for grp, gene in missing:
            print(f"      {gene:<12s} ({grp})")

    dot_df = pd.DataFrame(rows)
    dot_df["-log10fdr"] = -np.log10(dot_df["fdr"].clip(lower=1e-10))
    vmax = float(np.nanmax(np.abs(dot_df["log2fc"])))
    size_max = float(dot_df["-log10fdr"].max())
    if not (np.isfinite(vmax) and vmax > 0 and np.isfinite(size_max) and size_max > 0):
        raise ValueError("no finite log2fc or FDR in the panel; check that the gene "
                         "symbols match the values of gene_map")

    def area(v):
        """Marker area for -log10 FDR = v."""
        return dot_min ** 2 + (v / size_max) * (dot_max ** 2 - dot_min ** 2)

    n_genes, n_comps = len(protein_order), len(labels)
    fig, axm = create_figure(1, 1, figsize=fs.figsize(*size))
    ax = axm.next()

    for i, prot in enumerate(protein_order):
        for j, comp_label in enumerate(labels):
            r = dot_df[(dot_df["protein"] == prot) & (dot_df["comparison"] == comp_label)].iloc[0]
            # A cross means there is no test result for that gene in that comparison:
            # the symbol matched no protein group, or the t-test gave no p-value (fewer
            # than two values in a group). It does not mean "not significant".
            if np.isnan(r["fdr"]):
                ax.plot(j, i, marker="x", color=fs.BG_COLOR, markersize=fs.native_pt(2.0),
                        markeredgewidth=fs.MARKER_LW)
                continue
            sig = r["fdr"] < 0.05
            ax.scatter(j, i, s=area(r["-log10fdr"]),
                       c=[fs.DIRECTION_CMAP((r["log2fc"] / vmax + 1) / 2)],
                       edgecolors="black" if sig else "0.6",
                       linewidths=fs.AXES_LW if sig else fs.THIN_LW, zorder=3)

    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(protein_order, rotation=0)
    ax.set_xticks(range(n_comps))
    ax.set_xticklabels(labels, rotation=90)
    ax.set_xlim(-0.5, n_comps - 0.5)
    ax.set_ylim(n_genes - 0.5, -0.5)                 # first gene at the top
    ax.grid(alpha=0.15, zorder=0)

    # Dashed lines between the gene groups, and before column separator_before
    counts = [len(genes) for _, genes in groups]
    for b in np.cumsum(counts)[:-1]:
        ax.axhline(b - 0.5, color=fs.REF_COLOR, linewidth=fs.THIN_LW, linestyle="--", zorder=1)
    if separator_before is not None:
        ax.axvline(separator_before - 0.5, color=fs.REF_COLOR, linewidth=fs.THIN_LW,
                   linestyle="--", zorder=1)

    sm = plt.cm.ScalarMappable(cmap=fs.DIRECTION_CMAP, norm=plt.Normalize(-vmax, vmax))
    cbar = fig.colorbar(sm, ax=ax, shrink=0.25, pad=0.02)
    cbar.set_label(cbar_label)

    # A grey outline is not significant; the two black-outlined dots are size
    # examples for significant FDRs.
    legend_handles = [
        Line2D([], [], marker="o", ls="none", markersize=np.sqrt(area(1.3)),
               markerfacecolor="0.7", markeredgecolor="0.6", markeredgewidth=fs.THIN_LW,
               label="not significant"),
    ] + [
        Line2D([], [], marker="o", ls="none", markersize=np.sqrt(area(v)),
               markerfacecolor="0.7", markeredgecolor="black", markeredgewidth=fs.AXES_LW,
               label=label)
        for v, label in [(1.3, "FDR = 0.05"), (4, "FDR = 1e-4")]
    ]
    if dot_df["fdr"].isna().any():
        legend_handles.append(
            Line2D([], [], marker="x", ls="none", markersize=fs.native_pt(2.0), color=fs.BG_COLOR,
                   markeredgewidth=fs.MARKER_LW, label="not tested"))
    # The key's upper right corner is at the lower left corner of the axes, below
    # the gene names and left of the column names.
    ax.legend(handles=legend_handles, loc="upper right", bbox_to_anchor=(0, 0),
              borderaxespad=key_pad, frameon=False)

    # Group labels left of the gene names, rotated, staggered in two columns, so a
    # label may run past its own block of genes: each is wrapped to fit its block
    # plus half of each neighbour. Measured after the colourbar and the key are in
    # place, because those set the size of the axes.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px_to_pt = 72.0 / fig.dpi
    ab = ax.get_window_extent()
    gene_pitch_pt = ab.height / n_genes * px_to_pt

    def text_size_pt(txt):
        """(width, height) in points of txt set as a bold group label."""
        t = fig.text(0, 0, txt, fontsize=fs.TICK_FS, fontweight="bold")
        b = t.get_window_extent(renderer)
        t.remove()
        return b.width * px_to_pt, b.height * px_to_pt

    def fit_wrap(label, avail_pt):
        if "\n" in label:          # a line break written in the group name is kept
            return label
        for w in range(len(label), 3, -1):
            wrapped = textwrap.fill(label, width=w, break_long_words=False)
            if max(text_size_pt(line)[0] for line in wrapped.split("\n")) <= avail_pt:
                return wrapped
        return textwrap.fill(label, width=4, break_long_words=False)

    wrapped = []
    for k, (group_name, _) in enumerate(groups):
        span = (counts[k] + (counts[k - 1] / 2 if k > 0 else 0)
                + (counts[k + 1] / 2 if k + 1 < len(counts) else 0))
        wrapped.append(fit_wrap(group_name, span * gene_pitch_pt * 0.92))

    # Distances in points from the left axes edge. A rotated label extends to the
    # left by the height of its unrotated text. Groups 1, 3, 5, ... are in the
    # inner column, groups 2, 4, 6, ... in the outer one; each column is as deep
    # as its deepest label.
    genes_pt = (ab.x0 - min(t.get_window_extent(renderer).x0
                            for t in ax.get_yticklabels())) * px_to_pt
    depths = [text_size_pt(w)[1] for w in wrapped]
    inner_d = max(depths[0::2])
    outer_d = max(depths[1::2]) if len(depths) > 1 else 0.0
    inner_centre = genes_pt + group_gap + inner_d / 2
    outer_centre = genes_pt + group_gap + inner_d + stagger_gap + outer_d / 2

    mids = np.cumsum([0] + counts[:-1]) + (np.array(counts) - 1) / 2
    for k, (mid, txt) in enumerate(zip(mids, wrapped)):
        offset = -(outer_centre if k % 2 else inner_centre)
        ax.annotate(txt, xy=(0, mid),
                    xycoords=blended_transform_factory(ax.transAxes, ax.transData),
                    xytext=(offset, 0), textcoords="offset points",
                    rotation=90, va="center", ha="center",
                    fontsize=fs.TICK_FS, fontweight="bold", annotation_clip=False)

    fs.save_panel(fig, letter, name, folder)


def tested_dotplot(groups, table, contrasts, labels, letter, name, folder, size=(1, 2),
                   fdr_col="fdr", vmax=1.5, size_cap=4, dot_min_pt=1.0, dot_max_pt=4.0,
                   gene_fs=None, bracket_gap_pt=1.0, label_gap_pt=1.0, bracket_trim=0.25,
                   key_spacing=0.3):
    """Draw the dot plot and save it with figure_style.save_panel.

    groups           list of (group name, [gene symbols]), in row order; "\\n" in a
                     group name starts a second line
    table            long table with the columns gene, contrast, tested, log2fc and
                     fdr_col; one row per gene and contrast
    contrasts        values of table["contrast"], in column order
    labels           column names, one per contrast
    letter, name, folder   passed to save_panel
    size             panel size in units (columns, rows)
    fdr_col          the FDR column behind dot area and outline
    vmax             the colour scale ends at +-vmax; larger values take the end colour
    size_cap         dot area is linear in -log10 FDR from FDR = 1 up to
                     FDR = 10^-size_cap, and constant beyond
    dot_min_pt, dot_max_pt  printed dot diameters at FDR = 1 and at the cap
    gene_fs          font size of the gene names; default figure_style.TICK_FS
    bracket_gap_pt   printed gap from the gene names to the group brackets
    label_gap_pt     printed gap from the brackets to the group names
    bracket_trim     each bracket ends this many rows inside its group's boundary
    key_spacing      vertical space between the key entries, in font sizes
    """
    dot_min, dot_max = fs.native_pt(dot_min_pt), fs.native_pt(dot_max_pt)
    gene_fs = fs.TICK_FS if gene_fs is None else gene_fs
    bracket_gap, label_gap = fs.native_pt(bracket_gap_pt), fs.native_pt(label_gap_pt)
    lookup = table.drop_duplicates(["gene", "contrast"]).set_index(["gene", "contrast"])
    norm = mcolors.Normalize(-vmax, vmax, clip=True)

    def area(v):
        """Marker area for -log10 FDR = v."""
        return dot_min ** 2 + min(v, size_cap) / size_cap * (dot_max ** 2 - dot_min ** 2)

    gene_order = [g for _, genes in groups for g in genes]
    absent = [g for g in gene_order if g not in lookup.index.get_level_values("gene")]
    if absent:
        print(f"  not in the table, drawn as grey crosses: {absent}")

    fig, axm = create_figure(1, 1, figsize=fs.figsize(*size))
    ax = axm.next()
    for i, gene in enumerate(gene_order):
        for j, contrast in enumerate(contrasts):
            r = lookup.loc[(gene, contrast)] if (gene, contrast) in lookup.index else None
            # A cross means there is no test result: the protein was not tested in
            # this contrast. It does not mean "not significant".
            if r is None or not r["tested"] or not np.isfinite(r[fdr_col]):
                ax.plot(j, i, marker="x", color=fs.BG_COLOR, markersize=fs.native_pt(1.8),
                        markeredgewidth=fs.MARKER_LW)
                continue
            sig = r[fdr_col] < 0.05
            ax.scatter(j, i, s=area(-np.log10(r[fdr_col])),
                       c=[fs.DIRECTION_CMAP(norm(r["log2fc"]))],
                       edgecolors="black" if sig else "0.6",
                       linewidths=fs.AXES_LW if sig else fs.THIN_LW, zorder=3)

    n_genes, n_comps = len(gene_order), len(contrasts)
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(gene_order, fontsize=gene_fs)
    ax.set_xticks(range(n_comps))
    ax.set_xticklabels(labels, rotation=90)
    ax.set_xlim(-0.5, n_comps - 0.5)
    ax.set_ylim(n_genes - 0.5, -0.5)                 # first gene at the top
    ax.grid(alpha=0.15, zorder=0)

    counts = [len(genes) for _, genes in groups]
    for b in np.cumsum(counts)[:-1]:
        ax.axhline(b - 0.5, color=fs.REF_COLOR, linewidth=fs.THIN_LW, linestyle="--", zorder=1)

    sm = plt.cm.ScalarMappable(cmap=fs.DIRECTION_CMAP, norm=norm)
    cbar = fig.colorbar(sm, ax=ax, shrink=0.3, pad=0.02, extend="both")
    # Colourbar label under the bar, each line centred. It is a text in the
    # colourbar axes and not an axis label: constrained layout counts the width of
    # axes text, but not of an axis label. It starts below the arrow tip, which
    # reaches 0.05 of the bar length past its end.
    cbar.ax.text(0.0, -0.08, "log2\nFC", transform=cbar.ax.transAxes,
                 ha="left", va="top", multialignment="center", fontsize=fs.LABEL_FS)

    legend_handles = [
        Line2D([], [], marker="o", ls="none", markersize=np.sqrt(area(1.3)),
               markerfacecolor="0.7", markeredgecolor="0.6", markeredgewidth=fs.THIN_LW,
               label="not\nsignificant"),          # two lines: the key ends left of the contrast names
    ] + [
        Line2D([], [], marker="o", ls="none", markersize=np.sqrt(area(v)),
               markerfacecolor="0.7", markeredgecolor="black", markeredgewidth=fs.AXES_LW,
               label=label)
        for v, label in [(1.3, "FDR = 0.05"), (size_cap, f"FDR ≤ 1e-{size_cap}")]
    ] + [
        Line2D([], [], marker="x", ls="none", markersize=fs.native_pt(1.8), color=fs.BG_COLOR,
               markeredgewidth=fs.MARKER_LW, label="not tested"),
    ]
    # The key is placed in the lower left corner of the canvas, left of the
    # contrast names. It is a figure legend and not part of the layout, so it takes
    # no space from the axes; the checks at the end report any overlap with it.
    key = fig.legend(handles=legend_handles, loc="lower left", bbox_to_anchor=(0, 0),
                     borderaxespad=0, labelspacing=key_spacing, frameon=False,
                     handlelength=1.0, handletextpad=0.5)

    # Group labels left of the gene names, rotated, in one column. A bracket spans
    # the rows of each group, so a label longer than its group still shows which
    # rows it names. Measured after the colourbar and the key are in place,
    # because those set the size of the axes.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    genes_pt = (ax.get_window_extent().x0
                - min(t.get_window_extent(renderer).x0 for t in ax.get_yticklabels())) * 72 / fig.dpi
    bracket_pt = genes_pt + bracket_gap
    rows_y = blended_transform_factory(ax.transAxes, ax.transData)
    group_texts = []
    for (group_name, genes), start in zip(groups, np.cumsum([0] + counts[:-1])):
        end = start + len(genes) - 1
        ax.plot([0, 0], [start - 0.5 + bracket_trim, end + 0.5 - bracket_trim],
                color="black", linewidth=fs.THIN_LW,
                clip_on=False,
                transform=rows_y + ScaledTranslation(-bracket_pt / 72, 0, fig.dpi_scale_trans))
        group_texts.append(ax.annotate(
            group_name, xy=(0, (start + end) / 2), xycoords=rows_y,
            xytext=(-(bracket_pt + label_gap), 0), textcoords="offset points",
            rotation=90, ha="right", va="center", multialignment="center",
            fontsize=fs.TICK_FS, fontweight="bold", annotation_clip=False))

    # Capitals and digits are 0.72 font sizes tall; the CHECK warns when rows are closer than 0.75 font sizes.
    fig.canvas.draw()
    _pitch = ax.get_window_extent().height / fig.dpi * 72 / n_genes
    print(f"CHECK row pitch {_pitch / fs.DRAW_SCALE:.1f} pt for {gene_fs / fs.DRAW_SCALE:.1f} pt "
          f"gene names, {n_genes} rows (printed)"
          + ("   !! the names overlap" if _pitch < 0.75 * gene_fs else ""))
    renderer = fig.canvas.get_renderer()
    _boxes = [t.get_window_extent(renderer) for t in group_texts]
    _clash = [f"{groups[k][0]!r} / {groups[k + 1][0]!r}"
              for k in range(len(_boxes) - 1) if _boxes[k].overlaps(_boxes[k + 1])]
    print("CHECK neighbouring group labels overlap: "
          + (f"!! {_clash}; shorten or split them with \\n" if _clash else "none"))
    _key_box = key.get_window_extent(renderer)
    _gap = min(b.y0 for b in _boxes) - _key_box.y1
    _xgap = min(t.get_window_extent(renderer).x0 for t in ax.get_xticklabels()) - _key_box.x1
    print(f"CHECK key clear of the group labels: {_gap / fig.dpi * 25.4 / fs.DRAW_SCALE:.1f} mm, "
          f"of the contrast names: {_xgap / fig.dpi * 25.4 / fs.DRAW_SCALE:.1f} mm"
          + ("   !! overlap" if min(_gap, _xgap) < 0 else ""))

    fs.save_panel(fig, letter, name, folder)


def cross_system_dotplot(df_lfc, df_fdr, columns, groups, letter, name, folder,
                         size=(1, 3), section_breaks=None, group_layers=2,
                         legend_corner=False, dot_min_pt=1.0, dot_max_pt=4.0,
                         cross_pt=2.0, group_gap_pt=1.5, wrap=18):
    """Draw the dot plot and save it with figure_style.save_panel.

    df_lfc, df_fdr   DataFrames, gene symbols (index) x columns; a NaN FDR draws
                     the grey cross
    columns          column names of df_lfc and df_fdr, in drawing order
    groups           list of (group name, [gene symbols]), in row order
    letter, name, folder   passed to save_panel
    size             panel size in units (columns, rows)
    section_breaks   row indices where a new section starts (solid line)
    group_layers     number of staggered columns of group names; group k stands in
                     column k % group_layers
    legend_corner    the dot-size key in the corner below the gene names and left
                     of the column names; otherwise below the panel
    dot_min_pt, dot_max_pt  printed dot diameters at FDR = 1 and at the smallest
                     FDR in the panel; the area is linear in -log10 FDR
    cross_pt         printed size of the grey cross
    group_gap_pt     printed gap from the gene names to the first column of group
                     names, and between the columns
    wrap             characters per line of a group name
    """
    dot_min, dot_max = fs.native_pt(dot_min_pt), fs.native_pt(dot_max_pt)
    group_gap = fs.native_pt(group_gap_pt)
    section_breaks = section_breaks or []
    gene_order = [g for _, genes in groups for g in genes]

    rows = []
    for gene in gene_order:
        for col in columns:
            rows.append({
                "gene": gene, "column": col,
                "log2fc": df_lfc.loc[gene, col] if gene in df_lfc.index else np.nan,
                "fdr":    df_fdr.loc[gene, col] if gene in df_fdr.index else np.nan,
            })
    dot_df = pd.DataFrame(rows)
    dot_df["-log10fdr"] = -np.log10(dot_df["fdr"].clip(lower=1e-10))

    n_genes, n_cols = len(gene_order), len(columns)
    if n_genes == 0:
        print(f"{name}: nothing to plot.")
        return

    vmax = np.nanmax(np.abs(dot_df["log2fc"]))
    size_max = dot_df["-log10fdr"].replace([np.inf], np.nan).max() or 1.0

    def area(v):
        return dot_min ** 2 + (v / size_max) * (dot_max ** 2 - dot_min ** 2)

    fig, axm = create_figure(1, 1, figsize=fs.figsize(*size))
    ax = axm.next()

    cell = {(r["gene"], r["column"]): r for _, r in dot_df.iterrows()}
    for i, gene in enumerate(gene_order):
        for j, col in enumerate(columns):
            row = cell[(gene, col)]
            if not np.isfinite(row["fdr"]):
                ax.plot(j, i, marker="x", color=fs.BG_COLOR, markersize=fs.native_pt(cross_pt),
                        markeredgewidth=fs.MARKER_LW)
                continue
            sig = row["fdr"] < 0.05
            color = fs.DIRECTION_CMAP((row["log2fc"] / vmax + 1) / 2)
            ax.scatter(j, i, s=area(row["-log10fdr"]), c=[color],
                       edgecolors="black" if sig else "0.6",
                       linewidths=fs.AXES_LW if sig else fs.THIN_LW, zorder=3)

    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(gene_order)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(columns, rotation=90)
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(-0.5, n_genes - 0.5)
    ax.invert_yaxis()
    ax.grid(axis="both", alpha=0.15, zorder=0)
    ax.set_axisbelow(True)

    # dashed line between groups, solid line between sections
    running = 0
    group_mids, group_names = [], []
    for k, (label, members) in enumerate(groups):
        if k > 0 and running not in section_breaks:
            ax.axhline(running - 0.5, color=fs.REF_COLOR, linewidth=fs.THIN_LW,
                       linestyle="--", zorder=1)
        group_mids.append(running + (len(members) - 1) / 2)
        group_names.append(label)
        running += len(members)
    for r in section_breaks:
        ax.axhline(r - 0.5, color="black", linewidth=fs.AXES_LW, zorder=2)

    layer_axes = []
    for layer in range(group_layers):
        ax2 = ax.secondary_yaxis("left")
        ax2.set_yticks([m for k, m in enumerate(group_mids) if k % group_layers == layer])
        ax2.set_yticklabels([textwrap.fill(g, wrap)
                             for k, g in enumerate(group_names) if k % group_layers == layer],
                            fontweight="bold", rotation=90, va="center", ha="center")
        ax2.tick_params(length=0)
        for spine in ax2.spines.values():
            spine.set_visible(False)
        layer_axes.append(ax2)

    ax.set_ylabel("")
    ax.set_xlabel("")

    sm = plt.cm.ScalarMappable(cmap=fs.DIRECTION_CMAP, norm=plt.Normalize(-vmax, vmax))
    sm.set_array([])
    # fraction: the share of the dot axes' width that the colourbar takes
    cbar = plt.colorbar(sm, ax=ax, fraction=0.09, pad=0.05)
    cbar.set_label("log2 fold change")

    # dot-size key, from the same area formula as the dots
    for val, label in [(1.3, "FDR<0.05"), (4, "FDR<1e-4")]:
        ax.scatter([], [], s=area(val), c="0.7", edgecolors="black",
                   linewidths=fs.AXES_LW, label=label)
    if legend_corner:
        ax.legend(loc="upper right", bbox_to_anchor=(0, 0), title="Dot size",
                  frameon=False, handlelength=0.7, handletextpad=0.3)
    else:
        fig.legend(loc="outside lower center", title="Dot size", frameon=False,
                   ncol=2, borderaxespad=0)

    # The group-name columns start where the drawn gene names end (tick length +
    # pad + text); each column is as wide as its widest rotated name, and a name
    # is anchored at its centre.
    to_pt = 72 / fig.dpi
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    w_gene = (max(t.get_window_extent(rend).width for t in ax.get_yticklabels()) * to_pt
              + ax.yaxis.get_major_ticks()[0].get_tick_padding()
              + plt.rcParams["ytick.major.pad"])
    pad = w_gene
    for ax2 in layer_axes:
        w = max((t.get_window_extent(rend).width for t in ax2.get_yticklabels()
                 if t.get_text()), default=0) * to_pt
        ax2.tick_params(pad=pad + group_gap + w / 2)
        pad += group_gap + w
    # constrained layout places the axes once per draw: one more draw for the
    # final label columns before save_panel measures the panel
    fig.canvas.draw()

    fs.save_panel(fig, letter, name, folder)
