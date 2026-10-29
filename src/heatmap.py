"""Heatmap of one value per gene and comparison, in named gene groups.

group_heatmap() draws one row per gene and one column per comparison, coloured
by the value (log2 fold change, or row z-score in BrainTissue 4c) on
figure_style.DIRECTION_CMAP. Missing cells are grey. A
line separates neighbouring groups, and the group names are drawn rotated left
of the gene names. Sizes follow figure_style: printed points and panel units.
"""
import matplotlib.colors as mcolors
import numpy as np
import seaborn as sns
from alphapepttools.pl.figure import create_figure
from matplotlib.transforms import ScaledTranslation

import figure_style as fs


def group_heatmap(lfc_matrix, groups, letter, name, folder, size=(1, 1), norm=None,
                  gene_fs=None, label_gap_pt=1.5, stagger=1, align_ends=False,
                  block_after=None, label_above=None, cbar_shrink=0.4, col_linespacing=None,
                  cbar_label="log2\nFC", col_rotation=0):
    """Draw the heatmap and save it with figure_style.save_panel.

    lfc_matrix       DataFrame of values, gene symbols (index) x comparisons (columns, whose
                     names are drawn); genes that are in no group are left out
    groups           list of (group name, genes), in row order; a gene absent from
                     lfc_matrix is left out and reported. An empty name draws no label.
    letter, name, folder   passed to save_panel
    size             panel size in units (columns, rows)
    norm             colour normalisation; default symmetric around 0 at the
                     largest |value| drawn
    gene_fs          font size of the gene names; default figure_style.TICK_FS
    label_gap_pt     printed gap from the gene names to the group names
    stagger          number of depths the group names alternate between, so that a
                     name longer than its block does not run into its neighbour
    align_ends       align the first group name to the top edge of its block and
                     the last to the bottom edge, so neither runs past the panel
    block_after      number of groups after which a thicker line is drawn
    label_above      write the value into every cell whose |value|
                     exceeds this, for a colour scale capped below the data
    cbar_shrink      colourbar length as a fraction of the heatmap height
    col_linespacing  line spacing of multi-line column names; default matplotlib's
    cbar_label       text under the colourbar
    col_rotation     rotation of the column names, in degrees
    """
    gene_fs = fs.TICK_FS if gene_fs is None else gene_fs
    present = [[g for g in genes if g in lfc_matrix.index] for _, genes in groups]
    lost = [g for _, genes in groups for g in genes if g not in lfc_matrix.index]
    m = lfc_matrix.loc[[g for genes in present for g in genes]]
    if norm is None:
        vmax = float(np.nanmax(np.abs(m.values)))
        norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap = fs.DIRECTION_CMAP.with_extremes(bad=fs.BG_COLOR)     # missing, not zero

    fig, axm = create_figure(1, 1, figsize=fs.figsize(*size))
    ax = axm.next()
    ax.grid(False)
    sns.heatmap(m, ax=ax, cmap=cmap, norm=norm, linewidths=0,
                cbar_kws={"shrink": cbar_shrink, "pad": 0.02}, annot=False)

    # Colourbar label under the bar, each line centred. The block starts at the
    # bar's left edge, so it does not reach over the heatmap. It is a text in the
    # colourbar axes and not an x-axis label: constrained layout counts the width
    # of axes text, but not of an axis label.
    cbar = ax.collections[0].colorbar
    cbar.ax.text(0.0, -0.02, cbar_label, transform=cbar.ax.transAxes,
                 ha="left", va="top", multialignment="center", fontsize=fs.LABEL_FS)

    # Column names at col_rotation, each line centred. Set through tick_params, not
    # set_xticklabels: the latter makes the locator fixed and emits a warning
    # whenever the tick count is touched afterwards.
    ax.tick_params(axis="x", rotation=col_rotation)
    ax.tick_params(axis="both", length=0)       # every row and column is labelled; tick marks would take 0.9 mm
    for t in ax.get_xticklabels():
        t.set_multialignment("center")
        if col_linespacing is not None:
            t.set_linespacing(col_linespacing)

    # Every row labelled: seaborn thins the tick labels of a tall heatmap.
    ax.set_yticks(np.arange(len(m.index)) + 0.5)
    ax.set_yticklabels(m.index, rotation=0, fontsize=gene_fs)

    if label_above is not None:
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                v = m.iat[i, j]
                if np.isfinite(v) and abs(v) > label_above:
                    # The baseline is 0.36 font sizes below the row centre, and digits
                    # are 0.72 font sizes tall, so they are centred. va="center"
                    # centres the line box, whose space for descenders puts digits
                    # too high.
                    ax.annotate(f"{v:.1f}".replace("-", "−"), (j + 0.5, i + 0.5),
                                xytext=(0, -0.36 * gene_fs), textcoords="offset points",
                                ha="center", va="baseline",
                                fontsize=gene_fs, color="white", fontweight="bold")

    # Separators between groups, over the genes present: boundaries taken from the
    # full gene lists put them on the wrong rows when a listed gene is missing.
    sizes = [len(genes) for genes in present]
    bounds = np.cumsum(sizes)
    block = bounds[block_after - 1] if block_after else None
    for b in sorted(set(bounds[:-1])):
        ax.axhline(b, color="black", linewidth=fs.AXES_LW if b == block else fs.THIN_LW)

    # Group names left of the gene names. The gap is measured from the rendered
    # gene names: from the axis line to the left edge of the widest name.
    mids = [end - n / 2 for n, end in zip(sizes, bounds) if n > 0]
    names = [g for (g, _), n in zip(groups, sizes) if n > 0]
    fig.canvas.draw()
    names_pt = (ax.get_window_extent().x0
                - min(t.get_window_extent().x0 for t in ax.get_yticklabels())) / fig.dpi * 72
    ticks = list(mids)
    if align_ends:
        ticks[0], ticks[-1] = 0, len(m)          # the y axis of a heatmap runs downward
    ax2 = ax.secondary_yaxis("left")
    ax2.set_yticks(ticks)
    ax2.set_yticklabels(names, fontweight="bold", rotation=90, va="center", ha="right", ma="center")
    if align_ends:
        labels = ax2.get_yticklabels()
        labels[0].set_verticalalignment("top")
        labels[-1].set_verticalalignment("bottom")
    ax2.tick_params(length=0, pad=names_pt + fs.native_pt(label_gap_pt))
    # Axis lines are drawn with projecting caps, which carry the line half its
    # width (0.3 pt) past the top and bottom of the heatmap; butt caps end it there.
    ax2.spines["left"].set_capstyle("butt")

    # A name rotated 90 degrees runs along its block, so a long name is taller than
    # a short block. Neighbouring names are set at different depths, each depth the
    # width of the widest name plus a gap, measured.
    if stagger > 1:
        fig.canvas.draw()
        step = (max(t.get_window_extent().width for t in ax2.get_yticklabels() if t.get_text())
                / fig.dpi * 72 + fs.native_pt(1.0))
        for i, lbl in enumerate(ax2.get_yticklabels()):
            if i % stagger:
                lbl.set_transform(lbl.get_transform()
                                  + ScaledTranslation(-(i % stagger) * step / 72, 0, fig.dpi_scale_trans))

    ax.set_ylabel("")
    ax.set_xlabel("")

    # Capitals and digits are 0.72 font sizes tall; the CHECK warns when rows are closer than 0.75 font sizes.
    fig.canvas.draw()
    pitch = ax.get_window_extent().height / fig.dpi * 72 / max(len(m), 1)
    n_nan = int(np.isnan(m.values).sum())
    print(f"CHECK matrix {m.shape[0]} genes x {m.shape[1]} comparisons, {size[0]} x {size[1]} units")
    print(f"CHECK missing cells, drawn grey: {n_nan} ({100 * n_nan / m.size:.1f}%)")
    print(f"CHECK requested genes not shown: {len(lost)}" + (f" -> {lost[:10]}" if lost else ""))
    print(f"CHECK row pitch {pitch / fs.DRAW_SCALE:.1f} pt for {gene_fs / fs.DRAW_SCALE:.1f} pt "
          f"gene names, {len(m)} rows (printed)"
          + ("   !! the names overlap" if pitch < 0.75 * gene_fs else ""))

    fs.save_panel(fig, letter, name, folder)
