"""Box plot of one value per group: Disease units as coloured dots over a box, Control
units as black crosses, with two-level x labels.

fraction_boxplot() draws one box per group from the Disease units only. The Control
units are drawn as crosses on top, without a box. Sizes follow figure_style: printed
points and panel units.
"""
import numpy as np
from alphapepttools.pl.figure import create_figure
from matplotlib.lines import Line2D

import figure_style as fs

BOX_W     = 0.6
BOX_ALPHA = 0.55    # box fill behind the dots
DOT_ALPHA = 0.9
JITTER    = 0.18    # largest horizontal offset of a dot, in x units


def fraction_boxplot(data, group_col, value_col, order, colors, sub_labels, groups,
                     letter, name, folder, unit="donor", dot_pt=2.4, legend_dot_pt=None,
                     legend_handlelength=0.5, ylabel="log2 membrane / nuclear", title=None):
    """Draw the panel and save it with figure_style.save_panel.

    data            DataFrame, one row per unit (donor or sample), with group_col,
                    value_col and disease_binary ("Disease" or "Control")
    order           groups in x order; colors maps each group to its colour
    sub_labels      label under each box (DAB state, EA shape, or the group name)
    groups          list of (name, members): the bold second row of x labels, centred
                    under its members; a group in no entry is named by sub_labels only
    letter, name, folder   passed to save_panel
    unit            "donor" or "sample", named in the legend
    dot_pt          printed diameter of a data dot
    legend_dot_pt   printed diameter of the legend markers; default dot_pt
    legend_handlelength   legend handle length in font sizes; it must hold the marker
    title           panel title, or None
    """
    dot = fs.native_pt(dot_pt)
    leg_dot = fs.native_pt(dot_pt if legend_dot_pt is None else legend_dot_pt)
    status = data["disease_binary"].astype(str)
    disease = [data.loc[(status == "Disease") & (data[group_col] == g), value_col].dropna().values
               for g in order]
    control = [data.loc[(status == "Control") & (data[group_col] == g), value_col].dropna().values
               for g in order]
    x = np.arange(len(order))

    fig, axm = create_figure(1, 1, figsize=fs.figsize(1, 1))
    ax = axm.next()

    # Line widths come from the style: box edges, whiskers and caps MARKER_LW, medians LINE_LW.
    for xp, g, vals in zip(x, order, disease):
        if len(vals):
            bp = ax.boxplot([vals], positions=[xp], widths=BOX_W,
                            patch_artist=True, showfliers=False, zorder=2)
            bp["boxes"][0].set_facecolor(colors[g])
            bp["boxes"][0].set_alpha(BOX_ALPHA)

    # One generator for all dots, Disease first, then Control, so the jitter is fixed.
    rng = np.random.default_rng(42)
    for xp, g, vals in zip(x, order, disease):
        if len(vals):
            ax.scatter(xp + rng.uniform(-JITTER, JITTER, size=len(vals)), vals,
                       color=colors[g], s=dot ** 2, alpha=DOT_ALPHA, edgecolor="black",
                       linewidth=fs.MARKER_LW, zorder=3)
    for xp, vals in zip(x, control):
        if len(vals):
            ax.scatter(xp + rng.uniform(-JITTER, JITTER, size=len(vals)), vals,
                       marker="x", color="black", s=dot ** 2, alpha=DOT_ALPHA,
                       linewidth=fs.LINE_LW, zorder=4)

    ax.set_axisbelow(True)     # the grid goes below the data
    ax.set_xticks(x)
    ax.set_xticklabels([sub_labels[g] for g in order], rotation=90)
    ax.set_xlim(x[0] - 0.5, x[-1] + 0.5)

    # The group row hangs below the lowest rotated label of the grouped boxes: its pad
    # is measured from the axes edge to that label's far edge, plus 1 pt. A longer
    # label of a group in no entry (Neuropil) is left out of the measurement.
    fig.canvas.draw()
    grouped = {g for _, members in groups for g in members}
    sub_pt = (ax.get_window_extent().y0
              - min(t.get_window_extent().y0
                    for t, g in zip(ax.get_xticklabels(), order) if g in grouped)
              ) / fig.dpi * 72
    second = ax.secondary_xaxis("bottom")
    second.set_xticks([np.mean([order.index(g) for g in members]) for _, members in groups],
                      labels=[name for name, _ in groups], fontweight="bold")
    second.tick_params(length=0, pad=sub_pt + fs.native_pt(1.0))
    second.spines["bottom"].set_visible(False)   # the axes' own spine is drawn already
    ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title, fontsize=fs.TITLE_FS)

    # What a dot and a cross are, top left. The colours need no key: every box is named
    # by its x labels.
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=fs.BG_COLOR,
               markeredgecolor="black", markeredgewidth=fs.MARKER_LW, markersize=leg_dot,
               label=f"Disease {unit}"),
        Line2D([0], [0], marker="x", linestyle="none", color="black",
               markersize=leg_dot, markeredgewidth=fs.LINE_LW, label=f"Control {unit}"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", frameon=False,
              handlelength=legend_handlelength, handletextpad=0.3, borderaxespad=0.3)

    print(f"CHECK {unit}s drawn per group (Disease / Control):")
    for g, d, c in zip(order, disease, control):
        print(f"  {g:<32} {len(d):>3} / {len(c):>3}" + ("   <- no box drawn" if len(d) == 0 else ""))
    values = np.concatenate([v for v in disease + control if len(v)])
    print(f"CHECK y range of the plotted data: {values.min():.2f} to {values.max():.2f}; "
          f"axis {ax.get_ylim()[0]:.2f} to {ax.get_ylim()[1]:.2f}")
    fig.canvas.draw()
    lb = ax.get_legend().get_window_extent().transformed(ax.transData.inverted())
    n_under = 0
    for coll in ax.collections:
        off = coll.get_offsets()
        if len(off):
            n_under += int(((off[:, 0] >= lb.x0) & (off[:, 0] <= lb.x1)
                            & (off[:, 1] >= lb.y0) & (off[:, 1] <= lb.y1)).sum())
    print(f"CHECK {unit}s under the legend: {n_under}")

    fs.save_panel(fig, letter, name, folder)
