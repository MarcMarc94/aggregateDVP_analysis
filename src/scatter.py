"""Replication scatter: one effect per item in the discovery cohort (x) against the
replication cohort (y), items significant in the discovery cohort drawn prominent.

replication_scatter() draws the zero lines, the identity line and the points on equal
axes, and one text block with Pearson r, n and sign agreement: for the significant
items first, for all items in brackets. Sizes follow figure_style: printed points and
panel units.
"""
import numpy as np
from alphapepttools.pl.figure import create_figure, label_axes
from matplotlib.lines import Line2D
from scipy import stats

import figure_style as fs

PAD     = 1.08      # axis limit as a multiple of the largest |value|
LEG_DOT = 3.5       # printed diameter of the legend markers, larger than the data dots


def replication_scatter(x, y, sig, letter, name, folder, xlabel, ylabel, hit_color,
                        alpha=0.05, ns_open=True, ns_pt=1.6, hit_pt=2.2,
                        all_label="sets tested"):
    """Draw the panel and save it with figure_style.save_panel.

    x, y        effects in the discovery and the replication cohort, one per item
    sig         boolean, significant in the discovery cohort
    letter, name, folder   passed to save_panel
    hit_color   colour of the significant items
    alpha       significance threshold, named in the legend
    ns_open     True: items not significant as open circles; False: as small grey dots
    ns_pt, hit_pt   printed diameters of the two kinds of dots
    all_label   what the items are, for the text block when fewer than three are
                significant ("sets tested", "shared proteins")
    """
    x, y, sig = np.asarray(x, float), np.asarray(y, float), np.asarray(sig, bool)
    same = np.sign(x) == np.sign(y)
    lim = float(np.nanmax(np.abs(np.concatenate([x, y])))) * PAD
    ns_dot, hit_dot, leg_dot = fs.native_pt(ns_pt), fs.native_pt(hit_pt), fs.native_pt(LEG_DOT)

    fig, axm = create_figure(1, 1, figsize=fs.figsize(1, 1))
    ax = axm.next()

    # Zero lines and the identity line, behind the points.
    ax.axhline(0, color=fs.REF_COLOR, lw=fs.THIN_LW, zorder=1)
    ax.axvline(0, color=fs.REF_COLOR, lw=fs.THIN_LW, zorder=1)
    ax.plot([-lim, lim], [-lim, lim], ls="--", color=fs.REF_COLOR, lw=fs.AXES_LW, zorder=2)

    if ns_open:
        ax.scatter(x[~sig], y[~sig], s=ns_dot ** 2, facecolors="none",
                   edgecolors=fs.BG_COLOR, linewidths=fs.LINE_LW, zorder=3)
        ns_handle = dict(markerfacecolor="none", markeredgecolor=fs.BG_COLOR,
                         markeredgewidth=fs.LINE_LW)
    else:
        ax.scatter(x[~sig], y[~sig], s=ns_dot ** 2, color=fs.BG_COLOR, alpha=0.30,
                   linewidths=0, zorder=3)
        ns_handle = dict(markerfacecolor=fs.BG_COLOR, markeredgecolor="none")
    ax.scatter(x[sig], y[sig], s=hit_dot ** 2, color=hit_color,
               edgecolors="white", linewidths=fs.MARKER_LW, zorder=4)

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    label_axes(ax, xlabel=xlabel, ylabel=ylabel)

    # One left-aligned block, no column stops. The header names the two figures once,
    # significant first and all items in brackets. Below three points a correlation
    # is not defined (with two, r is always +/-1), so the block then shows all items only.
    r_all, a_all = float(stats.pearsonr(x, y)[0]), float(same.mean())
    if sig.sum() > 2:
        r_sub, a_sub = float(stats.pearsonr(x[sig], y[sig])[0]), float(same[sig].mean())
        txt = (f"significant (all)\n"
               f"r: {r_sub:+.2f} ({r_all:+.2f})\n"
               f"n: {int(sig.sum())} ({len(x)})\n"
               f"sign agreement: \n{100 * a_sub:.0f}% ({100 * a_all:.0f}%)")
    else:
        txt = (f"all {len(x)} {all_label}\n"
               f"r: {r_all:+.2f}\n"
               f"sign agreement: {100 * a_all:.0f}%")
        print(f"NOTE: {int(sig.sum())} significant items, too few for a subset "
              f"correlation; the block shows all items only.")
    ann = ax.text(0.03, 0.97, txt, transform=ax.transAxes, ha="left", va="top",
                  fontsize=fs.ANNOT_FS, color="black", linespacing=1.2, zorder=6)

    ax.legend(handles=[
        Line2D([0], [0], marker="o", linestyle="none", markersize=leg_dot,
               markerfacecolor=hit_color, markeredgecolor="white",
               markeredgewidth=fs.MARKER_LW, label=f"BH < {alpha}"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=leg_dot,
               label="n.s.", **ns_handle)],
        loc="lower right", frameon=False, handletextpad=0.5, labelspacing=0.3)

    # The block has no columns that could collide, so the check is whether it covers
    # data points. Measured on the rendered text.
    fig.canvas.draw()
    bb = ann.get_window_extent().transformed(ax.transAxes.inverted())
    fx, fy = (x + lim) / (2 * lim), (y + lim) / (2 * lim)
    hits = int(((fx >= bb.x0) & (fx <= bb.x1) & (fy >= bb.y0) & (fy <= bb.y1)).sum())
    print(f"CHECK annotation box spans x {bb.x0:.2f}-{bb.x1:.2f}, y {bb.y0:.2f}-{bb.y1:.2f} "
          f"in axes fractions; plotted points inside it: {hits}")
    if hits:
        print("  Points sit under the text. Move the block to the lower left "
              "(va='bottom', y=0.03), or split 'sign agreement' onto its own line.")
    if bb.x1 > 1.0 or bb.y0 < 0.0:
        print("  The block runs outside the axes; shorten a line.")
    mx = float(np.nanmax(np.abs(np.concatenate([x, y]))))
    print(f"CHECK axis limits {-lim:.2f} to {lim:.2f}; largest plotted {mx:.2f} "
          f"-> inside: {mx <= lim}")
    print(f"CHECK equal aspect: {ax.get_aspect()}")

    fs.save_panel(fig, letter, name, folder)
