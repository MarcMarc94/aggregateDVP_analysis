"""Panels of the cross-system notebook whose drawing code a reader does not need in
order to check the analysis. Sizes follow figure_style: printed points and panel
units. Every function saves its panel with figure_style.save_panel.

organelle_bars()   one panel per organelle; per module one bar per contrast, the
                   group difference of the module score, with the member count and
                   a star for a module in a called module pair
organelle_ylim()   the y limits organelle_bars() uses for a set of organelles
split_dumbbell()   pairs of assembly parts, one row per pair and one dumbbell per
                   contrast between the median fold changes of the two parts
property_levels()  rank correlation of each property with fold change, one subplot
                   per level of aggregation, shared rows and x axis
split_property_heatmap()  property difference between the two parts of each pair,
                   colour for the magnitude and a sign glyph per cell
concept_pair_heatmap()  log2 fold change per gene for pairs of columns, one pair
                   per concept, with significance asterisks and row groups
"""
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from alphapepttools.pl.figure import create_figure
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory, offset_copy

import figure_style as fs


def _sentence_case(label):
    """Sentence case, leaving names alone: ESCRT, Sec61, AP-2, 20S core, TRiC."""
    return label if label[:1].isupper() or label[:1].isdigit() else label[0].upper() + label[1:]


def _drawn_bars(bars, modules, organelles, contrasts, min_members):
    """The bar heights drawn for a set of organelles."""
    return [bars.loc[(c, g), "delta_module"]
            for g in modules[modules.organelle.isin(organelles)].go_id
            for c in contrasts
            if bars.loc[(c, g), "n_members"] >= min_members]


def organelle_ylim(bars, modules, organelles, contrasts, pad=0.30, min_members=4,
                   floor=0.5):
    """(low, high) axis limits for a set of organelles.

    Each side follows its own bars: an organelle whose modules are all depleted
    gets a small range above zero and the rest of the panel below it. `floor` is the
    smallest scale, so an organelle with only small differences is drawn on that
    scale. `pad` is the room for the count and the star beyond the last bar, as a
    fraction of the scale, added per side.
    """
    h = [v for v in _drawn_bars(bars, modules, organelles, contrasts, min_members)
         if np.isfinite(v)]
    hi = max([v for v in h if v > 0], default=0.0)
    lo = min([v for v in h if v < 0], default=0.0)
    scale = max(floor, hi, -lo)
    return (lo - pad * scale if lo else -0.18 * scale,
            hi + pad * scale if hi else 0.18 * scale)


def organelle_bars(bars, modules, organelles, contrasts, colors, labels, short, starred,
                   letter, name, folder, ncols=1, ylim=None, pad=0.30, min_members=4,
                   floor=0.5, wide=(), legend=True, l_in=1.80, r_in=0.35, t_in=0.75,
                   b_in=1.15, zero_eps=0.04, inside=0.45, star_ms=0.76, star_lift=0.58):
    """Draw one panel per organelle, save the figure, and return it.

    bars         DataFrame indexed by (contrast, go_id) with delta_module (bar
                 height) and n_members (the count printed on the bar)
    modules      DataFrame with organelle and go_id, one row per module, in drawing
                 order
    organelles   organelles to draw, one panel each, in reading order
    contrasts    contrast names, one bar each per module, in this order
    colors, labels  per contrast: bar colour, legend text
    short        go_id -> module label on the x axis
    starred      set of (contrast, go_id) that carry a star
    letter, name, folder   passed to save_panel
    ncols        panels per row
    ylim         None (each organelle its own limits from organelle_ylim), one
                 (low, high) pair for all, or a dict of pairs by organelle
    pad, floor   passed to organelle_ylim
    min_members  below this a module has no bar; "n<min" is written instead
    wide         organelles that need a two-unit canvas: module labels stay at
                 TICK_FS and the canvas widens instead
    legend       the colour key at the bottom; False keeps no band free for it
    l_in, r_in, t_in, b_in   margins in inches inside the panel grid: left of the
                 first column (y label and tick labels), right of every column,
                 above every row (title), below every row (module labels)
    zero_eps     a bar lower than this is invisible, so a tick is drawn at zero
    inside       a bar at least this tall, as a fraction of the axis, carries its
                 count inside
    star_ms      star marker size as a fraction of the count font size
    star_lift    star centre above the count, in count font sizes
    """
    fig = plt.figure()
    nrows = int(np.ceil(len(organelles) / ncols))
    slots = max(ncols, 2 if set(organelles) & set(wide) else 1)
    fig.set_size_inches(*fs.figsize(slots, nrows))
    fig_w, fig_h = fig.get_size_inches()

    def _lim(organelle):
        """(low, high) for one organelle: its own limits unless ylim gives them."""
        own = organelle_ylim(bars, modules, [organelle], contrasts, pad, min_members, floor)
        if isinstance(ylim, dict):
            return ylim.get(organelle) or own
        return ylim if ylim is not None else own

    # --- one legend across the whole figure, built before the axes so the bottom
    #     margin can be set from its measured height ----------------------------
    # The contrasts only, in one row; the star rule goes into the figure legend text.
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[c], label=labels[c])
               for c in contrasts]

    def _place_legend(loc, anchor, nc, size):
        return fig.legend(handles=handles, loc=loc, bbox_to_anchor=anchor, ncol=nc,
                          frameon=False, handlelength=0.8, handletextpad=0.3,
                          columnspacing=0.8, fontsize=size)

    # the fewest rows that fit the width, at the style size
    fits = False
    for nc, lfs in [(n, fs.LEGEND_FS) for n in (len(handles), 3, 2, 1)]:
        lgd = _place_legend("upper center", (0.5, 1.0), nc, lfs)
        fig.canvas.draw()
        bb = lgd.get_window_extent(fig.canvas.get_renderer())
        fits = bb.x0 >= 0 and bb.x1 <= fig_w * fig.dpi
        if fits:
            break
        lgd.remove()
    if not fits:      # no width fits: keep the one-column legend and print a warning
        nc, lfs = 1, fs.LEGEND_FS
        lgd = _place_legend("upper center", (0.5, 1.0), nc, lfs)
        fig.canvas.draw()
        bb = lgd.get_window_extent(fig.canvas.get_renderer())
        print("  legend is wider than the figure even at one column; shorten the "
              "contrast labels")
    legends = [lgd]
    # the key at the very bottom, under the last row of module labels
    leg_h = (bb.y1 - bb.y0) / fig.dpi + 0.25 if legend else 0.0
    # Columns after the first carry tick labels but no y label, so their left
    # margin is l_inner (the tick labels, 7 mm printed).
    l_inner = 1.25
    ax_w = (fig_w - l_in - (ncols - 1) * l_inner - ncols * r_in) / ncols
    ax_h = (fig_h - leg_h - nrows * (t_in + b_in)) / nrows

    panels = []      # (organelle, ax, bar labels)
    for k, organelle in enumerate(organelles):
        row, col = divmod(k, ncols)
        top = t_in + row * (ax_h + b_in + t_in)     # axes top, from the figure top
        ax = fig.add_axes([(l_in + col * (ax_w + r_in + l_inner)) / fig_w,
                           1 - (top + ax_h) / fig_h,
                           ax_w / fig_w, ax_h / fig_h])

        lo_ax, hi_ax = _lim(organelle)
        span = hi_ax - lo_ax
        mods = modules[modules.organelle == organelle]
        bar_labels = []
        x = np.arange(len(mods))
        w = 0.8 / len(contrasts)
        # a rotated count label is as wide on the page as the font's line height,
        # so its size follows from the bar slot, down to the 5 pt minimum
        slot_in = ax_w / (len(mods) * len(contrasts))
        annot_fs = min(fs.ANNOT_FS, max(fs.native_pt(5.0), slot_in * 60 * 0.85))
        ax.axhline(0, color="black", lw=fs.AXES_LW, zorder=1)

        for j, cname in enumerate(contrasts):
            for i, m in enumerate(mods.itertuples()):
                b = bars.loc[(cname, m.go_id)]
                xc = x[i] - 0.4 + w * (j + 0.5)
                if b.n_members < min_members:
                    # one marker per missing bar: a blank slot reads as zero. On the
                    # zero line, above it where the axis has room for the text,
                    # hanging below it otherwise.
                    txt = f"n<{min_members}"
                    up = hi_ax / span * ax_h * 72 >= (len(txt) * 0.6 + 1) * annot_fs
                    ax.text(xc, (0.025 if up else -0.025) * span, txt, ha="center",
                            va="bottom" if up else "top", rotation=90,
                            fontsize=annot_fs, color="0.35", zorder=4)
                    continue
                h = float(b.delta_module)
                ax.bar(xc, h, width=w * 0.88, color=colors[cname], edgecolor="none",
                       zorder=2)
                if abs(h) < zero_eps:
                    ax.plot([xc - w * 0.44, xc + w * 0.44], [0, 0], color=colors[cname],
                            lw=fs.LINE_LW, solid_capstyle="butt", zorder=3)
                pad_y, up = 0.025 * span, h >= 0
                is_inside = abs(h) >= inside * max(hi_ax, -lo_ax)
                if is_inside:
                    # dark fills need light type
                    lum = sum(to_rgb(colors[cname])) / 3
                    ty = h - (pad_y if up else -pad_y)
                    va = "top" if up else "bottom"
                    tcol = "white" if lum < 0.55 else "0.15"
                else:
                    ty = (h if abs(h) >= zero_eps else 0) + (pad_y if up else -pad_y)
                    va, tcol = ("bottom" if up else "top"), "0.25"
                bar_labels.append((ax.text(xc, ty, f"{int(b.n_members)}", ha="center",
                                           va=va, rotation=90, fontsize=annot_fs,
                                           color=tcol, zorder=4),
                                   (cname, m.go_id) in starred, up, xc,
                                   annot_fs, cname, h, is_inside))

        ax.set_xlim(-0.5, len(mods) - 0.5)
        ax.set_ylim(lo_ax, hi_ax)
        # ticks inside the limits, zero always among them
        step = next(v for v in (1.0, 0.5, 0.25, 0.2, 0.1, 0.05) if v <= span / 4)
        ax.set_yticks(np.round(np.arange(np.ceil(lo_ax / step) * step,
                                         hi_ax + 1e-9, step), 3))
        ax.set_xticks(x)
        # wrapped at 7 characters: at 9, "U2 snRNP" ran into its neighbour at 6 pt
        mem = ("GO:0005741", "GO:0005743")      # outer, inner mitochondrial membrane
        tick = {**short, mem[0]: "outer", mem[1]: "inner"}
        ax.set_xticklabels(["\n".join(textwrap.wrap(_sentence_case(tick[m.go_id]), 7,
                                                     break_long_words=False))
                            for m in mods.itertuples()])
        # The two membranes share one "membrane", centred under the pair on the
        # second label line. It is a minor tick label because the overlap check
        # below reads major labels only; minor ticks get the major size and pad.
        pair = [i for i, g in enumerate(mods.go_id) if g in mem]
        if pair:
            assert pair == [pair[0], pair[0] + 1], \
                f"{organelle}: outer and inner membrane are not neighbours"
            ax.set_xticks([x[pair].mean()], ["\nmembrane"], minor=True)
            ax.tick_params(axis="x", which="minor", color="none",
                           size=plt.rcParams["xtick.major.size"],
                           pad=plt.rcParams["xtick.major.pad"])
        if col == 0:
            ax.set_ylabel("log2 FC")
        ax.set_title(organelle[0].upper() + organelle[1:], pad=12)
        panels.append((organelle, ax, bar_labels))

    # --- the key at the very bottom of the figure ------------------------------
    # removed and drawn again at the bottom: matplotlib has no public way to
    # change the loc of a legend
    lgd.remove()
    if not legend:
        legends = []
    else:
        lgd = _place_legend("lower center", (0.5, 0.0), nc, lfs)
        legends = [lgd]

    # --- stars, placed against the rendered count labels -----------------------
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    stars = []
    for organelle, ax, bar_labels in panels:
        inv = ax.transData.inverted()
        sp = ax.get_ylim()[1] - ax.get_ylim()[0]
        for txt, sig, up, xc, size, cname, h, is_inside in bar_labels:
            if not sig:
                continue
            # past the bar end when the count sits inside it, past the count otherwise
            if is_inside:
                y = h
            else:
                bb = txt.get_window_extent(rend)
                y = inv.transform((0, bb.y1 if up else bb.y0))[1]
            # The star is a marker path, not the character U+2605: an SVG renderer
            # without DejaVu Sans draws the glyph as an empty box.
            tr = offset_copy(ax.transData, fig=fig, y=(1 if up else -1) * star_lift * size,
                             units="points")
            stars.append(ax.plot([xc], [y + (0.015 if up else -0.015) * sp],
                                 marker="*", ms=star_ms * size, ls="none",
                                 color=colors[cname], mec="none", transform=tr,
                                 clip_on=False, zorder=5)[0])

    # --- checks on the rendered canvas, printed ----------------------------------
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    W, H = fig.get_size_inches() * fig.dpi
    n_lab, problems = 0, []
    for organelle, ax, bar_labels in panels:
        for what, arts in (("x tick labels", ax.get_xticklabels()),
                           ("count labels", [t for t, *_ in bar_labels])):
            ov = [(a.get_text(), b.get_text()) for i, a in enumerate(arts)
                  for b in arts[i + 1:]
                  if a.get_window_extent(rend).overlaps(b.get_window_extent(rend))]
            if ov:
                problems.append(f"{organelle}: {what} overlap, first {ov[0]}")
        n_lab += len(bar_labels)
    for organelle, ax, _ in panels:
        ticks = [t.get_window_extent(rend) for t in ax.get_xticklabels()]
        for s in [s for s in stars if s.axes is ax]:
            bb = s.get_window_extent(rend)
            if bb.y1 > H + 1 or bb.y0 < -1:
                problems.append(f"{organelle}: a star falls off the canvas; raise the pad")
            if any(bb.overlaps(t) for t in ticks):
                problems.append(f"{organelle}: a star runs into a module label; raise the pad")
    for organelle, ax, _ in panels:
        tt = ax.title.get_window_extent(rend)
        if any(tt.overlaps(lg.get_window_extent(rend)) for lg in legends):
            problems.append(f"{organelle}: the title and the legend overlap")
    print(f"CHECK {len(panels)} panels, {n_lab} bar labels, {len(stars)} stars, "
          + (f"legend on {int(np.ceil(len(handles) / nc))} row(s), " if legend
             else "no legend, ")
          + "axis " + ", ".join(f"{a.get_ylim()[0]:+.2f} to {a.get_ylim()[1]:+.2f}"
                                for _, a, _ in panels) + ", count labels "
          f"{min(l[4] for _, _, ls in panels for l in ls) / fs.DRAW_SCALE:.1f} pt")
    for p in problems:
        print(f"  !! {p}")

    fs.save_panel(fig, letter, name, folder)
    return fig


def split_dumbbell(rows, contrasts, colors, missing_labels, legend_labels, row_labels,
                   letter, name, folder, size=(2, 1), ms_sig_pt=3.0, ms_ns_pt=2.6,
                   offset=0.24):
    """Draw the dumbbell panel and save it.

    rows            DataFrame, one row per pair, top to bottom; per contrast c the
                    columns c__measured, c__median_spared, c__median_lost, c__sig, c__q
    contrasts       the contrasts, one dumbbell each per row, the first above
    colors          per contrast: line colour of a significant dumbbell and its q value
    missing_labels  per contrast: name in "not measurable in ..."
    legend_labels   per contrast: legend entry of its line colour
    row_labels      one label per row
    letter, name, folder   passed to save_panel
    size            panel size in units (columns, rows)
    ms_sig_pt, ms_ns_pt    printed marker diameters: filled = significant, open = not
    offset          vertical distance of each dumbbell from the row centre, in rows;
                    two q labels of 6 pt need 2 x 0.24 rows
    """
    # The part with the higher log2 FC in an aggregate-vs-host contrast is the one
    # enriched in the aggregate, so UP_COLOR and DOWN_COLOR keep their meaning.
    enr_c, dep_c = fs.UP_COLOR, fs.DOWN_COLOR
    ns_line, ns_text, zero_color = fs.BG_COLOR, "0.45", "black"
    band_color = fs.shade(fs.BG_COLOR, 0.72)
    ms_sig, ms_ns = fs.native_pt(ms_sig_pt), fs.native_pt(ms_ns_pt)
    n = len(rows)

    fig, axm = create_figure(1, 1, figsize=fs.figsize(*size))
    ax = axm.next()
    ax.grid(False)

    v = pd.concat([rows[f"{c}__median_{p}"] for c in contrasts for p in ("spared", "lost")
                   if f"{c}__median_{p}" in rows]).dropna()
    pad = max(0.35, 0.18 * (v.max() - v.min()))
    xlo, xhi = v.min() - pad * 2.6, v.max() + pad      # left room for "not measurable"
    ys = np.arange(n)[::-1]

    for i, (_, r) in enumerate(rows.iterrows()):
        yb = ys[i]
        if i % 2:
            ax.axhspan(yb - 0.5, yb + 0.5, color=band_color, zorder=0)
        for k, c in enumerate(contrasts):
            yy = yb + (offset if k == 0 else -offset)
            if not r[f"{c}__measured"]:
                ax.text(xlo + 0.04 * (xhi - xlo), yy, f"not measurable in {missing_labels[c]}",
                        va="center", ha="left", fontsize=fs.ANNOT_FS, color=ns_text,
                        style="italic", zorder=4)
                continue
            ms, ml = r[f"{c}__median_spared"], r[f"{c}__median_lost"]
            sig = bool(r[f"{c}__sig"])
            ax.plot([ml, ms], [yy, yy], color=colors[c] if sig else ns_line,
                    lw=fs.LINE_LW if sig else fs.THIN_LW, solid_capstyle="round", zorder=2)
            for xv, col in ((ms, enr_c), (ml, dep_c)):
                if sig:
                    ax.plot(xv, yy, "o", ms=ms_sig, color=col, mec="black",
                            mew=fs.MARKER_LW, zorder=4)
                else:
                    ax.plot(xv, yy, "o", ms=ms_ns, color="white", mec=col,
                            mew=fs.LINE_LW, zorder=4)
            ax.text(1.02, yy, f"q = {r[f'{c}__q']:.3f}", va="center", ha="left",
                    fontsize=fs.ANNOT_FS, color=colors[c] if sig else ns_text,
                    fontweight="bold" if sig else "normal", clip_on=False,
                    transform=blended_transform_factory(ax.transAxes, ax.transData))

    ax.axvline(0, color=zero_color, lw=fs.AXES_LW, zorder=3)
    ax.grid(axis="x", color=fs.BG_COLOR, lw=fs.THIN_LW, alpha=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(-0.55, n - 0.45)
    ax.set_yticks(ys)
    ax.set_yticklabels(row_labels)
    ax.set_xlabel("median log2 FC, aggregate vs host")

    # One row: the two marker colours and one line colour per contrast. The figure
    # legend text carries the rest (grey line, open marker, grey q = not significant).
    handles = [
        Line2D([], [], marker="o", ls="none", ms=ms_sig, color=enr_c, mec="black",
               mew=fs.MARKER_LW, label="higher"),
        Line2D([], [], marker="o", ls="none", ms=ms_sig, color=dep_c, mec="black",
               mew=fs.MARKER_LW, label="lower"),
    ] + [Line2D([], [], color=colors[c], lw=fs.LINE_LW, label=legend_labels[c])
         for c in contrasts]
    fig.legend(handles=handles, loc="outside lower center", ncol=len(handles),
               frameon=False, borderaxespad=0, handlelength=1.0, handletextpad=0.3,
               columnspacing=0.8)

    fs.save_panel(fig, letter, name, folder)


def property_levels(table, levels, contrasts, order, prop_labels, colors, contrast_labels,
                    letter, name, folder, size=(2, 1), ms_pt=2.5, offset=0.18,
                    x_title="depleted ←   rank correlation of property with fold change"
                            "   → enriched"):
    """Draw one subplot per level and save the panel.

    table           long table with the columns contrast, level, prop, stat, q
    levels          levels in subplot order; "L1 subunits" is titled "L1" over
                    "subunits"
    contrasts       the contrasts drawn, one marker per property each, the first above
    order           properties, top to bottom
    prop_labels     property -> row label
    colors          per contrast: line and marker colour
    contrast_labels per contrast: legend entry
    letter, name, folder   passed to save_panel
    size            panel size in units (columns, rows)
    ms_pt           printed marker diameter
    offset          vertical distance of each contrast from the row centre, in rows

    Every subplot has the same rows, x limits and ticks. Marker fill is the only
    encoding of significance (q < 0.05).
    """
    tint_pos, tint_neg = fs.shade(fs.UP_COLOR, 0.88), fs.shade(fs.DOWN_COLOR, 0.88)
    ms = fs.native_pt(ms_pt)
    d = table[table.contrast.isin(contrasts)]
    labels = [prop_labels.get(p, p) for p in order]
    n = len(order)

    # shared x limits, and three labelled ticks: more do not fit four narrow subplots
    xlim = float(np.abs(d.stat).max()) * 1.18
    step = next((v for v in (1.0, 0.5, 0.4, 0.25, 0.2, 0.1, 0.05, 0.02, 0.01)
                 if v <= xlim * 0.85), round(xlim * 0.85, 3))
    xticks = [-step, 0.0, step]

    # one legend; the fill entries in neutral grey, because fill encodes significance
    # for every contrast and a coloured example would read as a third contrast
    handles = [Line2D([], [], marker="o", ls="-", ms=ms, lw=fs.LINE_LW, color=colors[c],
                      label=contrast_labels[c]) for c in contrasts]
    handles += [Line2D([], [], marker="o", ls="none", ms=ms, color="0.25", mec="0.25",
                       mew=fs.LINE_LW, label="q < 0.05"),
                Line2D([], [], marker="o", ls="none", ms=ms, color="white", mec="0.25",
                       mew=fs.LINE_LW, label="not significant")]

    fig, axm = create_figure(1, len(levels), figsize=fs.figsize(*size),
                             subplots_kwargs={"sharey": True})
    ys = np.arange(n)[::-1]
    axes = []
    for k, level in enumerate(levels):
        ax = axm.next()
        ax.grid(False)
        axes.append(ax)
        ax.axvspan(0, xlim, color=tint_pos, zorder=0)
        ax.axvspan(-xlim, 0, color=tint_neg, zorder=0)
        for i, prop in enumerate(order):
            for j, c in enumerate(contrasts):
                r = d[(d.contrast == c) & (d.level == level) & (d["prop"] == prop)]
                if not len(r):
                    continue
                val, q = float(r.stat.iloc[0]), float(r.q.iloc[0])
                yy = ys[i] + (offset if j == 0 else -offset)
                ax.plot([0, val], [yy, yy], color=colors[c], lw=fs.LINE_LW, zorder=2)
                ax.plot(val, yy, "o", ms=ms, zorder=4,
                        color=colors[c] if q < 0.05 else "white", mec=colors[c],
                        mew=fs.LINE_LW)
        ax.axvline(0, color="black", lw=fs.AXES_LW, zorder=3)
        ax.set_xlim(-xlim, xlim)
        ax.set_ylim(-0.6, n - 0.4)
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"{t:g}".replace("-", "−") for t in xticks])   # "0", not "0.0"
        ax.set_yticks(ys)
        if k == 0:
            ax.set_yticklabels(labels)
        else:   # the y axis is shared: hiding the labels here keeps them on the first
            ax.tick_params(axis="y", labelleft=False, length=0)
        num, unit = level.split(" ", 1)
        ax.set_title(f"{num}\n{unit}")
        ax.grid(axis="x", color=fs.BG_COLOR, lw=fs.THIN_LW, alpha=0.6, zorder=0)
        ax.set_axisbelow(True)

    fig.supxlabel(x_title, fontsize=fs.LABEL_FS)
    # above the subplots, in one row: below them it shares the bottom edge with
    # the x title, and in two rows it covers the property labels
    fig.legend(handles=handles, loc="outside upper center", ncol=4, frameon=False,
               borderaxespad=0, handlelength=1.2, handletextpad=0.3, columnspacing=1.0)

    fig.canvas.draw()
    assert tuple(t.get_text() for t in axes[0].get_yticklabels()) == tuple(labels), \
        "leftmost subplot labels are not the row order"
    assert all(not t.get_visible() for a in axes[1:] for t in a.get_yticklabels()), \
        "a subplot other than the first shows property labels"
    assert len({a.get_xlim() for a in axes}) == 1, "x limits differ between subplots"
    assert len({tuple(a.get_xticks()) for a in axes}) == 1, "x ticks differ between subplots"
    print(f"CHECK {len(axes)} subplots, {n} properties, x limits {axes[0].get_xlim()}, "
          f"ticks {xticks}, labels on the first subplot only")

    fs.save_panel(fig, letter, name, folder)


def split_property_heatmap(values, row_labels, headers, cbar_label, letter, name, folder,
                           size=(2, 1)):
    """Draw the heatmap and save it.

    values      DataFrame, one row per pair (top to bottom) and one column per
                property; NaN draws an empty cell
    row_labels  one label per row, on one line (two lines need twice the row height)
    headers     property -> column header, drawn above the grid
    cbar_label  colour bar label
    letter, name, folder   passed to save_panel
    size        panel size in units (columns, rows)

    The colour scale is symmetric at the largest |value|, rounded up to 0.1. A cell
    carries "+" or "−" from the same number, white on dark cells.
    """
    props = list(values.columns)
    n = len(values)
    fig, axm = create_figure(1, 1, figsize=fs.figsize(*size))
    ax = axm.next()
    ax.grid(False)
    ys = np.arange(n)[::-1]

    vmax = float(np.ceil(np.nanmax(np.abs(values.to_numpy(dtype=float))) * 10) / 10) or 1.0
    glyphs = []
    for i in range(n):
        y = ys[i]
        if i % 2:
            ax.axhspan(y - 0.5, y + 0.5, color=fs.shade(fs.BG_COLOR, 0.82), zorder=0)
        for j, p in enumerate(props):
            v = values.iloc[i][p]
            if not np.isfinite(v):
                ax.add_patch(Rectangle((j, y - 0.5), 1, 1, facecolor="white",
                                       edgecolor="0.85", lw=fs.THIN_LW, zorder=2))
                continue
            ax.add_patch(Rectangle((j, y - 0.5), 1, 1, zorder=2, lw=fs.THIN_LW,
                                   edgecolor="white",
                                   facecolor=fs.DIRECTION_CMAP(
                                       0.5 + 0.5 * np.clip(v / vmax, -1, 1))))
            g = "+" if v > 0 else ("−" if v < 0 else "0")
            ax.text(j + 0.5, y, g, va="center", ha="center", fontsize=fs.ANNOT_FS,
                    fontweight="bold", zorder=5,
                    color="white" if abs(v) > 0.55 * vmax else "black")
            glyphs.append((v, g))

    ax.set_xlim(0, len(props))
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_xticks(np.arange(len(props)) + 0.5)
    ax.set_xticklabels([headers[p] for p in props])
    ax.xaxis.set_ticks_position("top")
    ax.set_yticks(ys)
    ax.set_yticklabels(row_labels)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    # colour bar below the grid, placed by constrained layout
    cb = fig.colorbar(ScalarMappable(norm=Normalize(-vmax, vmax), cmap=fs.DIRECTION_CMAP),
                      ax=ax, orientation="horizontal", location="bottom",
                      shrink=0.45, aspect=30, pad=0.02)
    cb.set_ticks([-vmax, 0.0, vmax])
    cb.set_label(cbar_label)

    fig.canvas.draw()
    assert [t.get_text() for t in ax.get_yticklabels()] == list(row_labels)
    for v, g in glyphs:
        assert (v > 0 and g == "+") or (v < 0 and g == "−") or (v == 0 and g == "0"), \
            f"glyph {g} disagrees with its own value {v:+.3f}"
    print(f"CHECK {n} rows in the given order; all {len(glyphs)} glyphs agree in sign "
          f"with their own cell value")
    # column headers against the column width, one-line row labels against the row
    rend = fig.canvas.get_renderer()
    col_px = ax.get_window_extent().width / len(props)
    widths = []
    for p in props:
        probes = [fig.text(0, 0, line, fontsize=fs.TICK_FS) for line in headers[p].split("\n")]
        widths.append(max(t.get_window_extent(rend).width for t in probes))
        for t in probes:
            t.remove()
    print(f"CHECK header widths vs {col_px:.0f}px column: {[f'{w:.0f}' for w in widths]}")
    lh = max(t.get_window_extent().height for t in ax.get_yticklabels())
    row_px = ax.get_window_extent().height / n
    print(f"CHECK row label height {lh:.0f}px in a {row_px:.0f}px row")
    if lh > row_px:
        print("  !! a row label is taller than its row; shorten the names or drop rows")

    fs.save_panel(fig, letter, name, folder)


def concept_pair_heatmap(lfc, annot, concepts, subtitles, column_labels, row_groups,
                         labelled, letter, name, folder, size=(2, 2), norm=None,
                         group_gap_pt=1.5):
    """Draw the heatmap and save it.

    lfc            DataFrame, genes (rows, top to bottom) x columns; two columns per
                   concept, in concept order
    annot          DataFrame like lfc with the text of each cell (e.g. "*")
    concepts       concept names, one per column pair; a space starts a new line
    subtitles      concept -> parenthetical, drawn once over each run of neighbouring
                   concepts that share it
    column_labels  one label per column
    row_groups     list of (group name, genes), in row order; a solid line separates
                   the groups
    labelled       group names that get a rotated label left of the gene names
    letter, name, folder   passed to save_panel
    size           panel size in units (columns, rows)
    norm           colour normalisation
    group_gap_pt   printed gap from the gene names to the group labels
    """
    group_gap = fs.native_pt(group_gap_pt)
    fig, axm = create_figure(1, 1, figsize=fs.figsize(*size))
    ax = axm.next()
    ax.grid(False)

    sns.heatmap(lfc, ax=ax, cmap=fs.DIRECTION_CMAP, norm=norm, linewidths=0,
                cbar_kws={"shrink": 0.4, "pad": 0.02}, annot=annot, fmt="",
                annot_kws={"fontsize": fs.ANNOT_FS, "color": "black"})
    ax.set_xticklabels(column_labels, rotation=90, ha="center", fontsize=fs.TICK_FS)
    ax.set_yticks(np.arange(len(lfc.index)) + 0.5)
    ax.set_yticklabels(lfc.index, rotation=0, fontsize=fs.TICK_FS)
    cbar = ax.collections[0].colorbar
    cbar.set_label("log2 fold change", fontsize=fs.LABEL_FS)
    cbar.ax.tick_params(labelsize=fs.TICK_FS)

    for i in range(1, len(concepts)):
        ax.axvline(i * 2, color="black", linewidth=fs.AXES_LW)

    # Concept headers on two lines over their column pair, the italic parenthetical
    # once over each run of concepts that share it. Offsets are in points, so the
    # two rows cannot overlap whatever the axes height.
    trans = ax.get_xaxis_transform()
    for i, c in enumerate(concepts):
        ax.annotate(c.replace(" ", "\n"), xy=(i * 2 + 1, 1.0), xycoords=trans,
                    xytext=(0, 8 + fs.TITLE_FS + 4), textcoords="offset points",
                    ha="center", va="bottom", annotation_clip=False, fontsize=fs.TITLE_FS)
    runs = []   # [parenthetical, first concept, last concept]
    for i, c in enumerate(concepts):
        if runs and subtitles[c] == runs[-1][0]:
            runs[-1][2] = i
        else:
            runs.append([subtitles[c], i, i])
    for sub, first, last in runs:
        ax.annotate(sub, xy=(first + last + 1, 1.0), xycoords=trans,   # centre of the run
                    xytext=(0, 8), textcoords="offset points", ha="center", va="bottom",
                    annotation_clip=False, fontsize=fs.TITLE_FS, fontstyle="italic",
                    color="grey")

    # row groups: a solid line between groups, labels for the labelled ones
    starts, running = [], 0
    for group, genes in row_groups:
        n = sum(g in lfc.index for g in genes)
        if n:
            starts.append((group, running, n))
            running += n
    for _, start, _ in starts[1:]:
        ax.axhline(start, color="black", linewidth=fs.AXES_LW)
    ax2 = ax.secondary_yaxis("left")
    ax2.set_yticks([start + n / 2 for group, start, n in starts if group in labelled])
    ax2.set_yticklabels([group for group, _, _ in starts if group in labelled],
                        fontsize=fs.TICK_FS, rotation=90, va="center", ha="center")
    ax2.tick_params(length=0)
    for spine in ax2.spines.values():
        spine.set_visible(False)

    # seaborn draws no cell for a missing value, so the axes background shows:
    # grey, because white is log2 FC = 0 on this colour map
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_facecolor(fs.BG_COLOR)
    print(f"CHECK {int(lfc.isna().sum().sum())} of {lfc.size} cells have no value "
          f"(drawn grey)")

    # group labels start where the drawn gene names end; a label is anchored at
    # its centre
    fig.canvas.draw()
    rend, to_pt = fig.canvas.get_renderer(), 72 / fig.dpi
    w_gene = (max(t.get_window_extent(rend).width for t in ax.get_yticklabels()) * to_pt
              + ax.yaxis.get_major_ticks()[0].get_tick_padding()
              + plt.rcParams["ytick.major.pad"])
    w_grp = max(t.get_window_extent(rend).width for t in ax2.get_yticklabels()
                if t.get_text()) * to_pt
    ax2.tick_params(pad=w_gene + group_gap + w_grp / 2)
    fig.canvas.draw()   # one more layout pass with the final pad before save_panel

    fs.save_panel(fig, letter, name, folder)
