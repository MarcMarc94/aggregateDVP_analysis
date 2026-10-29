"""Figure style shared by the aggregate DVP analysis notebooks.

The single definition of the page grid, the fonts and line widths, the colour
roles, the matplotlib and alphapepttools settings, save_panel and
compose_figure. Each notebook keeps only what belongs to its experimental
system: the condition order, colours and short labels, and its figure folders.

Nature: a figure is at most 180 x 170 mm (full page) or half that width
(single column); text 5-7 pt; lines 0.25-1 pt; panel letters 8 pt bold
lowercase; 300 dpi RGB.

Every size here is stated as it must PRINT. matplotlib draws a panel at
DRAW_SCALE times its printed size, so the value set in matplotlib is the
printed value times DRAW_SCALE; native_pt() does that conversion.

Importing this module sets matplotlib's rcParams and alphapepttools' config.
A notebook with src/ on sys.path imports the names it uses, e.g.

    from figure_style import figsize, native_pt, save_panel, compose_figure
"""
import io
import re
from datetime import date
from pathlib import Path

import cairosvg
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import svgutils.compose as svgutils
import svgutils.transform as st
from IPython.display import SVG as _SVG, display as _display
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.text import Text
from PIL import Image

import alphapepttools as at
from alphapepttools.pl import figure as _at_fig, plots as _at_plt, colors as _at_col

__all__ = [
    # page grid
    "MM_PER_IN", "DRAW_SCALE", "SCALE", "PAGE_W_MM", "PAGE_H_MM", "MARGIN_MM", "GAP_MM",
    "GRID_COLS", "GRID_ROWS", "UNIT_W_MM", "UNIT_H_MM", "HALF_W_MM",
    "figsize", "native_pt",
    # limits
    "TEXT_MIN_PT", "TEXT_MAX_PT", "STROKE_MIN_PT", "STROKE_MAX_PT", "LABEL_PT", "DPI", "SCALE_TOL",
    # fonts and lines
    "ANNOT_FS", "TICK_FS", "LEGEND_FS", "LABEL_FS", "TITLE_FS",
    "THIN_LW", "MARKER_LW", "AXES_LW", "LINE_LW", "PANEL_FONTS", "PANEL_STROKES",
    # colour roles
    "UP_COLOR", "DOWN_COLOR", "DIRECTION_CMAP", "BG_COLOR", "HIGHLIGHT_COLOR", "REF_COLOR",
    "shade",
    # functions
    "style_summary", "audit_panel", "save_panel", "compose_figure",
]

# ---------------------------------------------------------------- PAGE GRID
# A full page is 4 x 4 units, a single-column figure 2 x 4 units. The unit
# follows from the page: 180 mm = margin + 4 units + 3 gaps, and the same for
# 170 mm in height. The margin is the strip left of and above the panels that
# holds the panel letters.
MM_PER_IN  = 25.4
DRAW_SCALE = 4.5                  # matplotlib canvas = 4.5 x the printed size
SCALE      = 1 / DRAW_SCALE       # printed size / matplotlib size

PAGE_W_MM, PAGE_H_MM = 180.0, 170.0
MARGIN_MM, GAP_MM    = 3.0, 3.0
GRID_COLS, GRID_ROWS = 4, 4

UNIT_W_MM = (PAGE_W_MM - MARGIN_MM - (GRID_COLS - 1) * GAP_MM) / GRID_COLS   # 42.0
UNIT_H_MM = (PAGE_H_MM - MARGIN_MM - (GRID_ROWS - 1) * GAP_MM) / GRID_ROWS   # 39.5
HALF_W_MM = MARGIN_MM + 2 * UNIT_W_MM + GAP_MM                               # 90.0


def figsize(cols=1, rows=1):
    """Canvas in inches for a panel cols x rows units big, inner gaps included."""
    w_mm = cols * UNIT_W_MM + (cols - 1) * GAP_MM
    h_mm = rows * UNIT_H_MM + (rows - 1) * GAP_MM
    return (w_mm * DRAW_SCALE / MM_PER_IN, h_mm * DRAW_SCALE / MM_PER_IN)


def native_pt(printed_pt):
    """Value to set in matplotlib so that it prints at printed_pt."""
    return printed_pt * DRAW_SCALE


# ---------------------------------------------------------------- LIMITS
TEXT_MIN_PT, TEXT_MAX_PT     = 5.0, 7.0
STROKE_MIN_PT, STROKE_MAX_PT = 0.25, 1.0
LABEL_PT  = 8.0       # panel letters
DPI       = 300
SCALE_TOL = 0.05      # largest rescale compose_figure applies without a warning

# ---------------------------------------------------------------- FONTS
# Nature: all text 5-7 pt, panel letters 8 pt. Set near the middle, 6.0-6.5 pt,
# so a panel rescaled by SCALE_TOL still prints at 5.7-6.8 pt.
ANNOT_FS  = native_pt(6.0)    # text inside the axes
TICK_FS   = native_pt(6.0)
LEGEND_FS = native_pt(6.0)
LABEL_FS  = native_pt(6.5)    # axis labels
TITLE_FS  = native_pt(6.5)

# ---------------------------------------------------------------- LINES
# Nature: every line and stroke 0.25-1 pt. Set at 0.4-0.8 pt.
THIN_LW   = native_pt(0.4)    # ticks, grid lines, separators
MARKER_LW = native_pt(0.4)    # edges of markers, bars, boxes and violins
AXES_LW   = native_pt(0.6)    # axis lines, reference lines
LINE_LW   = native_pt(0.8)    # data lines

_font = "Arial" if "Arial" in {f.name for f in mpl.font_manager.fontManager.ttflist} else "Helvetica"

plt.rcParams.update({
    "svg.fonttype":          "none",            # text stays editable text in the SVG
    "font.family":           "sans-serif",
    "font.sans-serif":       [_font, "Helvetica", "DejaVu Sans"],
    "axes.titlesize":        TITLE_FS,
    "axes.labelsize":        LABEL_FS,
    "xtick.labelsize":       TICK_FS,
    "ytick.labelsize":       TICK_FS,
    "legend.fontsize":       LEGEND_FS,
    "legend.title_fontsize": LEGEND_FS,
    # Legend layout, in units of the legend font size. frameon, borderpad,
    # handlelength, handleheight, handletextpad and columnspacing take the values
    # that legendkit sets when decoupler imports it; labelspacing and borderaxespad
    # are matplotlib's defaults. Stating them here gives every notebook the same
    # legends in any import order.
    "legend.frameon":        False,
    "legend.borderpad":      0.0,
    "legend.labelspacing":   0.5,
    "legend.handlelength":   1.0,
    "legend.handleheight":   1.0,
    "legend.handletextpad":  0.5,
    "legend.borderaxespad":  0.5,
    "legend.columnspacing":  1.0,
    # matplotlib's default lengths in points assume a figure drawn at print
    # size; at DRAW_SCALE a 3.5 pt tick prints 0.27 mm long. Stated as printed:
    "axes.titlepad":         native_pt(3.0),
    "axes.labelpad":         native_pt(2.0),
    "xtick.major.size":      native_pt(2.5),
    "ytick.major.size":      native_pt(2.5),
    "xtick.minor.size":      native_pt(1.5),
    "ytick.minor.size":      native_pt(1.5),
    "xtick.major.pad":       native_pt(1.5),
    "ytick.major.pad":       native_pt(1.5),
    "axes.linewidth":        AXES_LW,
    "xtick.major.width":     THIN_LW,
    "ytick.major.width":     THIN_LW,
    "xtick.minor.width":     THIN_LW,
    "ytick.minor.width":     THIN_LW,
    "lines.linewidth":       LINE_LW,
    "lines.markeredgewidth": MARKER_LW,
    "patch.linewidth":       MARKER_LW,
    "grid.linewidth":        THIN_LW,
    # Light, so that a 0.4 pt grid stays behind the data.
    "grid.color":            "0.95",
    "hatch.linewidth":       THIN_LW,
    # box plots take their widths from their own settings, not from patch.linewidth;
    # the default median colour is orange, which is a condition colour in the
    # HEK notebook (CW)
    "boxplot.boxprops.linewidth":     MARKER_LW,
    "boxplot.whiskerprops.linewidth": MARKER_LW,
    "boxplot.capprops.linewidth":     MARKER_LW,
    "boxplot.medianprops.linewidth":  LINE_LW,
    "boxplot.medianprops.color":      "black",
    "axes.spines.top":       False,
    "axes.spines.right":     False,
    # a saved panel is exactly its canvas, never cropped to the drawn content
    "savefig.bbox":          "standard",
    "savefig.transparent":   False,
    "savefig.facecolor":     "white",
    "savefig.dpi":           DPI,
})

# alphapepttools keeps its own sizes. stylize(), which axm.next() calls on every
# axes, sets tick and grid widths from linewidths["small"]; the plot functions
# use "medium" for reference lines and "large" for box, violin and bar outlines.
# At their defaults (0.25 / 0.5 / 1.25) these print at 0.06 / 0.11 / 0.28 pt.
for _m in (_at_fig, _at_plt, _at_col):
    _m.config["axes"].update(title_size=TITLE_FS, label_size=LABEL_FS, tick_size=TICK_FS)
    _m.config["legend"].update(font_size=LEGEND_FS, title_size=LEGEND_FS)
    _m.config["linewidths"].update(small=THIN_LW, medium=AXES_LW, large=MARKER_LW)
    _m.config["font_sizes"].update(small=ANNOT_FS, medium=TICK_FS, large=LABEL_FS)

# ---------------------------------------------------------------- COLOUR ROLES
# The condition colours belong to each notebook. These roles are the same in all.

# Direction of change: red up, blue down, discrete and continuous. Never used
# in a panel together with the condition colours.
UP_COLOR, DOWN_COLOR = "#b2182b", "#2166ac"
DIRECTION_CMAP = LinearSegmentedColormap.from_list("direction", [DOWN_COLOR, "white", UP_COLOR])

# Other roles. Call-outs differ by lightness, because any hue would collide
# with a condition colour or with up/down.
BG_COLOR        = "#d0d0d0"   # background points and lines
HIGHLIGHT_COLOR = "#882255"   # call-outs only
REF_COLOR       = "0.5"       # label leader lines and reference lines (zero, thresholds)


def shade(color, f):
    """`color` blended toward white by the fraction f (0 = unchanged, 1 = white), as RGB."""
    r, g, b = mcolors.to_rgb(color)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


# ---------------------------------------------------------------- CHECKS
_fonts   = {"ANNOT": ANNOT_FS, "TICK": TICK_FS, "LEGEND": LEGEND_FS,
            "LABEL": LABEL_FS, "TITLE": TITLE_FS}
_strokes = {"THIN": THIN_LW, "MARKER": MARKER_LW, "AXES": AXES_LW, "LINE": LINE_LW}
for _n, _v in _fonts.items():       # must stay inside the limits after a SCALE_TOL rescale
    _lo, _hi = _v * SCALE * (1 - SCALE_TOL), _v * SCALE * (1 + SCALE_TOL)
    assert TEXT_MIN_PT <= _lo and _hi <= TEXT_MAX_PT, f"{_n}_FS prints at {_lo:.2f}-{_hi:.2f} pt"
for _n, _v in _strokes.items():
    _lo, _hi = _v * SCALE * (1 - SCALE_TOL), _v * SCALE * (1 + SCALE_TOL)
    assert STROKE_MIN_PT <= _lo and _hi <= STROKE_MAX_PT, f"{_n}_LW prints at {_lo:.2f}-{_hi:.2f} pt"

PANEL_FONTS   = (min(_fonts.values()), max(_fonts.values()))      # read by compose_figure
PANEL_STROKES = (min(_strokes.values()), max(_strokes.values()))


def style_summary():
    """Typeface, unit, and the printed font and line sizes, as text."""
    return "\n".join([
        f"typeface {_font}",
        f"unit {UNIT_W_MM:.1f} x {UNIT_H_MM:.1f} mm, canvas {figsize()[0]:.2f} x {figsize()[1]:.2f} in; "
        f"full page {GRID_COLS} x {GRID_ROWS} units, single column 2 x {GRID_ROWS} units ({HALF_W_MM:.0f} mm)",
        "fonts   " + ", ".join(f"{n} {v:.2f} -> {v * SCALE:.1f} pt" for n, v in _fonts.items()),
        "lines   " + ", ".join(f"{n} {v:.2f} -> {v * SCALE:.2f} pt" for n, v in _strokes.items()),
    ])


# ═══════════════════════════════════════════════════════════════════════
# save_panel — the one way a panel is written for a composed figure.
#
# Writes <letter>_<date>_<name>.svg and .png into the figure's folder at exactly the
# canvas the figure was created with. bbox_inches is never passed: "tight"
# crops the canvas to the drawn content, after which the panel is no longer a
# whole number of units and compose_figure has to rescale it.
#
# The file is always written. Before that, each problem found is printed as a
# "!!" line:
#   - a canvas more than SCALE_TOL away from a whole number of units
#   - anything drawn past the canvas edge, which is cut off in the file
#   - text or lines that print outside Nature's limits. The sizes are read
#     from the SVG being written, so text and lines drawn by seaborn,
#     adjustText and alphapepttools are checked as well.
# ═══════════════════════════════════════════════════════════════════════

def _svg_sizes(svg):
    """Font sizes and drawn stroke widths in a matplotlib SVG, in its own pt.

    Returns (fonts, strokes): a list of sizes and a list of (width, colour). A
    stroke counts only if it is drawn: colour not none, width and opacity above 0.
    """
    fonts, strokes = [], []
    for style in re.findall(r'style="([^"]*)"', svg):
        decl = {}
        for item in style.split(";"):
            if ":" in item:
                key, val = item.split(":", 1)
                decl[key.strip()] = val.strip()
        # matplotlib writes "font-size: 24.75px", older versions the "font:" shorthand
        m = re.search(r"([\d.]+)px", decl.get("font-size", decl.get("font", "")))
        if m:
            fonts.append(float(m.group(1)))
        colour = decl.get("stroke", "none")
        if colour == "none" or float(decl.get("stroke-opacity", 1)) == 0:
            continue
        width = float(decl.get("stroke-width", 1))      # SVG default width is 1
        if width > 0:
            strokes.append((width, colour))
    return fonts, strokes


def audit_panel(fig):
    """Every way `fig` leaves its canvas or prints outside the limits. Empty if none."""
    problems = []

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # Tick labels are left out: matplotlib keeps labels for ticks outside the axis
    # limits that are never drawn.
    ticks = {id(t) for ax in fig.axes
             for t in ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both")}
    texts = [t for t in fig.findobj(Text)
             if id(t) not in ticks and t.get_visible() and t.get_text().strip()]
    # get_tightbbox leaves out the width of x-axis labels and the height of y-axis
    # labels (matplotlib computes axis extents "for layout only"), so an axis label
    # longer than the canvas is not in it. The extents of all text are added.
    drawn = mpl.transforms.Bbox.union([fig.get_tightbbox(renderer).transformed(fig.dpi_scale_trans)]
                                      + [t.get_window_extent(renderer) for t in texts])   # pixels
    c, px_per_pt = fig.bbox, fig.dpi / 72
    over = {"left": c.x0 - drawn.x0, "bottom": c.y0 - drawn.y0,
            "right": drawn.x1 - c.x1, "top": drawn.y1 - c.y1}
    over = {side: d / px_per_pt for side, d in over.items() if d / px_per_pt > 1.0}  # 1 pt tolerance
    if over:
        names = sorted({t.get_text().replace("\n", " ")[:30] for t in texts
                        for b in [t.get_window_extent(renderer)]
                        if b.x0 < c.x0 - px_per_pt or b.y0 < c.y0 - px_per_pt
                        or b.x1 > c.x1 + px_per_pt or b.y1 > c.y1 + px_per_pt})
        problems.append("drawn past the canvas edge by "
                        + ", ".join(f"{d:.0f} pt {side}" for side, d in over.items())
                        + (f"; text outside: {names}" if names else ""))

    buf = io.StringIO()
    fig.savefig(buf, format="svg")
    fonts, strokes = _svg_sizes(buf.getvalue())

    lo_f, hi_f = TEXT_MIN_PT / SCALE, TEXT_MAX_PT / SCALE
    for size in sorted({f for f in fonts if not lo_f - 1e-6 <= f <= hi_f + 1e-6}):
        examples = sorted({t.get_text().replace("\n", " ")[:20] for t in fig.findobj(Text)
                           if t.get_visible() and t.get_text().strip()
                           and abs(t.get_fontsize() - size) < 0.05})[:4]
        problems.append(f"text {size:g} prints at {size * SCALE:.2f} pt "
                        f"(limits {TEXT_MIN_PT:g}-{TEXT_MAX_PT:g}): {examples}")

    lo_s, hi_s = STROKE_MIN_PT / SCALE, STROKE_MAX_PT / SCALE
    bad = {}
    for width, colour in strokes:
        if not lo_s - 1e-6 <= width <= hi_s + 1e-6:
            bad.setdefault(width, set()).add(colour)
    for width in sorted(bad):
        problems.append(f"line {width:g} prints at {width * SCALE:.2f} pt "
                        f"(limits {STROKE_MIN_PT:g}-{STROKE_MAX_PT:g}), colours {sorted(bad[width])}")
    return problems


def save_panel(fig, letter, name, folder, also_png=True):
    """Write `fig` as <letter>_<date>_<name>.svg and .png in `folder`; print size and problems."""
    w_in, h_in = fig.get_size_inches()
    w_mm, h_mm = w_in * MM_PER_IN * SCALE, h_in * MM_PER_IN * SCALE
    cols = (w_mm + GAP_MM) / (UNIT_W_MM + GAP_MM)
    rows = (h_mm + GAP_MM) / (UNIT_H_MM + GAP_MM)

    problems = []
    for n, side in ((cols, "wide"), (rows, "tall")):
        if abs(n / max(1, round(n)) - 1) > SCALE_TOL:
            problems.append(f"canvas is {n:.2f} units {side}, more than {SCALE_TOL:.0%} "
                            f"from a whole number; create the figure with figsize(cols, rows)")
    problems += audit_panel(fig)

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    stem = f"{letter}_{date.today():%Y%m%d}_{name}"
    at.pl.save_figure(fig, f"{stem}.svg", output_dir=folder)          # no bbox_inches
    if also_png:
        at.pl.save_figure(fig, f"{stem}.png", output_dir=folder, dpi=DPI)

    print(f"panel {letter}: {stem}.svg   {cols:.2f} x {rows:.2f} units, "
          f"prints {w_mm:.1f} x {h_mm:.1f} mm")
    for p in problems:
        print(f"  !! {p}")
    return folder / f"{stem}.svg"


def compose_figure(folder, layout, output="composed_figure.svg",
                   total_width_mm=PAGE_W_MM, max_height_mm=PAGE_H_MM, gap_mm=GAP_MM,
                   gap_x_mm=None, gap_y_mm=None,
                   panel_scale=None, panel_offset_mm=None,
                   panel_files=None, show_box=False,
                   font_range=PANEL_FONTS, stroke_range=PANEL_STROKES):
    """
    Compose a multi-panel figure from SVGs in a folder, within Nature's size limits.

    Panels are matched by the letter before the first underscore in the filename, and
    each one is scaled to the width of its COLUMN span. A row span never enlarges a
    panel; it only divides that panel's height over the rows it covers when row heights
    are computed.

    A panel must be an exact number of units wide, which is what save_panel guarantees.
    One that is not is rejected by name, because stretching it changes the size its text
    prints at and nothing downstream would show that. panel_scale overrides this per
    panel: a letter listed there is exempt from the check and is rendered at
    column_width * its multiplier. Use it for a panel that cannot be regenerated, such
    as an Illustrator figure or an older file saved with bbox_inches="tight"; the value
    that restores a uniform scale is native_units / colspan.

    gap_x_mm / gap_y_mm set the horizontal and the vertical gap between panels; each
    defaults to gap_mm. A panel that spans several rows already contains the gaps
    between them (figsize counts them), so they are subtracted before its height is
    shared out over the rows.

    panel_offset_mm = {letter: (dx, dy)} moves a panel's content right and down
    inside its cell, in mm; the panel letter stays at the corner of the cell.

    Each letter takes the NEWEST SVG (by modification time) whose name starts with
    "<letter>_"; the other candidates are listed. panel_files = {letter: filename}
    picks a file explicitly, e.g. one of two versions saved under the same letter.

    The saved SVG, PNG and TIFF are tight to the content and have a white background.
    show_box=True adds the page box to the DISPLAYED preview only.
    """
    folder = Path(folder)
    if panel_scale is None:
        panel_scale = {}
    if panel_offset_mm is None:
        panel_offset_mm = {}
    if panel_files is None:
        panel_files = {}
    gap_x_mm = gap_mm if gap_x_mm is None else gap_x_mm
    gap_y_mm = gap_mm if gap_y_mm is None else gap_y_mm

    rows = [row.split() for row in layout.strip().split("\n")]
    n_rows = len(rows)
    n_cols = len(rows[0])

    panels = {}
    for r, row in enumerate(rows):
        for c, letter in enumerate(row):
            if letter == ".":
                continue
            if letter not in panels:
                panels[letter] = {"r_min": r, "r_max": r, "c_min": c, "c_max": c}
            else:
                panels[letter]["r_min"] = min(panels[letter]["r_min"], r)
                panels[letter]["r_max"] = max(panels[letter]["r_max"], r)
                panels[letter]["c_min"] = min(panels[letter]["c_min"], c)
                panels[letter]["c_max"] = max(panels[letter]["c_max"], c)

    # The newest file per letter, by modification time.
    candidates = {}
    for f in folder.glob("*.svg"):
        panel_id = f.stem.split("_")[0]
        if panel_id in panels:
            candidates.setdefault(panel_id, []).append(f)
    svg_files = {}
    for letter, files in sorted(candidates.items()):
        if letter in panel_files:
            chosen = folder / panel_files[letter]
            if not chosen.exists():
                raise FileNotFoundError(f"panel_files[{letter!r}] = {panel_files[letter]!r} "
                                        f"is not in {folder}")
        else:
            files.sort(key=lambda p: p.stat().st_mtime, reverse=True)
            chosen = files[0]
        svg_files[letter] = chosen
        others = sorted(p.name for p in files if p != chosen)
        if others:
            print(f"{letter}: using {chosen.name}; also in the folder: {others}")

    missing = set(panels.keys()) - set(svg_files.keys())
    if missing:
        raise FileNotFoundError(f"No SVGs found for panels: {missing}. "
                                f"Available: {[f.name for f in folder.glob('*.svg')]}")

    unknown = set(panel_scale) - set(panels)
    if unknown:
        raise KeyError(f"panel_scale names panels that are not in the layout: {sorted(unknown)}")

    px_per_mm = 3.7795        # SVG px are 1/96 in: 96 / 25.4
    pt_per_px = 0.75   # panels are declared in pt, the composed canvas in px

    svg_data = {}
    for letter, fpath in svg_files.items():
        svg = st.fromfile(str(fpath))
        w, h = svg.get_size()
        if not w or not h:
            raise ValueError(
                f"{fpath.name} has no width/height on its root <svg> element, only a "
                f"viewBox. Illustrator writes this when 'Responsive' is checked on "
                f"export. Re-export with it unchecked, or add the attributes.")
        w = float(re.sub(r"[a-zA-Z%]+", "", w))
        h = float(re.sub(r"[a-zA-Z%]+", "", h))
        svg_data[letter] = {"w": w, "h": h}

    # ---- every panel is its column span in units, unless panel_scale overrides it ----
    bad = []
    for letter, span in panels.items():
        if letter in panel_scale:
            continue
        colspan = span["c_max"] - span["c_min"] + 1
        want_pt = figsize(colspan, 1)[0] * 72
        got_pt = svg_data[letter]["w"]
        if abs(got_pt - want_pt) / want_pt > 0.02:
            units = got_pt / want_pt * colspan
            bad.append(f"  {svg_files[letter].name}: {got_pt/72:.2f} in wide = {units:.2f} "
                       f"units, {colspan} columns need figsize({colspan}, rows) = "
                       f"{want_pt/72:.2f} in wide. Either create it with figsize({colspan}, "
                       f"rows) and save it with save_panel, or pass panel_scale={{'{letter}': "
                       f"{units/colspan:.3f}}} to keep it at the scale of the others.")
    if bad:
        raise ValueError("panels are not an exact number of units wide, so composing "
                         "them would change the size their text prints at:\n"
                         + "\n".join(bad))

    label_size = LABEL_PT / pt_per_px     # panel letters, 8 pt bold lowercase
    label_offset_x = -2
    label_offset_y = label_size - 2
    pad_top = label_size
    pad_left = label_size                 # taken OUT of the box, not added to it

    def _geometry(total_width_px, gap_x_px, gap_y_px):
        col_width = (total_width_px - pad_left - gap_x_px * (n_cols - 1)) / n_cols
        specs = {}
        for letter, span in panels.items():
            colspan = span["c_max"] - span["c_min"] + 1
            cell_w = colspan * col_width + (colspan - 1) * gap_x_px
            target_w = cell_w * panel_scale.get(letter, 1.0)
            scale = target_w / svg_data[letter]["w"]
            specs[letter] = {"scale": scale, "cell_w": cell_w, "target_w": target_w,
                             "cell_h": svg_data[letter]["h"] * scale, **span}

        # A panel spanning R rows also covers the R - 1 gaps between them, so only
        # the rest of its height is shared out over the rows.
        row_heights = [0] * n_rows
        for letter, spec in specs.items():
            rowspan = spec["r_max"] - spec["r_min"] + 1
            per_row = max(0.0, (spec["cell_h"] - (rowspan - 1) * gap_y_px) / rowspan)
            for r in range(spec["r_min"], spec["r_max"] + 1):
                row_heights[r] = max(row_heights[r], per_row)

        row_y = [0] * n_rows
        for r in range(1, n_rows):
            row_y[r] = row_y[r - 1] + row_heights[r - 1] + gap_y_px

        pos = {}
        content_x_max = content_y_max = 0
        for letter, spec in specs.items():
            x = spec["c_min"] * (col_width + gap_x_px) + pad_left
            y = row_y[spec["r_min"]] + pad_top
            pos[letter] = (x, y)
            dx, dy = (v * px_per_mm for v in panel_offset_mm.get(letter, (0.0, 0.0)))
            content_x_max = max(content_x_max, x + dx + spec["target_w"])
            content_y_max = max(content_y_max, y + dy + spec["cell_h"])

        return specs, pos, content_x_max, content_y_max + 5

    total_width_px = total_width_mm * px_per_mm
    gap_x_px, gap_y_px = gap_x_mm * px_per_mm, gap_y_mm * px_per_mm
    specs, pos, content_w, content_h = _geometry(total_width_px, gap_x_px, gap_y_px)

    binding = "width"
    if max_height_mm is not None:
        max_h_px = max_height_mm * px_per_mm
        for _ in range(10):
            if content_h <= max_h_px:
                break
            binding = "height"
            shrink = max_h_px / content_h
            total_width_px *= shrink
            gap_x_px *= shrink
            gap_y_px *= shrink
            specs, pos, content_w, content_h = _geometry(total_width_px, gap_x_px, gap_y_px)
        assert content_h <= max_h_px + 1e-6, "height cap not reached"

    def _build_panels():
        """Fresh Panel objects. Adding an element to a Figure re-parents it, so the
        saved figure and the preview each need their own set."""
        built = []
        for letter, spec in specs.items():
            x, y = pos[letter]
            # The offset moves the content only; the letter stays at the cell corner.
            dx, dy = (v * px_per_mm for v in panel_offset_mm.get(letter, (0.0, 0.0)))
            built.append(svgutils.Panel(
                svgutils.SVG(str(svg_files[letter])).scale(spec["scale"]).move(dx, dy),
                svgutils.Text(letter, x=label_offset_x, y=label_offset_y,
                              size=label_size, weight="bold", font="Arial"),
            ).move(x, y))
        return built

    canvas_w = max(content_w, total_width_px)
    svgutils.Figure(f"{canvas_w}px", f"{content_h}px", *_build_panels()).save(output)

    # 300 dpi RGB raster. The canvas is px at 96 dpi, so 300 dpi is a scale of 300/96.
    # background_color is what keeps the gaps between panels white: without it the
    # transparent pixels carry RGB (0, 0, 0) and Word shows them as black.
    stem = output.rsplit(".", 1)[0]
    cairosvg.svg2png(url=output, write_to=f"{stem}.png", scale=DPI / 96,
                     background_color="white")
    _im = Image.open(f"{stem}.png").convert("RGB")
    _im.save(f"{stem}.tif", compression="tiff_lzw", dpi=(DPI, DPI))
    _mb = Path(f"{stem}.tif").stat().st_size / 1e6

    if show_box and max_height_mm is not None:
        box_w_px, box_h_px = total_width_mm * px_per_mm, max_height_mm * px_per_mm
        # svgutils.Line sets stroke only, and an unset fill paints the closed
        # rectangle solid black, so fill has to be turned off explicitly.
        box = svgutils.Line([(0.5, 0.5), (box_w_px - 0.5, 0.5),
                             (box_w_px - 0.5, box_h_px - 0.5), (0.5, box_h_px - 0.5),
                             (0.5, 0.5)], width=1, color="#b0b0b0")
        box.root.set("fill", "none")
        _display(svgutils.Figure(f"{box_w_px}px", f"{box_h_px}px", *_build_panels(), box))
    else:
        _display(_SVG(output))

    print(f"Saved: {output}   (tight, {canvas_w/px_per_mm:.1f} x {content_h/px_per_mm:.1f} mm)")
    print(f"Saved: {stem}.tif   {_im.size[0]} x {_im.size[1]} px, {DPI} dpi, RGB, white, "
          f"{_mb:.1f} Mb" + ("" if _mb <= 10 else "   OVER the 10 Mb limit"))
    print(f"Panels: {', '.join(f'{k}={svg_files[k].name}' for k in sorted(svg_files))}")
    if panel_scale:
        print(f"panel_scale overrides in force: {panel_scale}")
    if panel_offset_mm:
        print(f"panel_offset_mm in force: {panel_offset_mm}")
    print(f"Gaps: {gap_x_mm:g} mm between columns, {gap_y_mm:g} mm between rows")

    box_w = total_width_mm
    box_h = max_height_mm if max_height_mm is not None else content_h / px_per_mm
    used_w, used_h = content_w / px_per_mm, content_h / px_per_mm
    print(f"Box      {box_w:6.1f} x {box_h:6.1f} mm")
    print(f"Content  {used_w:6.1f} x {used_h:6.1f} mm   "
          f"({used_w/box_w*100:.0f}% of width, {used_h/box_h*100:.0f}% of height)")
    print(f"Unused   {box_w-used_w:6.1f} x {box_h-used_h:6.1f} mm   "
          f"limited by {binding}; grid is {n_cols} x {n_rows}")
    print(f"Printed text ({font_range[0]:g}-{font_range[1]:g} pt native) and strokes "
          f"({stroke_range[0]:g}-{stroke_range[1]:g} pt native). "
          f"Nature: text {TEXT_MIN_PT:g}-{TEXT_MAX_PT:g} pt, "
          f"strokes {STROKE_MIN_PT:g}-{STROKE_MAX_PT:g} pt.")
    for letter in sorted(specs):
        s = specs[letter]["scale"] * pt_per_px
        t_lo, t_hi = font_range[0] * s, font_range[1] * s
        k_lo, k_hi = stroke_range[0] * s, stroke_range[1] * s
        why = []
        if t_lo < TEXT_MIN_PT:     why.append(f"text {t_lo:.1f} < {TEXT_MIN_PT:g} pt")
        if t_hi > TEXT_MAX_PT:     why.append(f"text {t_hi:.1f} > {TEXT_MAX_PT:g} pt")
        if k_lo < STROKE_MIN_PT:   why.append(f"stroke {k_lo:.2f} < {STROKE_MIN_PT:g} pt")
        if k_hi > STROKE_MAX_PT:   why.append(f"stroke {k_hi:.2f} > {STROKE_MAX_PT:g} pt")
        mark = " *" if letter in panel_scale else "  "
        print(f"  {letter}{mark}{specs[letter]['target_w']/px_per_mm:6.1f} x "
              f"{specs[letter]['cell_h']/px_per_mm:6.1f} mm   scale {s:.3f}   "
              f"text {t_lo:4.1f}-{t_hi:4.1f} pt   strokes {k_lo:4.2f}-{k_hi:4.2f} pt"
              + ("" if not why else "   <-- " + "; ".join(why)))
