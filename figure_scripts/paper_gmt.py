#!/usr/bin/env python3
"""Shared pyGMT style for the manuscript figures.

Everything the four figure scripts need to look like one set: the typeface and
type sizes, the printed widths Geology uses, a data-to-centimetre mapper, and a
text-collision check built on Helvetica's own metrics.

Geology is a three-column journal. Its artwork guidance asks for Helvetica or
Arial, body lettering of 7-12 pt and part labels (A, B, C) of 13-16 pt, all at
the size the figure is printed. pyGMT lays out in real centimetres, so a figure
built here is already at printed size and the point sizes below are what the
reader sees - no reduction factor to reason about.

Collision detection works on Helvetica's character widths (the Adobe metrics,
reproduced below), which are the same numbers GMT uses to set the text, so a
label's box can be computed without rendering. Use `Panel.label_box()` to get a
box for every string a script places by hand, then pass the list to
`collisions()`; anything the frame draws (axis annotations and labels) GMT keeps
clear by construction.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# ---- paths -----------------------------------------------------------------
# Works in the working tree (this file in Paper/, the workflow in a sibling
# CCD_workflow_clean/) and in the public repo (this file in figure_scripts/,
# the workflow at the repo root).
_HERE = Path(__file__).resolve().parent


def workflow_root() -> Path:
    for p in (_HERE.parent / "CCD_workflow_clean", _HERE.parent, _HERE):
        if (p / "steps").is_dir() and (p / "ccdworkflow").is_dir():
            return p
    raise SystemExit("cannot locate the CCD workflow root (expected steps/ + ccdworkflow/)")


CW = workflow_root()
FIGDIR = _HERE / "Figures" if (_HERE / "Figures").is_dir() else CW / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

# ---- journal geometry and type ---------------------------------------------
W_1COL = 5.9      # cm, Geology single column
W_2COL = 12.28    # cm, two columns
W_PAGE = 18.5     # cm, full page

FONT = "Helvetica"
PT_ANNOT = 7.0    # axis annotations; the journal floor
PT_LABEL = 8.0    # axis labels
PT_LEG = 7.0      # legend and in-panel annotation
PT_TAG = 13.0     # part labels A, B, C, D; the journal floor for these
PT_MIN = 7.0      # nothing on the page may be smaller

# GMT clips anything drawn left of the page origin, and axis annotations, axis labels
# and part letters all live left of their panel, so every script starts by shifting the
# origin right by this much. Wide enough for a five-digit annotation, a rotated label
# and the part letter.
MARGIN_L = 1.75   # cm
MARGIN_B = 1.10   # cm: the bottom annotation row and axis label are below the
                  # origin too, and GMT clips them at the page edge just the same

PEN_FRAME = "0.5p,black"
PEN_TICK = "0.5p,black"

CM_PER_PT = 2.54 / 72.0


def defaults() -> dict:
    """GMT settings shared by every figure. Pass to `pygmt.config(**defaults())`."""
    return dict(
        FONT_ANNOT_PRIMARY=f"{PT_ANNOT}p,{FONT},black",
        FONT_ANNOT_SECONDARY=f"{PT_ANNOT}p,{FONT},black",
        FONT_LABEL=f"{PT_LABEL}p,{FONT},black",
        FONT_TAG=f"{PT_TAG}p,{FONT}-Bold,black",
        FONT_TITLE=f"{PT_LABEL}p,{FONT},black",
        MAP_FRAME_TYPE="plain",
        MAP_FRAME_PEN=PEN_FRAME,
        MAP_TICK_PEN_PRIMARY=PEN_TICK,
        MAP_TICK_LENGTH_PRIMARY="0.10c/0.05c",
        MAP_ANNOT_OFFSET_PRIMARY="0.08c",
        MAP_LABEL_OFFSET="0.10c",
        MAP_TITLE_OFFSET="0.10c",
        MAP_VECTOR_SHAPE="0.5",
        PS_CHAR_ENCODING="ISOLatin1+",
        PROJ_LENGTH_UNIT="c",
        PS_LINE_JOIN="round",
        PS_LINE_CAP="round",
        FORMAT_GEO_MAP="dddF",
    )


# ---- Helvetica metrics -----------------------------------------------------
# Character widths in 1/1000 em for codes 32-126, from the Adobe AFM files that
# ship with GMT's PostScript library. Anything outside that range falls back to
# the mean width, which is close enough for a clearance test.
_W_REG = [278, 278, 355, 556, 556, 889, 667, 222, 333, 333, 389, 584, 278, 333, 278, 278,
          556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 278, 278, 584, 584, 584, 556,
          1015, 667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722, 778,
          667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 278, 278, 278, 469, 556,
          222, 556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, 556, 556,
          556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500, 334, 260, 334, 584]
_W_BOLD = [278, 333, 474, 556, 556, 889, 722, 278, 333, 333, 389, 584, 278, 333, 278, 278,
           556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 333, 333, 584, 584, 584, 611,
           975, 722, 722, 722, 722, 667, 611, 778, 722, 278, 556, 722, 611, 833, 722, 778,
           667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 333, 278, 333, 584, 556,
           278, 556, 611, 556, 611, 556, 333, 611, 611, 278, 278, 556, 278, 889, 611, 611,
           611, 611, 389, 556, 333, 611, 556, 778, 556, 556, 500, 389, 280, 389, 584]
_ASCENT, _DESCENT = 0.718, 0.207   # em, Helvetica


def text_width_cm(text: str, size_pt: float, bold: bool = False) -> float:
    w = _W_BOLD if bold else _W_REG
    mean = sum(w) / len(w)
    total = sum(w[ord(c) - 32] if 32 <= ord(c) <= 126 else mean for c in text)
    return total / 1000.0 * size_pt * CM_PER_PT


def text_height_cm(size_pt: float, lines: int = 1, leading: float = 1.2) -> float:
    one = (_ASCENT + _DESCENT) * size_pt * CM_PER_PT
    return one if lines <= 1 else one + (lines - 1) * leading * size_pt * CM_PER_PT


@dataclass
class Box:
    """A rectangle on the page, in centimetres from the current origin."""
    name: str
    x0: float
    y0: float
    x1: float
    y1: float

    def overlaps(self, other: "Box", pad: float = 0.0) -> bool:
        return (self.x0 - pad < other.x1 and other.x0 - pad < self.x1
                and self.y0 - pad < other.y1 and other.y0 - pad < self.y1)


def collisions(boxes, pad_cm: float = 0.03) -> list[str]:
    """Pairs of boxes that touch, as 'a / b' strings. Empty means clear."""
    bad = []
    for i, a in enumerate(boxes):
        for b in boxes[i + 1:]:
            if a.overlaps(b, pad_cm):
                bad.append(f"{a.name} / {b.name}")
    return bad


def report(bad, what: str = "text") -> None:
    print(f"  {what} collisions: " + (", ".join(bad) if bad else "none"))


def too_small(*sizes_pt) -> list[float]:
    """Any type size below the journal floor."""
    return [s for s in sizes_pt if s < PT_MIN - 1e-9]


# ---- linear panels ---------------------------------------------------------
@dataclass
class Panel:
    """A linear pyGMT panel: its region, its size in cm, and where it sits.

    `x_reversed` covers the age axes, which run from old on the left to 0 on the
    right - GMT does that with a negative projection width, so the mapping from
    data to centimetres has to know about it.
    """
    region: tuple
    width: float
    height: float
    origin: tuple = (0.0, 0.0)     # lower-left corner of the frame, in cm
    x_reversed: bool = False
    y_reversed: bool = False

    @property
    def projection(self) -> str:
        w = -self.width if self.x_reversed else self.width
        h = -self.height if self.y_reversed else self.height
        return f"X{w}c/{h}c"

    def x_cm(self, x) -> float:
        x0, x1 = self.region[0], self.region[1]
        f = (float(x) - x0) / (x1 - x0)
        if self.x_reversed:
            f = 1.0 - f
        return self.origin[0] + f * self.width

    def y_cm(self, y) -> float:
        y0, y1 = self.region[2], self.region[3]
        f = (float(y) - y0) / (y1 - y0)
        if self.y_reversed:
            f = 1.0 - f
        return self.origin[1] + f * self.height

    def frame_box(self, name: str = "frame") -> Box:
        return Box(name, self.origin[0], self.origin[1],
                   self.origin[0] + self.width, self.origin[1] + self.height)

    def label_box(self, name: str, text: str, x, y, size_pt: float,
                  justify: str = "LB", bold: bool = False,
                  dx: float = 0.0, dy: float = 0.0, lines: int = 1) -> Box:
        """Box of a string placed at data coordinates (x, y) with GMT justify.

        dx/dy are the offsets given to pygmt.text(offset=...), in cm.
        """
        w = max(text_width_cm(t, size_pt, bold) for t in text.split("\n"))
        h = text_height_cm(size_pt, lines=max(lines, len(text.split("\n"))))
        cx, cy = self.x_cm(x) + dx, self.y_cm(y) + dy
        hx, vy = justify[0].upper(), justify[1].upper()
        x0 = cx if hx == "L" else (cx - w if hx == "R" else cx - w / 2)
        y0 = cy if vy == "B" else (cy - h if vy == "T" else cy - h / 2)
        return Box(name, x0, y0, x0 + w, y0 + h)

    def data_in_box(self, box: Box, xs, ys) -> bool:
        """True if any data point falls inside the box - used to keep legends
        and in-panel annotation off the curves they describe."""
        for x, y in zip(xs, ys):
            try:
                cx, cy = self.x_cm(x), self.y_cm(y)
            except (TypeError, ValueError):
                continue
            if cx != cx or cy != cy:     # NaN
                continue
            if box.x0 <= cx <= box.x1 and box.y0 <= cy <= box.y1:
                return True
        return False


def cm_frame(panel: Panel, top: float = 0.0, right: float = 0.0) -> dict:
    """region/projection that make the panel's own centimetres the coordinates,
    so anything measured by the collision check can be drawn where it was
    measured.

    `top` and `right` extend the drawable area beyond the frame, anchored at the
    same lower-left corner. GMT clips lines and polygons at the region boundary
    whatever -N says - only text escapes it - so anything drawn above or beside
    a panel has to be inside an enlarged region rather than outside a small one.
    """
    w, h = panel.width + right, panel.height + top
    return dict(region=[0, w, 0, h], projection=f"X{w}c/{h}c")


def draw_legend(fig, panel: Panel, entries, x0: float, y0: float, ncol: int = 1,
                pt: float = PT_LEG, sym_w: float = 0.32, gap: float = 0.10,
                row_lead: float = 1.30, col_gap: float = 0.25,
                frame: bool = False, pad: float = 0.06) -> Box:
    """Draw a legend in the panel's centimetre coordinates, from the top-left
    corner (x0, y0), and return the box it occupies.

    `entries` are (text, colour, kind, pen) with kind 'line' or 'patch'. Built
    by hand rather than with GMT's own legend so the figure scripts can measure
    the box first and prove it clears the data.
    """
    rows = -(-len(entries) // ncol)
    row = text_height_cm(pt) * row_lead
    col_w = []
    for c in range(ncol):
        part = entries[c * rows:(c + 1) * rows]
        col_w.append(sym_w + gap + max(text_width_cm(t, pt) for t, *_ in part) if part else 0.0)
    total_w = sum(col_w) + col_gap * (ncol - 1)
    total_h = row * rows
    box = Box("legend", x0 - pad, y0 - total_h - pad, x0 + total_w + pad, y0 + pad)
    cf = cm_frame(panel, top=max(0.0, box.y1 - panel.height),
                  right=max(0.0, box.x1 - panel.width))
    if frame:
        fig.plot(x=[box.x0, box.x1, box.x1, box.x0], y=[box.y0, box.y0, box.y1, box.y1],
                 pen="0.25p,gray50", fill="white@15", close=True, no_clip=True, **cf)
    for i, (text, colour, kind, pen) in enumerate(entries):
        c, r = divmod(i, rows)
        cx = x0 + sum(col_w[:c]) + col_gap * c
        cy = y0 - row * (r + 0.5)
        if kind == "line":
            # A GMT pen is width,colour,style, so an entry that carries a dash pattern
            # has to have the colour spliced into the middle rather than appended.
            _w, *_st = str(pen).split(",")
            _pen = ",".join([_w, colour] + _st)
            fig.plot(x=[cx, cx + sym_w], y=[cy, cy], pen=_pen, no_clip=True, **cf)
        else:
            h = 0.045
            fig.plot(x=[cx, cx + sym_w, cx + sym_w, cx], y=[cy - h, cy - h, cy + h, cy + h],
                     fill=colour, pen=pen or "0.2p,gray60", close=True, no_clip=True, **cf)
        fig.text(x=cx + sym_w + gap, y=cy, text=text, justify="LM",
                 font=f"{pt}p,{FONT},black", no_clip=True, **cf)
    return box


def rotated_label(fig, panel: Panel, lines, side: str = "left", pt: float = PT_LABEL,
                  colour: str = "black", gap: float = 0.10) -> Box:
    """A two-line axis label, rotated, beside a panel.

    GMT axis labels are one line, and these panels carry a quantity and its units that
    together run longer than the panel is tall. Drawn by hand, line by line, and
    measured so the caller knows how much room the label takes.
    """
    cf = cm_frame(panel)
    line_h = text_height_cm(pt)
    width = len(lines) * line_h
    x_in = panel.origin[0] - gap if side == "left" else panel.origin[0] + panel.width + gap
    ys = panel.origin[1] + panel.height / 2
    for i, text in enumerate(lines):
        off = (i + 0.5) * line_h
        x = x_in - off if side == "left" else x_in + off
        fig.text(x=x, y=ys, text=text, angle=90 if side == "left" else 270,
                 justify="CM", font=f"{pt}p,{FONT},{colour}", no_clip=True, **cf)
    x0 = x_in - width if side == "left" else x_in
    return Box(f"{side} label", x0, panel.origin[1], x0 + width, panel.origin[1] + panel.height)


# ---- page and colour -------------------------------------------------------
def begin(fig, left: float = None, bottom: float = None) -> None:
    """Reserve the left and bottom margins. Call once, before anything is drawn."""
    fig.shift_origin(xshift=f"{MARGIN_L if left is None else left}c",
                     yshift=f"{MARGIN_B if bottom is None else bottom}c")


def rgb(colour: str) -> str:
    """'#1f77b4' -> '31/119/180'. GMT's pen parser mistakes a hex colour for a dash
    pattern when a style follows it, so pens carry colours in r/g/b form."""
    c = colour.strip()
    if not c.startswith("#"):
        return c
    c = c[1:]
    if len(c) == 3:
        c = "".join(ch * 2 for ch in c)
    return "/".join(str(int(c[i:i + 2], 16)) for i in (0, 2, 4))



def _int_rgb(c: str) -> str:
    """GMT rejects fractional r/g/b and falls back to red, which is how a foreground
    colour of '134.5/45/6.5' painted every deep-ocean cell bright red."""
    if "/" not in c:
        return c
    try:
        return "/".join(str(int(round(float(v)))) for v in c.split("/"))
    except ValueError:
        return c


def pale_rgb(colour: str, chroma: float = 0.45, l_floor: float = 55.0,
             l_scale: float = 0.42) -> str:
    """Wash a colour out the way the project's pale thickness scale was made: in CIE
    LCh, chroma is cut to `chroma` of its value and lightness remapped
    L -> l_floor + l_scale * L. Applying it to a colour map keeps the hue order and the
    class boundaries while lifting everything towards the paper.
    """
    import math

    r, g, b = (int(v) for v in _int_rgb(colour).split("/"))

    def lin(c):
        c = c / 255.0
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4

    R, G, B = lin(r), lin(g), lin(b)
    X = (0.4124 * R + 0.3576 * G + 0.1805 * B) / 0.95047
    Y = (0.2126 * R + 0.7152 * G + 0.0722 * B) / 1.00000
    Z = (0.0193 * R + 0.1192 * G + 0.9505 * B) / 1.08883

    def f(t):
        return t ** (1 / 3) if t > 216 / 24389 else (841 / 108) * t + 4 / 29

    fx, fy, fz = f(X), f(Y), f(Z)
    L, a, bb = 116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz)
    C, h = math.hypot(a, bb), math.atan2(bb, a)
    L, C = l_floor + l_scale * L, C * chroma
    a, bb = C * math.cos(h), C * math.sin(h)

    fy = (L + 16) / 116
    fx, fz = fy + a / 500, fy - bb / 200

    def finv(t):
        return t ** 3 if t ** 3 > 216 / 24389 else (t - 4 / 29) * (108 / 841)

    X, Y, Z = finv(fx) * 0.95047, finv(fy), finv(fz) * 1.08883
    R = 3.2406 * X - 1.5372 * Y - 0.4986 * Z
    G = -0.9689 * X + 1.8758 * Y + 0.0415 * Z
    B = 0.0557 * X - 0.2040 * Y + 1.0570 * Z

    def enc(c):
        c = max(0.0, min(1.0, c))
        c = 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055
        return str(int(round(c * 255)))

    return "/".join(enc(c) for c in (R, G, B))


def stepped_cpt(master: str, bounds, out: Path, reverse: bool = False,
                background: str | None = None, foreground: str | None = None,
                nan: str = "189/189/189", header: str = "", pale: float | None = None) -> Path:
    """Write a CPT whose colours step evenly through `master` while its class
    boundaries stay at `bounds`.

    Sampling a colour map across uneven classes by z compresses every colour
    change into the wide classes at the top; stepping by class index instead
    gives each class the same share of the colour range, which is what keeps the
    thin low-thickness bands apart. Returns the path written.
    """
    import pygmt

    n = len(bounds) - 1
    tmp = out.with_suffix(".master.cpt")
    pygmt.makecpt(cmap=master, series=[0, n, 1], reverse=reverse, output=str(tmp))
    rows = []
    for line in tmp.read_text().splitlines():
        s = line.split()
        if not s or s[0].startswith("#") or s[0] in ("B", "F", "N"):
            continue
        rows.append(s[1])
    if len(rows) < n:
        raise SystemExit(f"{master}: got {len(rows)} colours for {n} classes")
    rows = [_int_rgb(c) for c in rows[:n]]
    if pale is not None:
        rows = [pale_rgb(c, chroma=pale) for c in rows]
    lines = [f"# {header}" if header else "# stepped CPT",
             f"# {master}{' reversed' if reverse else ''}, one colour class per interval",
             "# COLOR_MODEL = RGB"]
    for i in range(n):
        lines.append(f"{bounds[i]}\t{rows[i]}\t{bounds[i + 1]}\t{rows[i]}")
    lines.append(f"B\t{background or rows[0]}")
    lines.append(f"F\t{foreground or rows[-1]}")
    lines.append(f"N\t{nan}")
    out.write_text("\n".join(lines) + "\n")
    tmp.unlink(missing_ok=True)
    return out
