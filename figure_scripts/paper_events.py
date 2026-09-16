#!/usr/bin/env python3
"""
Canonical list of the events marked on the paper's time-series panels.

Defined once and imported, so that Fig. 2b and Fig. 3d cannot drift apart. They did
drift before: Fig. 2b and Fig. 4a carry different sets, with the Valanginian event at
134 Ma in one and 133 Ma in the other, and an early Aptian event at 115 Ma in one
against 120 Ma in the other. Fig. 4a is left as it stands for now, since its panel
also names ages in the labels; if it is ever brought into line, it should import from
here rather than keep a second list.
"""
from __future__ import annotations

# (age in Ma, label). Ages on GTS2020, matching the rest of the workflow.
EVENTS = [
    (134.0, "WE"),      # Weissert
    (115.0, "LACI"),    # late Aptian carbonate crisis / OAE1b interval
    (93.9, "OAE2"),
    (56.0, "PETM"),
    (33.9, "EOT"),
    (15.2, "MMCO"),
]

# Drawing style, shared so the two panels look identical.
LINE_KW = dict(color="0.45", ls=":", lw=1.0, zorder=2)
TEXT_KW = dict(color="0.30", ha="center", va="bottom", zorder=6)


def draw_events(ax, fontsize: float = 8.5, y: float = 1.005, pad_pt: float = 2.0,
                lw_scale: float = 1.0) -> None:
    """Dotted verticals with labels just above the axes.

    A label that would touch its neighbour is lifted onto a second line instead. On a
    narrow panel EOT and MMCO are 18.7 Myr apart, which is about 6 mm, while the two
    strings need nearly 12 - they used to run together into EOTMMCO, and a check that
    only fires on real overlap called that clear.
    """
    fig = ax.figure
    # lw_scale lets the calling figure put these verticals on its own line-weight
    # scale; left at 1.0 they are drawn exactly as LINE_KW specifies.
    line_kw = dict(LINE_KW, lw=LINE_KW["lw"] * lw_scale)
    drawn = []
    for age, label in EVENTS:
        ax.axvline(age, **line_kw)
        drawn.append((age, ax.text(age, y, label, transform=ax.get_xaxis_transform(),
                                   fontsize=fontsize, **TEXT_KW)))
    fig.canvas.draw()
    rr = fig.canvas.get_renderer()
    ax_h_in = ax.get_window_extent(renderer=rr).height / fig.dpi
    # 1.6 line heights, not 1.0: a lift of exactly one line leaves the two rows
    # separated by a fraction of a point, which reads as one crowded block.
    line = fontsize * 1.60 / 72.0 / ax_h_in        # in axes fractions
    pad_px = pad_pt / 72.0 * fig.dpi
    rows = [None, None]                            # right edge of the last label per row
    for age, t in sorted(drawn, key=lambda p: -p[0]):
        e = t.get_window_extent(renderer=rr)
        row = 0 if rows[0] is None or e.x0 - pad_px > rows[0] else 1
        if row:
            t.set_y(y + line)
            fig.canvas.draw()
            e = t.get_window_extent(renderer=rr)
        rows[row] = e.x1


def draw_events_gmt(fig, panel, pt: float = 7.0, colour: str = "40",
                    label_colour: str = "30", dy: float = 0.16):
    """The same events on a pyGMT panel: dotted verticals inside the frame and
    labels just above it. Returns the label boxes so the caller's collision
    check can see them.

    `panel` is a paper_gmt.Panel; `dy` lifts the labels off the frame, in cm.
    """
    import paper_gmt as S

    region = list(panel.region)
    boxes = []
    for age, label in EVENTS:
        if not (min(region[0], region[1]) <= age <= max(region[0], region[1])):
            continue
        fig.plot(x=[age, age], y=[region[2], region[3]], pen=f"0.5p,{colour},.",
                 region=region, projection=panel.projection)
        boxes.append(panel.label_box(f"event {label}", label, age, region[3], pt,
                                     justify="CB", dy=dy))
    cf = S.cm_frame(panel, top=max(b.y1 for b in boxes) - panel.height if boxes else 0.0)
    for (age, label), box in zip([e for e in EVENTS
                                  if min(region[0], region[1]) <= e[0] <= max(region[0], region[1])],
                                 boxes):
        fig.text(x=panel.x_cm(age), y=panel.height + dy, text=label, justify="CB",
                 font=f"{pt}p,{S.FONT},{label_colour}", no_clip=True, **cf)
    return boxes
