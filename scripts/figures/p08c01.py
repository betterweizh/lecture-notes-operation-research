"""Figures of chapter 8.1, LP Modeling and Geometry (Part VIII).

* ``08-01-capital-activities``: the two two-period production activities of
  the capital/consumer goods example (hand-written SVG, words translated);
* ``08-01-oilco-region``: graphical solution of the minimum-capacity
  refinery LP (unbounded region, optimum (55, 30));
* ``08-01-basic-solutions-2d``: max 2x1 + 3x2 with its four basic feasible
  and two infeasible basic solutions;
* ``08-01-feasible-segment``: projection of the feasible set of
  2x1 + 3x2 + x3 = 5, x1 + 2x2 - x3 = 3, x >= 0 on the (x1, x2) plane;
* ``08-01-lad-fit``: typing-speed data with the least-absolute-deviation line
  and the least-squares line.
"""

from fractions import Fraction as F

import matplotlib.pyplot as plt
import numpy as np

from .common import (
    BLUE,
    INK,
    MUTED,
    ORANGE,
    draw_line,
    lp_region,
    save_mpl,
    svg_line,
    svg_text,
    write_svg,
)


def _halo(fig):
    """White outline behind every text in the axes so that labels stay legible
    where they cross constraint lines (inverted with the figure in dark mode).
    """
    import matplotlib.patheffects as pe

    for ax in fig.axes:
        for t in ax.texts:
            t.set_path_effects(
                [pe.withStroke(linewidth=3, foreground="white")]
            )


# Words of the capital-activities diagram, per edition.
_WORDS = {
    "en": {
        "cap": "capital",
        "con": "consumer goods",
        "t": "period t",
        "t2": "period t + 2",
    },
    "zh": {
        "cap": "资本品",
        "con": "消费品",
        "t": "第 t 期",
        "t2": "第 t + 2 期",
    },
}


def _dot(ax, x, y, *, filled=True, color=BLUE, size=6.5, z=6):
    """Mark a point: filled dot (feasible) or hollow dot (infeasible).

    Not clipped: points on the axes, e.g. (0, 0) in the corner, stay whole.
    """
    ax.plot(
        [float(x)], [float(y)], marker="o", markersize=size, zorder=z,
        color=color, markerfacecolor=color if filled else "white", markeredgewidth=1.6,
        clip_on=False,
    )  # fmt: skip


def capital_activities(outdir, lang="en"):
    """Input-output diagrams of the activities x_t1 and x_t2 (ex-08-01-capital-consumer).

    Left: 1 unit of capital plus 3 consumer goods in period t give 4 units
    of capital in period t + 2. Right: 1 unit of capital in period t gives
    1 unit of capital and 10 consumer goods in period t + 2.
    """
    w = _WORDS.get(lang, _WORDS["en"])
    body = ""

    def box(x):
        # The process box (2 periods), 30 px wide, 100 px tall.
        return (f'<rect x="{x}" y="45" width="30" height="100" rx="3" '
                f'fill="none" stroke="{INK}" stroke-width="1.6"/>')  # fmt: skip

    def arrow(x1, y, x2, label, color=INK):
        # A horizontal flow arrow with its label above the middle.
        s = svg_line(x1, y, x2, y, arrow=True, color=color, width=1.8)
        return s + svg_text(
            (x1 + x2) / 2, y - 13, label, size=13, math=False, color=color
        )

    # Panel (a): activity x_t1, capital makes capital. Box at x = 160..190.
    body += box(160)
    body += svg_text(175, 28, "x_{t1}", size=14)
    body += arrow(10, 70, 160, f"1 {w['cap']}", BLUE)
    body += arrow(10, 120, 160, f"3 {w['con']}", ORANGE)
    body += arrow(190, 95, 300, f"4 {w['cap']}", BLUE)
    body += svg_text(85, 168, w["t"], size=13, math=False, color=MUTED)
    body += svg_text(245, 168, w["t2"], size=13, math=False, color=MUTED)

    # Panel (b): activity x_t2, capital makes consumer goods. Box at ox + 120.
    ox = 330
    body += box(ox + 120)
    body += svg_text(ox + 135, 28, "x_{t2}", size=14)
    body += arrow(ox + 10, 95, ox + 120, f"1 {w['cap']}", BLUE)
    body += arrow(ox + 150, 70, ox + 300, f"1 {w['cap']}", BLUE)
    body += arrow(ox + 150, 120, ox + 300, f"10 {w['con']}", ORANGE)
    body += svg_text(ox + 65, 168, w["t"], size=13, math=False, color=MUTED)
    body += svg_text(ox + 225, 168, w["t2"], size=13, math=False, color=MUTED)
    return write_svg("08-01-capital-activities", body, 640, 185, outdir)


def oilco_region(outdir, lang="en"):
    """Graphical solution of the minimum-capacity refinery LP (ex-08-01-oilco).

    min x1 + x2 (thousand bbl/day) s.t. x2 <= 1.5 x1 (Iraqi share),
    diesel 0.2x1 + 0.1x2 >= 14, gasoline 0.25x1 + 0.6x2 >= 30,
    lubricants 0.1x1 + 0.15x2 >= 10, jet 0.15x1 + 0.1x2 >= 8.
    Coefficients are scaled to integers for exact vertices; the labels show
    the original rows.
    """
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    xlim, ylim = (0, 130), (0, 100)
    cons = [
        {"a": (-3, 2), "op": "<=", "b": 0, "label": "$x_2 = 1.5x_1$", "label_at": 0.83},
        {"a": (2, 1), "op": ">=", "b": 140, "label": "$0.2x_1 + 0.1x_2 = 14$", "label_at": 0.2},
        {"a": (5, 12), "op": ">=", "b": 600, "label": "$0.25x_1 + 0.6x_2 = 30$", "label_at": 0.76},
        {"a": (2, 3), "op": ">=", "b": 200, "label": "$0.1x_1 + 0.15x_2 = 10$", "label_at": 0.13},
        {"a": (3, 2), "op": ">=", "b": 160, "label": "$0.15x_1 + 0.1x_2 = 8$", "label_at": 0.66,
         "linestyle": ":"},
    ]  # fmt: skip
    vertices = lp_region(ax, cons, xlim, ylim, label_size=9, exact=True)
    # Corners of the unbounded region: (40,60), (55,30), (200/3,200/9), (120,0).
    for v in vertices:
        _dot(ax, *v, color=INK, size=4.5)
    # Optimal contour x1 + x2 = 85 and the optimum.
    draw_line(
        ax, 1, 1, 85, xlim, ylim, color=INK, linestyle="--", linewidth=1.3
    )
    ax.text(3, 85, "$z = 85$", color=INK, fontsize=10, va="bottom")
    _dot(ax, 55, 30, color=ORANGE)
    ax.annotate(r"$\mathbf{x}^* = (55,\ 30)$", xy=(55, 30), xytext=(10, 6),
                textcoords="offset points", fontsize=10)  # fmt: skip
    _halo(fig)
    return save_mpl(fig, "08-01-oilco-region", outdir)


def basic_solutions_2d(outdir, lang="en"):
    """max 2x1 + 3x2, x1 + 3x2 <= 12, 3x1 + 2x2 <= 12 (ex-08-01-basic-solutions-2d).

    Filled dots: the four basic feasible solutions; hollow dots: the
    infeasible basic solutions (12, 0) and (0, 6); dashed: the optimal
    contour 2x1 + 3x2 = 96/7 through (12/7, 24/7).
    """
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    xlim, ylim = (0, 13), (0, 7)
    lp_region(ax, [{"a": (1, 3), "op": "<=", "b": 12, "label_at": 0.75},
                   {"a": (3, 2), "op": "<=", "b": 12, "label_at": 0.75}],
              xlim, ylim)  # fmt: skip
    draw_line(ax, 2, 3, F(96, 7), xlim, ylim, color=INK, linestyle="--", linewidth=1.3,
              label="$z = 96/7$", label_at=0.72, label_color=INK)  # fmt: skip
    # The four basic feasible solutions (corners); only the optimum is labeled.
    for x, y in [(0, 0), (4, 0), (0, 4), (F(12, 7), F(24, 7))]:
        _dot(ax, x, y)
    ax.annotate(r"$(\frac{12}{7},\ \frac{24}{7})$", xy=(12 / 7, 24 / 7), xytext=(8, 6),
                textcoords="offset points")  # fmt: skip
    for (x, y), lab, off in [
        ((12, 0), "$(12,\\ 0)$", (-20, 8)),
        ((0, 6), "$(0,\\ 6)$", (8, 0)),
    ]:
        _dot(ax, x, y, filled=False, color=ORANGE)
        ax.annotate(
            lab,
            xy=(x, y),
            xytext=off,
            textcoords="offset points",
            color=ORANGE,
        )
    _halo(fig)
    return save_mpl(fig, "08-01-basic-solutions-2d", outdir)


def feasible_segment(outdir, lang="en"):
    """Feasible set of 2x1 + 3x2 + x3 = 5, x1 + 2x2 - x3 = 3, x >= 0 (ex-08-01-feasible-segment).

    Projected on (x1, x2) it is the segment of 3x1 + 5x2 = 8 between the
    basic feasible solutions (1, 1, 0) and (0, 8/5, 1/5); beyond (1, 1) the
    line has x3 < 0 and ends at the infeasible basic solution (8/3, 0, -1/3).
    """
    fig, ax = plt.subplots(figsize=(4.6, 3.3))
    xlim, ylim = (0, 3.5), (0, 2)
    # The whole line 3x1 + 5x2 = 8, thin and muted.
    draw_line(ax, 3, 5, 8, xlim, ylim, color=MUTED, linewidth=1.2)
    ax.text(0.05, 0.78, "$3x_1 + 5x_2 = 8$", color=MUTED)  # below the segment
    # The feasible segment, thick.
    ax.plot(
        [0, 1],
        [1.6, 1],
        color=BLUE,
        linewidth=4,
        solid_capstyle="round",
        zorder=4,
    )
    _dot(ax, 1, 1)
    _dot(ax, 0, F(8, 5))
    _dot(ax, F(8, 3), 0, filled=False, color=ORANGE)
    ax.annotate(
        "$(1,\\ 1,\\ 0)$", xy=(1, 1), xytext=(8, 6), textcoords="offset points"
    )
    ax.annotate(r"$(0,\ \frac{8}{5},\ \frac{1}{5})$", xy=(0, 1.6), xytext=(10, 6),
                textcoords="offset points")  # fmt: skip
    ax.annotate(r"$(\frac{8}{3},\ 0,\ -\frac{1}{3})$", xy=(8 / 3, 0), xytext=(8, 10),
                textcoords="offset points", color=ORANGE)  # fmt: skip
    ax.text(1.75, 0.62, "$x_3 < 0$", color=MUTED)  # above the infeasible part
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    _halo(fig)
    return save_mpl(fig, "08-01-feasible-segment", outdir)


def lad_fit(outdir, lang="en"):
    """Typing speed y (wpm) against week x with two fitted lines (ex-08-01-lad-regression).

    Solid: least-absolute-deviation line y = 20/7 x + 45/7, which passes
    through the points of weeks 3 and 10 (ringed). Dashed: least-squares
    line (np.polyfit), about y = 3.170 x + 4.067.
    """
    x = np.arange(1, 11)
    y = np.array([5, 9, 15, 19, 21, 24, 26, 30, 31, 35])
    a_lad, b_lad = 20 / 7, 45 / 7
    a_ls, b_ls = np.polyfit(x, y, 1)

    fig, ax = plt.subplots(figsize=(5.0, 3.5))
    xx = np.array([0.5, 10.5])
    ax.plot(xx, a_lad * xx + b_lad, color=BLUE, linewidth=2,
            label=r"$\hat y = \frac{20}{7}x + \frac{45}{7}$")  # fmt: skip
    ax.plot(xx, a_ls * xx + b_ls, color=ORANGE, linestyle="--", linewidth=1.6,
            label=rf"$\hat y = {a_ls:.3f}x + {b_ls:.3f}$")  # fmt: skip
    ax.plot(x, y, "o", color=INK, markersize=5, zorder=5)
    # The two data points on the LAD line.
    ax.plot([3, 10], [15, 35], "o", markersize=11, markerfacecolor="none",
            markeredgecolor=BLUE, markeredgewidth=1.6, zorder=6)  # fmt: skip
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 40)
    ax.set_xticks(range(1, 11))
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.legend(loc="upper left")
    _halo(fig)
    return save_mpl(fig, "08-01-lad-fit", outdir)


FIGURES = {
    "08-01-capital-activities": capital_activities,
    "08-01-oilco-region": oilco_region,
    "08-01-basic-solutions-2d": basic_solutions_2d,
    "08-01-feasible-segment": feasible_segment,
    "08-01-lad-fit": lad_fit,
}
