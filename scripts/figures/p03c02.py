"""Figures of chapter 3.2, Geometry of Linear Programs (all matplotlib).

* ``03-02-product-mix-region``: feasible region of the product-mix LP with
  three objective contours and the optimal vertex (graphical method);
* ``03-02-graphical-cases``: four small LPs showing a unique optimum,
  alternative optima, an unbounded LP and an infeasible LP;
* ``03-02-half-space``: the closed half-space x1 + x2 >= 1 and its boundary;
* ``03-02-product-mix-basic-solutions``: all basic solutions of the product
  mix in standard form, feasible (filled) and infeasible (hollow).

Labels are the book's symbols and numbers (language neutral).
"""

from fractions import Fraction as F

import matplotlib.pyplot as plt

from .common import (
    BLUE,
    INK,
    MUTED,
    ORANGE,
    draw_line,
    lp_region,
    save_mpl,
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


# ---------------------------------------------------------------------------
# Data: the product-mix LP  max 500 x1 + 450 x2
#   6 x1 + 5 x2 <= 60 (machine h), 10 x1 + 20 x2 <= 150 (storage m3),
#   x1 <= 8 (demand, hundreds of cases), x >= 0.
# ---------------------------------------------------------------------------
PM_C = (500, 450)
PM_X_OPT = (F(45, 7), F(30, 7))  # optimal vertex, z* = 36000/7


def _pm_constraints(label_at, machine_label="$6x_1 + 5x_2 = 60$"):
    """Product-mix constraints in lp_region's dict form; label positions vary.

    ``machine_label=None`` leaves the machine line unlabeled (it is then
    labeled by hand, clear of a nearly parallel objective contour).
    """
    return [
        {
            "a": (6, 5),
            "op": "<=",
            "b": 60,
            "label": machine_label,
            "label_at": label_at[0],
        },
        {"a": (10, 20), "op": "<=", "b": 150, "label_at": label_at[1]},
        {
            "a": (1, 0),
            "op": "<=",
            "b": 8,
            "label": "$x_1 = 8$",
            "label_at": label_at[2],
        },
    ]


def _dot(ax, x, y, *, filled=True, color=ORANGE, size=6.5, z=6):
    """Mark a point: filled dot (feasible) or hollow dot (infeasible).

    Not clipped, so that points on the axes show as whole dots.
    """
    ax.plot(
        [float(x)], [float(y)], marker="o", markersize=size, zorder=z,
        color=color, markerfacecolor=color if filled else "white", markeredgewidth=1.6,
        clip_on=False,
    )  # fmt: skip


def product_mix_region(outdir, lang="en"):
    """Graphical solution of the product mix (ex-03-02-graphical).

    Shades the pentagon, draws the contours z = 2250, 4000 and 36000/7 and
    marks the optimal vertex x* = (45/7, 30/7).
    """
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    xlim, ylim = (0, 13), (0, 10)
    vertices = lp_region(ax, _pm_constraints((0.72, 0.85, 0.93), machine_label=None),
                         xlim, ylim, mark_vertices=True, exact=True)  # fmt: skip
    # The machine line and the optimal contour are nearly parallel, so both
    # get hand-placed labels: the machine line right of its x1-intercept
    # (10, 0), the contour just right of where it leaves the top edge.
    ax.text(10.4, 0.45, "$6x_1 + 5x_2 = 60$", color=BLUE, fontsize=10)

    # Objective contours 500 x1 + 450 x2 = z (slope -10/9), dashed and muted.
    for z, at in [(2250, 0.5), (4000, 0.5)]:
        draw_line(ax, *PM_C, z, xlim, ylim, color=MUTED, linestyle="--",
                  linewidth=1.2, label=f"$z = {z}$", label_at=at, label_color=MUTED)  # fmt: skip
    # The optimal contour, darker.
    z_opt = PM_C[0] * PM_X_OPT[0] + PM_C[1] * PM_X_OPT[1]
    assert z_opt == F(36000, 7) and PM_X_OPT in vertices
    draw_line(
        ax, *PM_C, z_opt, xlim, ylim, color=INK, linestyle="--", linewidth=1.4
    )
    ax.text(2.25, 9.45, "$z = 36000/7$", color=INK, fontsize=10, ha="left")

    # Optimal vertex.
    _dot(ax, *PM_X_OPT)
    ax.annotate(r"$\mathbf{x}^* = (\frac{45}{7},\ \frac{30}{7})$", xy=(45 / 7, 30 / 7),
                xytext=(14, 10), textcoords="offset points")  # fmt: skip
    _halo(fig)
    return save_mpl(fig, "03-02-product-mix-region", outdir)


def graphical_cases(outdir, lang="en"):
    """The four outcomes of an LP, one small example each (sec-03-02-graphical-cases).

    (a) max x1 + x2 over x1 + 2x2 <= 6, 2x1 + x2 <= 6: unique optimum (2, 2).
    (b) max x1 + 2x2 over the same region: the edge (0,3)-(2,2) is optimal.
    (c) max x1 + x2 over -x1 + x2 <= 2, x1 - 2x2 <= 2: unbounded along (2, 1).
    (d) x1 + x2 <= 2 and x1 + 2x2 >= 6: no common point (infeasible).
    """
    fig, axs = plt.subplots(2, 2, figsize=(8.0, 6.6))
    small = 9  # label size inside the panels
    region = [{"a": (1, 2), "op": "<=", "b": 6, "label_at": 0.75},
              {"a": (2, 1), "op": "<=", "b": 6, "label_at": 0.15}]  # fmt: skip
    win = ((0, 6.5), (0, 6.5))

    # (a) unique optimum: contour x1 + x2 = 4 touches only (2, 2).
    ax = axs[0, 0]
    lp_region(ax, region, *win, label_size=small)
    draw_line(ax, 1, 1, 4, *win, color=MUTED, linestyle="--", linewidth=1.3,
              label="$z = 4$", label_at=0.85, label_color=MUTED, label_size=small)  # fmt: skip
    _dot(ax, 2, 2)
    ax.annotate("$(2,\\ 2)$", xy=(2, 2), xytext=(-8, -16), textcoords="offset points",
                ha="right", fontsize=small)  # fmt: skip
    ax.set_title("(a)", loc="left")

    # (b) alternative optima: contour x1 + 2x2 = 6 lies on the first edge.
    ax = axs[0, 1]
    lp_region(ax, region, *win, label_size=small)
    ax.plot(
        [0, 2],
        [3, 2],
        color=ORANGE,
        linewidth=4,
        solid_capstyle="round",
        zorder=5,
    )
    _dot(ax, 0, 3, size=5.5)
    _dot(ax, 2, 2, size=5.5)
    ax.annotate("$z = 6$", xy=(1, 2.5), xytext=(-22, 14), textcoords="offset points",
                color=ORANGE, fontsize=small)  # fmt: skip
    ax.set_title("(b)", loc="left")

    # (c) unbounded: the region continues to the upper right; z = x1 + x2
    # grows along the ray t (2, 1) from the origin.
    ax = axs[1, 0]
    win_c = ((0, 8), (0, 6))
    lp_region(ax, [{"a": (-1, 1), "op": "<=", "b": 2, "label_at": 0.3},
                   {"a": (1, -2), "op": "<=", "b": 2, "label_at": 0.55}],
              *win_c, label_size=small)  # fmt: skip
    for z in (2, 5, 8):
        draw_line(
            ax, 1, 1, z, *win_c, color=MUTED, linestyle="--", linewidth=1.0
        )
    ax.annotate("", xy=(6, 3), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color=ORANGE, linewidth=2))  # fmt: skip
    ax.text(
        1.0, 1.25, "$t\\,(2,\\ 1)$", color=ORANGE, fontsize=small, ha="left"
    )
    ax.set_title("(c)", loc="left")

    # (d) infeasible: shade each half-plane in its own color; they do not meet.
    ax = axs[1, 1]
    win_d = ((0, 6.5), (0, 4.5))
    lp_region(ax, [{"a": (1, 1), "op": "<=", "b": 2, "label_at": 0.35}], *win_d,
              fill=BLUE, colors=[BLUE], label_size=small)  # fmt: skip
    lp_region(ax, [{"a": (1, 2), "op": ">=", "b": 6, "label_at": 0.62}], *win_d,
              fill=ORANGE, colors=[ORANGE], label_size=small)  # fmt: skip
    ax.set_title("(d)", loc="left")

    fig.tight_layout()
    _halo(fig)
    return save_mpl(fig, "03-02-graphical-cases", outdir)


def half_space(outdir, lang="en"):
    """The closed half-space x1 + x2 >= 1 (shaded) and the hyperplane x1 + x2 = 1.

    Used in Sec. "Half-spaces and polyhedra" of chapter 3.2. The arrow is the
    normal vector (1, 1), drawn from a point of the line; it points into
    the half-space.
    """
    fig, ax = plt.subplots(figsize=(4.4, 3.8))
    xlim, ylim = (0, 2.2), (0, 2.2)
    lp_region(
        ax, [{"a": (1, 1), "op": ">=", "b": 1, "label": None}], xlim, ylim
    )
    # Label the boundary below the line, in the unshaded corner.
    ax.text(0.06, 0.28, "$x_1 + x_2 = 1$", color=BLUE)
    # The two points that fix the line.
    for p, off, ha in [((0, 1), (8, 4), "left"), ((1, 0), (4, 8), "left")]:
        ax.plot(*p, marker="o", color=INK, markersize=5, zorder=5)
        ax.annotate(f"$({p[0]},\\ {p[1]})$", xy=p, xytext=off, textcoords="offset points",
                    ha=ha)  # fmt: skip
    # Normal vector (1, 1) scaled to length ~0.5, from (0.5, 0.5).
    ax.annotate("", xy=(0.85, 0.85), xytext=(0.5, 0.5),
                arrowprops=dict(arrowstyle="-|>", color=ORANGE, linewidth=2))  # fmt: skip
    ax.text(0.9, 0.78, "$(1,\\ 1)$", color=ORANGE)
    ax.text(1.35, 1.5, "$x_1 + x_2 \\geq 1$", color=BLUE, fontsize=12)
    _halo(fig)
    return save_mpl(fig, "03-02-half-space", outdir)


def product_mix_basic_solutions(outdir, lang="en"):
    """All nine basic solutions of the product mix in standard form.

    Each basic solution is the intersection of two of the five lines
    6x1 + 5x2 = 60, 10x1 + 20x2 = 150, x1 = 8, x1 = 0, x2 = 0. Filled dots:
    the five basic feasible solutions (the corners); hollow dots: the four
    infeasible ones. Used with tab-03-02-product-mix-basic.
    """
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    xlim, ylim = (0, 16.5), (0, 13.5)
    lp_region(ax, _pm_constraints((0.3, 0.72, 0.93)), xlim, ylim, alpha=0.12)

    # Feasible basic solutions (corners) with label offsets in points.
    feasible = [
        ((0, 0), "$(0,\\ 0)$", (6, 6)),
        ((8, 0), "$(8,\\ 0)$", (-36, 6)),
        ((8, F(12, 5)), "$(8,\\ \\frac{12}{5})$", (7, -10)),
        (
            (F(45, 7), F(30, 7)),
            "$(\\frac{45}{7},\\ \\frac{30}{7})$",
            (-58, -16),
        ),
        ((0, F(15, 2)), "$(0,\\ \\frac{15}{2})$", (7, 4)),
    ]
    # Infeasible basic solutions: intersections outside the region.
    infeasible = [
        ((8, F(7, 2)), "$(8,\\ \\frac{7}{2})$", (6, 13)),
        ((10, 0), "$(10,\\ 0)$", (-6, 10)),
        ((15, 0), "$(15,\\ 0)$", (5, 6)),  # right of the dot
        ((0, 12), "$(0,\\ 12)$", (8, -4)),
    ]
    for (x, y), lab, off in feasible:
        _dot(ax, x, y, filled=True, color=BLUE)
        ax.annotate(lab, xy=(float(x), float(y)), xytext=off, textcoords="offset points",
                    fontsize=10)  # fmt: skip
    for (x, y), lab, off in infeasible:
        _dot(ax, x, y, filled=False, color=ORANGE)
        ax.annotate(lab, xy=(float(x), float(y)), xytext=off, textcoords="offset points",
                    fontsize=10, color=ORANGE)  # fmt: skip
    _halo(fig)
    return save_mpl(fig, "03-02-product-mix-basic-solutions", outdir)


FIGURES = {
    "03-02-product-mix-region": product_mix_region,
    "03-02-graphical-cases": graphical_cases,
    "03-02-half-space": half_space,
    "03-02-product-mix-basic-solutions": product_mix_basic_solutions,
}
