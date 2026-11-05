"""Figures of chapter 3.3 (The Simplex Method).

03-03-simplex-path: feasible region of the product-mix LP with the path of the
simplex method (Dantzig's rule from the slack basis), the optimal objective
contour and the infeasible basic solution met in the matrix-form example.
"""

from fractions import Fraction

import matplotlib.pyplot as plt

from .common import (
    BLUE,
    GREEN,
    INK,
    MUTED,
    ORANGE,
    PALETTE,
    draw_line,
    lp_region,
    save_mpl,
)


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


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


def simplex_path(outdir, lang="en"):
    """Product-mix LP max 500 x1 + 450 x2 with the simplex path.

    Used in 03-03-simplex-method (fig-03-03-simplex-path). Shows the feasible
    region of 6 x1 + 5 x2 <= 60, 10 x1 + 20 x2 <= 150, x1 <= 8, x >= 0,
    the vertices visited by the tableau method, x(0) = (0, 0) -> x(1) =
    (8, 0) -> x(2) = (8, 12/5) -> x(3) = (45/7, 30/7), arrows along the
    edges, the optimal contour z = 36000/7, and the infeasible basic
    solution (8, 7/2) obtained when the wrong variable leaves.
    """
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    xlim, ylim = (0, 12), (0, 10)

    # Constraints a1*x1 + a2*x2 <= b; labels placed along each line.
    constraints = [
        {"a": (6, 5), "op": "<=", "b": 60, "label_at": 0.13},  # machine hours
        {
            "a": (10, 20),
            "op": "<=",
            "b": 150,
            "label_at": 0.2,
        },  # storage (m^3)
        {
            "a": (1, 0),
            "op": "<=",
            "b": 8,
            "label": "$x_1 = 8$",
            "label_at": 0.88,
        },
    ]
    # Constraint lines in blue, green and purple; orange is kept for the path.
    lp_region(
        ax,
        constraints,
        xlim,
        ylim,
        exact=True,
        colors=[BLUE, GREEN, PALETTE[3]],
    )

    # Optimal objective contour 500 x1 + 450 x2 = 36000/7 (dashed, muted).
    z_opt = Fraction(36000, 7)
    draw_line(
        ax,
        500,
        450,
        float(z_opt),
        xlim,
        ylim,
        color=MUTED,
        linestyle="--",
        linewidth=1.1,
    )  # fmt: skip  (explained in the caption; a label would crowd the edge)

    # Vertices visited by the simplex method, in order.
    path = [(0, 0), (8, 0), (8, 12 / 5), (45 / 7, 30 / 7)]
    names = [
        r"$\mathbf{x}^{(0)} = (0, 0)$",
        r"$\mathbf{x}^{(1)} = (8, 0)$",
        r"$\mathbf{x}^{(2)} = (8, \frac{12}{5})$",
        r"$\mathbf{x}^{(3)} = (\frac{45}{7}, \frac{30}{7})$",
    ]
    # Label offsets (points) chosen to keep labels clear of lines.
    offsets = [(8, 8), (-8, 8), (8, -2), (-10, -26)]
    aligns = ["left", "right", "left", "right"]

    # Arrows along the edges between successive vertices.
    for p, q in zip(path, path[1:]):
        ax.annotate(
            "", xy=q, xytext=p,
            arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2.2,
                            shrinkA=4, shrinkB=6, mutation_scale=14),
            zorder=4,
        )  # fmt: skip

    # Visited vertices (filled) and their labels.
    for k, (p, name) in enumerate(zip(path, names)):
        color = ORANGE if k == len(path) - 1 else INK
        # not clipped: vertices on the axes show as whole dots
        ax.plot(*p, marker="o", markersize=6.5, color=color, zorder=5,
                clip_on=False)  # fmt: skip
        ax.annotate(
            name, xy=p, xytext=offsets[k], textcoords="offset points",
            ha=aligns[k], va="bottom", fontsize=10.5, color=INK, zorder=6,
        )  # fmt: skip

    # Infeasible basic solution (8, 7/2): storage line meets x1 = 8.
    ax.plot(
        8,
        3.5,
        marker="o",
        markersize=6.5,
        markerfacecolor="none",
        markeredgecolor=BLUE,
        markeredgewidth=1.4,
        zorder=5,
    )
    ax.annotate(
        _t(lang, r"$(8, \frac{7}{2})$ infeasible", r"$(8, \frac{7}{2})$ 不可行"),
        xy=(8, 3.5), xytext=(8, 4), textcoords="offset points",
        ha="left", va="bottom", fontsize=10, color=BLUE,
    )  # fmt: skip

    # The unvisited vertex (0, 15/2) is marked lightly for reference.
    ax.plot(0, 7.5, marker="o", markersize=5, color=MUTED, zorder=5,
            clip_on=False)  # fmt: skip
    _halo(fig)
    return save_mpl(fig, "03-03-simplex-path", outdir)


FIGURES = {
    "03-03-simplex-path": simplex_path,
}
