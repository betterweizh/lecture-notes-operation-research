"""Figures of chapter 8.4, Revised and Dual Simplex, Sensitivity (Part VIII).

* ``08-04-reddy-mikks-objectives``: optimality cones of the paint-mix LP
  max c1 x1 + c2 x2 s.t. 6x1 + 4x2 <= 24, x1 + 2x2 <= 6, -x1 + x2 <= 1,
  x2 <= 2, x >= 0, with the objective gradients (5, 4), (3, 2), (8, 10) and
  (2, 5). Used in the cost-ranging example of the chapter.
"""

import math

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from .common import (
    BLUE,
    GREEN,
    INK,
    MUTED,
    ORANGE,
    PALETTE,
    lp_region,
    save_mpl,
)


def _unit(v, length):
    """Vector v rescaled to the given length (data units)."""
    n = math.hypot(*v)
    return (v[0] / n * length, v[1] / n * length)


def reddy_mikks_objectives(outdir, lang="en"):
    """Feasible region, optimality cones at (3, 1.5) and (2, 2), and gradients.

    The cone at a corner is spanned by the outward normals of its binding
    constraints: (6, 4) and (1, 2) at (3, 1.5); (1, 2) and (0, 1) at (2, 2).
    A gradient inside a cone makes that corner optimal.
    """
    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    # Constraints (all <=): paint-mix LP. Labels placed where they are clear.
    constraints = [
        {"a": (6, 4), "op": "<=", "b": 24, "label": r"$6x_1+4x_2=24$",
         "label_at": 0.86, "color": MUTED},
        {"a": (1, 2), "op": "<=", "b": 6, "label": r"$x_1+2x_2=6$",
         "label_at": 0.80, "color": MUTED},
        {"a": (-1, 1), "op": "<=", "b": 1, "label": r"$-x_1+x_2=1$",
         "label_at": 0.22, "color": MUTED},
        {"a": (0, 1), "op": "<=", "b": 2, "label": r"$x_2=2$",
         "label_at": 0.94, "color": MUTED},
    ]  # fmt: skip
    lp_region(ax, constraints, (0, 5.6), (0, 3.6), alpha=0.12)
    ax.set_aspect("equal")

    L = 1.15  # length of the cone edges and gradient arrows (data units)

    # Optimality cone at (3, 1.5): normals (6, 4) and (1, 2).
    p = (3.0, 1.5)
    e1, e2 = _unit((6, 4), L), _unit((1, 2), L)
    ax.add_patch(
        Polygon(
            [p, (p[0] + e1[0], p[1] + e1[1]), (p[0] + e2[0], p[1] + e2[1])],
            closed=True,
            color=BLUE,
            alpha=0.18,
            lw=0,
        )
    )
    # Optimality cone at (2, 2): normals (1, 2) and (0, 1).
    q = (2.0, 2.0)
    f1, f2 = _unit((1, 2), L), _unit((0, 1), L)
    ax.add_patch(
        Polygon(
            [q, (q[0] + f1[0], q[1] + f1[1]), (q[0] + f2[0], q[1] + f2[1])],
            closed=True,
            color=GREEN,
            alpha=0.18,
            lw=0,
        )
    )

    # Gradients from (3, 1.5): original (5, 4), (a) (3, 2), (b) (8, 10);
    # (c) (2, 5) drawn from (3, 1.5) dashed (outside the cone) and from (2, 2)
    # solid (inside the cone of the new optimum).
    arrows = [
        (p, (5, 4), INK, "-", r"$(5,4)$", (5, -2)),
        (p, (3, 2), PALETTE[3], "-", r"(a) $(3,2)$", (4, -18)),
        (p, (8, 10), ORANGE, "-", r"(b) $(8,10)$", (5, 2)),
        (p, (2, 5), MUTED, "--", None, None),
        # (c): right of the tip, clear of the lines through (2, 3)
        (q, (2, 5), GREEN, "-", r"(c) $(2,5)$", (8, -6)),
    ]
    for start, v, col, ls, label, off in arrows:
        d = _unit(v, L)
        end = (start[0] + d[0], start[1] + d[1])
        ax.annotate(
            "", xy=end, xytext=start,
            arrowprops=dict(arrowstyle="-|>", color=col, lw=1.6, linestyle=ls,
                            mutation_scale=12, shrinkA=0, shrinkB=0),
        )  # fmt: skip
        if label:
            ax.annotate(
                label,
                xy=end,
                xytext=off,
                textcoords="offset points",
                fontsize=9.5,
                color=col,
            )

    # Corners.
    ax.plot(*p, "o", color=INK, markersize=5, zorder=5)
    ax.plot(*q, "o", color=INK, markersize=5, zorder=5)
    ax.annotate(
        r"$(3,1.5)$",
        xy=p,
        xytext=(-46, -14),
        textcoords="offset points",
        fontsize=9.5,
    )
    ax.annotate(
        r"$(2,2)$",
        xy=q,
        xytext=(-40, -14),
        textcoords="offset points",
        fontsize=9.5,
    )
    return save_mpl(fig, "08-04-reddy-mikks-objectives", outdir)


FIGURES = {
    "08-04-reddy-mikks-objectives": reddy_mikks_objectives,
}
