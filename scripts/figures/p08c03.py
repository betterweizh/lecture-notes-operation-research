"""Figures of chapter 8.3, Duality (Part VIII).

* ``08-03-big-m-dual-region-1``: graphical solution of the dual of
  max 5x1 + 2x2 + 3x3 s.t. x1 + 5x2 + 2x3 = 30, x1 - 5x2 - 6x3 <= 40,
  x >= 0 (example ex-08-03-a4-q3).
* ``08-03-big-m-dual-region-2``: graphical solution of the dual of
  max x1 + 5x2 + 3x3 s.t. x1 + 2x2 + x3 = 3, 2x1 - x2 = 4, x >= 0
  (example ex-08-03-a4-q4); both dual variables are free.
"""

import matplotlib.pyplot as plt

from .common import INK, MUTED, ORANGE, draw_line, lp_region, save_mpl


def _mark_optimum(ax, p, text, xytext):
    """Mark the optimal dual point p; label it at xytext (data coordinates)
    with a thin leader line."""
    px, py = float(p[0]), float(p[1])
    ax.plot(px, py, "o", color=ORANGE, markersize=7, zorder=6, clip_on=False)
    ax.annotate(text, xy=(px, py), xytext=xytext, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=INK, lw=0.8, shrinkB=5))  # fmt: skip


def big_m_dual_region_1(outdir, lang="en"):
    """Dual region of ex-08-03-a4-q3: min w = 30 pi1 + 40 pi2.

    Dual constraints (one per primal variable x1, x2, x3):
    pi1 + pi2 >= 5, 5 pi1 - 5 pi2 >= 2, 2 pi1 - 6 pi2 >= 3, with pi1 free
    (equation row) and pi2 >= 0 (<= row of a max problem). The region lies
    in pi1 >= 33/8 > 0, so drawing only the first quadrant loses nothing.
    Vertices: (5, 0) with w = 150 (optimal) and (33/8, 7/8) with w = 635/4.
    """
    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    xlim, ylim = (0, 9), (0, 3.2)
    constraints = [
        {"a": (1, 1), "op": ">=", "b": 5, "label": "($x_1$)", "label_at": 0.62},
        {"a": (5, -5), "op": ">=", "b": 2, "label": "($x_2$)", "label_at": 0.3,
         "color": MUTED},
        {"a": (2, -6), "op": ">=", "b": 3, "label": "($x_3$)", "label_at": 0.75},
    ]  # fmt: skip
    vertices = lp_region(
        ax, constraints, xlim, ylim, mark_vertices=True, exact=True,
        axis_labels=(r"$\pi_1$", r"$\pi_2$"),
    )  # fmt: skip
    # The optimum minimizes w over the vertices (the region is unbounded,
    # but w grows along its recession direction (3, 1)).
    p_opt = min(vertices, key=lambda p: 30 * p[0] + 40 * p[1])
    assert p_opt == (5, 0)
    # Objective contour through the optimum, w = 150.
    draw_line(
        ax, 30, 40, 150, xlim, ylim, color=INK, linestyle="--",
        label=None,  # unlabeled: it runs almost parallel to line (x1);
        # the caption names it
    )  # fmt: skip
    _mark_optimum(ax, p_opt, r"$\boldsymbol{\pi}^* = (5,\ 0)$", (5.6, 0.3))
    # The other vertex, for comparison; label inside the region, below the
    # line (x3) that leaves the vertex to the upper right.
    ax.annotate(r"$(\frac{33}{8},\ \frac{7}{8})$, $w = 158.75$", xy=(33 / 8, 7 / 8),
                xytext=(12, -6), textcoords="offset points", ha="left",
                va="top", color=INK)  # fmt: skip
    return save_mpl(fig, "08-03-big-m-dual-region-1", outdir)


def big_m_dual_region_2(outdir, lang="en"):
    """Dual region of ex-08-03-a4-q4: min w = 3 pi1 + 4 pi2, pi free.

    Dual constraints: pi1 + 2 pi2 >= 1 (x1), 2 pi1 - pi2 >= 5 (x2),
    pi1 >= 3 (x3). Both primal rows are equations, so pi1 and pi2 are free
    and the window includes negative pi2. Vertices: (3, -1) with w = 5
    (optimal) and (3, 1) with w = 13.
    """
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    xlim, ylim = (0, 7.5), (-3, 4)
    constraints = [
        {
            "a": (1, 2),
            "op": ">=",
            "b": 1,
            "label": "($x_1$)",
            "label_at": 0.85,
        },
        {
            "a": (2, -1),
            "op": ">=",
            "b": 5,
            "label": "($x_2$)",
            "label_at": 0.72,
        },
        {
            "a": (1, 0),
            "op": ">=",
            "b": 3,
            "label": "($x_3$)",
            "label_at": 0.88,
        },
    ]
    vertices = lp_region(
        ax, constraints, xlim, ylim, nonneg=False, mark_vertices=True,
        exact=True, axis_labels=(r"$\pi_1$", r"$\pi_2$"),
    )  # fmt: skip
    # pi2 = 0 drawn as a thin reference line (pi2 is free).
    ax.axhline(0, color=MUTED, lw=0.8, zorder=1)
    p_opt = min(vertices, key=lambda p: 3 * p[0] + 4 * p[1])
    assert p_opt == (3, -1)
    draw_line(
        ax, 3, 4, 5, xlim, ylim, color=INK, linestyle="--",
        label=r"$w = 3\pi_1 + 4\pi_2 = 5$", label_at=0.17,
    )  # fmt: skip
    _mark_optimum(ax, p_opt, r"$\boldsymbol{\pi}^* = (3,\ -1)$", (4.4, -0.55))
    ax.annotate(r"$(3,\ 1)$, $w = 13$", xy=(3, 1), xytext=(-10, 10),
                textcoords="offset points", ha="right", va="bottom",
                color=INK)  # fmt: skip
    return save_mpl(fig, "08-03-big-m-dual-region-2", outdir)


FIGURES = {
    "08-03-big-m-dual-region-1": big_m_dual_region_1,
    "08-03-big-m-dual-region-2": big_m_dual_region_2,
}
