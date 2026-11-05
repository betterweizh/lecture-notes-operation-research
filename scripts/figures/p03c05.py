"""Figures of chapter 3.5 (Duality Theory).

* ``03-05-weak-duality-line``: objective values of feasible primal and dual
  solutions of the diet problem on one number line (weak duality).
* ``03-05-diet-dual-region``: feasible region of the pill-maker's (dual)
  problem with its optimal prices; the tight constraints are those of the
  foods used in the optimal diet (complementary slackness).
* ``03-05-cs-cone``: the product-mix LP at its optimal vertex; the objective
  gradient lies in the cone spanned by the normals of the active
  constraints, with the optimal dual values as weights.
"""

import math
from fractions import Fraction

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from .common import (
    BLUE,
    GREEN,
    INK,
    MUTED,
    ORANGE,
    draw_line,
    lp_region,
    save_mpl,
)


def weak_duality_line(outdir, lang="en"):
    """Number line of primal and dual objective values (diet problem).

    Used in section "Weak duality" of chapter 3.5. Every dual feasible price
    vector gives a value w = b^T pi to the left of every primal feasible
    cost z = c^T x; the two sets meet at the common optimum 179 cents.
    All values were checked in scratch/E1/verify_ch11.py.
    """
    # Dual feasible prices pi -> w = 9 pi1 + 19 pi2 (cents):
    # (11, 0) -> 99, (5, 6) -> 159, (0, 9) -> 171.
    dual_vals = [99, 159, 171]
    # Primal feasible diets x -> z = c^T x (cents):
    # x5 = 3, x6 = 5 -> 191; only food 6 (9.5 kg) -> 209; only food 5 (9 kg) -> 243.
    primal_vals = [191, 209, 243]
    opt = 179  # z* = w* = 179 cents

    fig, ax = plt.subplots(figsize=(6.4, 1.9))
    ax.axis("off")
    lo, hi = 85, 257
    ax.set_xlim(lo - 3, hi + 3)
    ax.set_ylim(-1.15, 1.35)

    # The number line itself, with an arrow head on the right (values grow).
    ax.annotate(
        "", xy=(hi, 0), xytext=(lo, 0),
        arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2),
    )  # fmt: skip

    # Dual values (blue circles) and primal values (orange squares),
    # each labeled with its value below the line.
    ax.plot(dual_vals, [0] * 3, "o", color=BLUE, markersize=7, zorder=3)
    ax.plot(primal_vals, [0] * 3, "s", color=ORANGE, markersize=7, zorder=3)
    for v in dual_vals + primal_vals:
        ax.text(
            v, -0.32, str(v), ha="center", va="top", fontsize=10, color=INK
        )

    # The common optimum: green diamond, label above.
    ax.plot([opt], [0], "D", color=GREEN, markersize=8, zorder=4)
    ax.plot([opt, opt], [0, 0.55], color=GREEN, lw=1.0)
    ax.text(opt, 0.62, r"$z^* = w^* = 179$", ha="center", va="bottom",
            color=GREEN, fontsize=11)  # fmt: skip

    # Brackets above each set: dual values increase toward the optimum,
    # primal values decrease toward it.
    ax.annotate(
        "", xy=(opt - 4, 0.32), xytext=(dual_vals[0], 0.32),
        arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.4),
    )  # fmt: skip
    ax.text(127, 0.40, r"$w = \mathbf{b}^\mathsf{T}\boldsymbol{\pi}$",
            ha="center", va="bottom", color=BLUE, fontsize=11)  # fmt: skip
    ax.annotate(
        "", xy=(opt + 4, 0.32), xytext=(primal_vals[-1], 0.32),
        arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.4),
    )  # fmt: skip
    ax.text(222, 0.40, r"$z = \mathbf{c}^\mathsf{T}\mathbf{x}$",
            ha="center", va="bottom", color=ORANGE, fontsize=11)  # fmt: skip
    return save_mpl(fig, "03-05-weak-duality-line", outdir)


def diet_dual_region(outdir, lang="en"):
    """Feasible region of the dual of the diet problem, optimum (3, 8).

    Used in section "Complementary slackness" of chapter 3.5. Constraint
    line (j) is the price limit of food j: pi^T A_j <= c_j. Foods 1 and 2
    (pi1 <= 35, pi2 <= 30) lie outside the window. Only the lines of
    foods 5 and 6 pass through pi* = (3, 8), the foods bought in the
    optimal diet x* = (0, 0, 0, 0, 5, 2).
    """
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    xlim, ylim = (0, 26), (0, 21)
    # Dual constraints of foods 3-6 (foods 1, 2 fall outside the window).
    constraints = [
        {
            "a": (2, 3),
            "op": "<=",
            "b": 60,
            "label": "(3)",
            "label_at": 0.55,
            "color": MUTED,
        },
        {
            "a": (2, 1),
            "op": "<=",
            "b": 50,
            "label": "(4)",
            "label_at": 0.35,
            "color": MUTED,
        },
        {
            "a": (1, 3),
            "op": "<=",
            "b": 27,
            "label": "(5)",
            "label_at": 0.62,
            "color": ORANGE,
        },
        {
            "a": (1, 1),
            "op": "<=",
            "b": 11,
            "label": "(6)",
            "label_at": 0.7,
            "color": BLUE,
        },
    ]
    vertices = lp_region(
        ax, constraints, xlim, ylim, mark_vertices=True, exact=True,
        axis_labels=(r"$\pi_1$", r"$\pi_2$"),
    )  # fmt: skip
    # Optimal prices maximise w = 9 pi1 + 19 pi2 over the vertices.
    p_opt = max(vertices, key=lambda p: 9 * p[0] + 19 * p[1])
    w_opt = 9 * p_opt[0] + 19 * p_opt[1]  # 179
    assert p_opt == (3, 8) and w_opt == 179

    # Objective contour through the optimum and the direction of increase b.
    draw_line(
        ax, 9, 19, w_opt, xlim, ylim, color=GREEN, linestyle="--",
        label=r"$w = 9\pi_1 + 19\pi_2 = 179$", label_at=0.83,
    )  # fmt: skip
    px, py = float(p_opt[0]), float(p_opt[1])
    s = 0.22  # arrow scale for the gradient b = (9, 19)
    ax.annotate(
        "", xy=(px + s * 9, py + s * 19), xytext=(px, py),
        arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.5),
    )  # fmt: skip
    ax.text(px + s * 9 + 0.4, py + s * 19, r"$\mathbf{b}$", color=GREEN,
            va="center", fontsize=12)  # fmt: skip

    # Mark pi* and label it in the empty area between lines (3) and (5),
    # with a thin leader line (the area left of pi* holds the axis ticks).
    ax.plot(px, py, "o", color=ORANGE, markersize=7, zorder=5)
    ax.annotate(
        r"$\boldsymbol{\pi}^* = (3,\ 8)$", xy=(px, py), xytext=(8.5, 9.4),
        ha="left", va="center",
        arrowprops=dict(arrowstyle="-", color=INK, lw=0.8, shrinkB=5),
    )  # fmt: skip
    # Integer ticks are enough on both axes.
    ax.set_xticks(range(0, 27, 5))
    ax.set_yticks(range(0, 21, 5))
    return save_mpl(fig, "03-05-diet-dual-region", outdir)


def cs_cone(outdir, lang="en"):
    """Objective gradient in the cone of active-constraint normals.

    Used in section "Complementary slackness" of chapter 3.5. Product-mix
    LP max 500 x1 + 450 x2 s.t. 6x1 + 5x2 <= 60, 10x1 + 20x2 <= 150,
    x1 <= 8. At x* = (45/7, 30/7) rows 1 and 2 are active and
    c = (550/7) a1 + (20/7) a2: c lies in the cone spanned by a1, a2.
    """
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    xlim, ylim = (0, 11.5), (0, 9)
    constraints = [
        {"a": (6, 5), "op": "<=", "b": 60, "label": "(1)", "label_at": 0.2},
        {"a": (10, 20), "op": "<=", "b": 150, "label": "(2)", "label_at": 0.2},
        {"a": (1, 0), "op": "<=", "b": 8, "label": "(3)", "label_at": 0.12},
    ]
    vertices = lp_region(
        ax, constraints, xlim, ylim, mark_vertices=True, exact=True
    )
    x_opt = max(vertices, key=lambda p: 500 * p[0] + 450 * p[1])
    assert x_opt == (Fraction(45, 7), Fraction(30, 7))
    x0, y0 = float(x_opt[0]), float(x_opt[1])

    def unit(v):
        n = math.hypot(*v)
        return v[0] / n, v[1] / n

    # Normals of the active rows and the objective gradient (scaled for display).
    a1, a2, c = unit((6, 5)), unit((10, 20)), unit((500, 450))
    r = 3.3  # radius of the shaded cone
    # Shaded cone {lambda1 a1 + lambda2 a2 : lambda >= 0}, drawn as a wedge.
    ang1, ang2 = math.atan2(a1[1], a1[0]), math.atan2(a2[1], a2[0])
    arc = [ang1 + (ang2 - ang1) * k / 30 for k in range(31)]
    wedge = [(x0, y0)] + [
        (x0 + r * math.cos(t), y0 + r * math.sin(t)) for t in arc
    ]
    ax.add_patch(
        Polygon(wedge, closed=True, color=GREEN, alpha=0.18, lw=0, zorder=2)
    )

    # Arrows: a1 and a2 (ink), c (orange, a little longer).
    for vec, length, color, label, off in [
        (a1, 2.4, INK, r"$\mathbf{a}_1$", (0.15, -0.35)),
        (a2, 2.4, INK, r"$\mathbf{a}_2$", (-0.55, 0.05)),
        (c, 3.0, ORANGE, r"$\mathbf{c}$", (0.18, -0.1)),
    ]:
        tip = (x0 + length * vec[0], y0 + length * vec[1])
        ax.annotate("", xy=tip, xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=1.6), zorder=5)  # fmt: skip
        ax.text(
            tip[0] + off[0], tip[1] + off[1], label, color=color, fontsize=12
        )

    # Objective contour through x*.
    draw_line(ax, 500, 450, Fraction(36000, 7), xlim, ylim, color=MUTED,
              linestyle="--", linewidth=1.2)  # fmt: skip
    ax.plot(x0, y0, "o", color=ORANGE, markersize=7, zorder=6)
    ax.annotate(
        r"$\mathbf{x}^* = (\frac{45}{7},\ \frac{30}{7})$", xy=(x0, y0), xytext=(-10, -14),
        textcoords="offset points", ha="right", va="top",
    )  # fmt: skip
    # The decomposition c = pi1 a1 + pi2 a2 (unscaled normals).
    ax.text(0.3, 0.55,
            r"$\mathbf{c} = \frac{550}{7}\,\mathbf{a}_1 + \frac{20}{7}\,\mathbf{a}_2$",
            fontsize=11, color=INK)  # fmt: skip
    return save_mpl(fig, "03-05-cs-cone", outdir)


FIGURES = {
    "03-05-weak-duality-line": weak_duality_line,
    "03-05-diet-dual-region": diet_dual_region,
    "03-05-cs-cone": cs_cone,
}
