"""Figures of chapter 7.3, Piecewise-Linear Approximation.

* ``07-03-interpolation``: a convex function (BPR total travel time of route 1,
  Example ex-07-03-congestion) and a concave one (plant cost q^0.6, Example
  ex-07-03-economies) with their piecewise-linear interpolants;
* ``07-03-sos2``: two representations of q = 50 as weights on breakpoints,
  adjacent (SOS2) and non-adjacent, with a bar chart of the weights;
* ``07-03-triangles``: a triangulated grid with the six triangles around one
  grid point (triangle method for functions of two variables).
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

from .common import BLUE, GREEN, INK, MUTED, ORANGE, save_mpl

# Route 1 of the congestion example: free-flow time 10 min, capacity 3000 veh/h
T0, CAP = 10.0, 3000.0
# Breakpoints of the plant-cost example (ML/d)
Q_BREAKS = np.array([0, 20, 40, 60, 80], float)


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


def total_time(v):
    """T_1(v) = v t_1(v) / 60 in veh-h/h with the BPR function."""
    return v * T0 * (1 + 0.15 * (v / CAP) ** 4) / 60.0


def plant_cost(q):
    """C_1(q) = q^0.6 in $ million per year (q in ML/d)."""
    return q**0.6


def _q_label(lang):
    """Axis title of the plant capacity q (megalitres per day)."""
    return _t(lang, "$q$ (ML/d)", "$q$（百万升/天）")


def _cost_label(lang):
    """Axis title of the plant cost."""
    return _t(lang, "\\$ million / year", "百万美元/年")


def interpolation(outdir, lang="en"):
    """Convex and concave functions with their interpolants.

    Used in Section sec-07-03-pwl. Left: T_1 with breakpoints every 1000 veh/h
    (interpolant above the curve). Right: C_1 with breakpoints every
    20 ML/d (interpolant below the curve) and the chord from 0 to 80, the
    convex envelope used by an LP without adjacency.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 3.6))
    interpolant = _t(lang, "interpolant", "插值函数")

    # --- convex: BPR total travel time --------------------------------
    v = np.linspace(0, 6000, 400)
    vb = np.arange(0, 6001, 1000.0)
    ax1.plot(v, total_time(v), color=BLUE, linewidth=2.2, label="$T_1(v)$")
    ax1.plot(
        vb,
        total_time(vb),
        color=ORANGE,
        linewidth=1.6,
        marker="o",
        markersize=4.5,
        label=interpolant,
    )
    ax1.set_xlabel(_t(lang, "$v$ (veh/h)", "$v$（辆/小时）"))
    ax1.set_ylabel(_t(lang, "veh-h/h", "辆·小时/小时"))
    ax1.set_title(_t(lang, "convex", "凸函数"), fontsize=11)
    ax1.set_xlim(0, 6200)
    ax1.set_ylim(0, 3600)
    ax1.legend(loc="upper left", frameon=False, fontsize=10)

    # --- concave: plant cost ------------------------------------------
    q = np.linspace(0, 80, 400)
    ax2.plot(q, plant_cost(q), color=BLUE, linewidth=2.2, label="$C_1(q)$")
    ax2.plot(
        Q_BREAKS,
        plant_cost(Q_BREAKS),
        color=ORANGE,
        linewidth=1.6,
        marker="o",
        markersize=4.5,
        label=interpolant,
    )
    # chord from (0, 0) to (80, C(80)): convex envelope of the interpolant
    ax2.plot(
        [0, 80],
        [0, plant_cost(80.0)],
        color=MUTED,
        linestyle="--",
        linewidth=1.4,
        label=_t(lang, "convex envelope", "凸包络"),
    )
    ax2.set_xlabel(_q_label(lang))
    ax2.set_ylabel(_cost_label(lang))
    ax2.set_title(_t(lang, "concave", "凹函数"), fontsize=11)
    ax2.set_xlim(0, 84)
    ax2.set_ylim(0, 15.5)
    ax2.legend(loc="lower right", frameon=False, fontsize=10)

    fig.tight_layout()
    return save_mpl(fig, "07-03-interpolation", outdir)


def sos2(outdir, lang="en"):
    """Adjacent (SOS2) and non-adjacent weights for q = 50.

    Used in Section sec-07-03-examples. Adjacent weights alpha_3 = alpha_4 = 0.5
    give f^a(50) = 10.41 (true 10.46); non-adjacent weights alpha_1 = 0.375,
    alpha_5 = 0.625 give the chord value 8.66.
    """
    fq = plant_cost(Q_BREAKS)
    adj = np.array([0, 0, 0.5, 0.5, 0])  # SOS2 weights of q = 50
    non = np.array([0.375, 0, 0, 0, 0.625])  # non-adjacent weights of q = 50
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(8.6, 3.6), gridspec_kw={"width_ratios": [1.25, 1]}
    )

    # --- left: curve, interpolant and the two represented points ---------
    q = np.linspace(0, 80, 400)
    ax1.plot(q, plant_cost(q), color=BLUE, linewidth=2.0, label="$C_1(q)$")
    ax1.plot(
        Q_BREAKS,
        fq,
        color=INK,
        linewidth=1.0,
        marker="o",
        markersize=4,
        label=_t(lang, "interpolant", "插值函数"),
    )
    ax1.plot([0, 80], [0, fq[-1]], color=MUTED, linestyle="--", linewidth=1.2)
    p_adj = (50.0, adj @ fq)  # (50, 10.41)
    p_non = (50.0, non @ fq)  # (50, 8.66)
    ax1.plot(*p_adj, marker="o", markersize=7, color=GREEN, zorder=5)
    ax1.plot(*p_non, marker="s", markersize=7, color=ORANGE, zorder=5)
    ax1.annotate(
        _t(lang, "adjacent: ", "相邻：") + f"{p_adj[1]:.2f}",
        xy=p_adj,
        xytext=(-8, 10),
        textcoords="offset points",
        ha="right",
        va="bottom",
        fontsize=10,
        color=GREEN,
    )
    ax1.annotate(
        _t(lang, "non-adjacent: ", "不相邻：") + f"{p_non[1]:.2f}",
        xy=p_non,
        xytext=(8, -6),
        textcoords="offset points",
        ha="left",
        va="top",
        fontsize=10,
        color=ORANGE,
    )
    ax1.axvline(
        50, color=MUTED, linewidth=0.8, linestyle=(0, (2, 3)), zorder=0
    )
    ax1.set_xticks(Q_BREAKS)
    ax1.set_xlabel(_q_label(lang))
    ax1.set_ylabel(_cost_label(lang))
    ax1.set_xlim(0, 84)
    ax1.set_ylim(0, 15.5)
    ax1.legend(loc="upper left", frameon=False, fontsize=10)

    # --- right: the weights alpha_i -------------------------------------
    # The two weight vectors are never positive at the same index, so each
    # bar is centered on its own tick (no side-by-side offset); a zero weight
    # draws nothing.
    idx = np.arange(1, 6)
    assert not np.any((adj > 0) & (non > 0))
    ax2.bar(idx, adj, width=0.5, color=GREEN,
            label=_t(lang, "SOS2 (adjacent)", "SOS2（相邻）"))  # fmt: skip
    ax2.bar(idx, non, width=0.5, color=ORANGE,
            label=_t(lang, "non-adjacent", "不相邻"))  # fmt: skip
    ax2.set_xticks(idx)
    ax2.set_xticklabels([f"$\\alpha_{i}$" for i in idx])
    ax2.set_ylim(0, 1.0)
    ax2.set_ylabel(_t(lang, "weight", "权重"))
    ax2.legend(loc="upper center", frameon=False, fontsize=10)
    ax2.grid(axis="x", visible=False)

    fig.tight_layout()
    return save_mpl(fig, "07-03-sos2", outdir)


def triangles(outdir, lang="en"):
    """Triangulated grid and the six triangles around one grid point.

    Used in Section sec-07-03-two-variables (triangle method). Rectangles
    (i, j) are split by the diagonal from (a_i, b_j) to (a_{i+1}, b_{j+1});
    the grid point (a_2, b_2) (0-based (1, 1)) belongs to six triangles.
    """
    n, m = 4, 4  # grid points per axis
    fig, ax = plt.subplots(figsize=(4.4, 3.9))
    # grid lines and diagonals
    for i in range(n):
        ax.plot([i, i], [0, m - 1], color=MUTED, linewidth=1.0)
    for j in range(m):
        ax.plot([0, n - 1], [j, j], color=MUTED, linewidth=1.0)
    for i in range(n - 1):
        for j in range(m - 1):
            ax.plot([i, i + 1], [j, j + 1], color=MUTED, linewidth=1.0)

    # six triangles around the grid point (1, 1); upper triangle of a
    # rectangle (i, j) has vertices (i,j), (i+1,j+1), (i,j+1); lower has
    # (i,j), (i+1,j), (i+1,j+1)
    def upper(i, j):
        return [(i, j), (i + 1, j + 1), (i, j + 1)]

    def lower(i, j):
        return [(i, j), (i + 1, j), (i + 1, j + 1)]

    c = (1, 1)
    six = [
        upper(1, 1),
        lower(1, 1),
        upper(0, 0),
        lower(0, 0),
        upper(1, 0),
        lower(0, 1),
    ]
    for tri in six:
        ax.add_patch(
            Polygon(
                tri,
                closed=True,
                facecolor=BLUE,
                alpha=0.25,
                edgecolor=BLUE,
                linewidth=1.2,
            )
        )
    ax.plot(*c, marker="o", markersize=8, color=ORANGE, zorder=5)
    ax.annotate(
        r"$(a_i,b_j)$",
        xy=c,
        xytext=(8, -10),
        textcoords="offset points",
        ha="left",
        va="top",
        fontsize=11,
    )
    # axis ticks named by grid coordinates
    ax.set_xticks(range(n))
    ax.set_xticklabels(["$a_{i-1}$", "$a_i$", "$a_{i+1}$", "$a_{i+2}$"])
    ax.set_yticks(range(m))
    ax.set_yticklabels(["$b_{j-1}$", "$b_j$", "$b_{j+1}$", "$b_{j+2}$"])
    ax.set_xlim(-0.2, n - 0.8)
    ax.set_ylim(-0.2, m - 0.8)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    return save_mpl(fig, "07-03-triangles", outdir)


FIGURES = {
    "07-03-interpolation": interpolation,
    "07-03-sos2": sos2,
    "07-03-triangles": triangles,
}
