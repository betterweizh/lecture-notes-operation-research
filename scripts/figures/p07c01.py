"""Figures of chapter 7.1, Logical Constraints and the Big-M Method.

* ``07-01-either-or``: the either-or set of Example ex-07-01-either-or
  (triangle T = {x1 + x2 <= 2} or square S = {x1 >= 3, x2 >= 3} in the box
  [0, 4]^2),
  its convex hull, and the projections of two big-M LP relaxations;
* ``07-01-fixed-charge``: a fixed-charge cost C(x) and the linear costs that
  its LP relaxations charge (M = u and M = 10 u).

Titles, legend entries and axis units are translated in the zh edition.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from .common import BLUE, GREEN, INK, MUTED, ORANGE, save_mpl

# ---------------------------------------------------------------------------
# Data of the either-or example (vertices listed counter-clockwise)
# ---------------------------------------------------------------------------
TRIANGLE = [(0, 0), (2, 0), (0, 2)]  # T: x1 + x2 <= 2, x >= 0
SQUARE = [(3, 3), (4, 3), (4, 4), (3, 4)]  # S: x1 >= 3, x2 >= 3, x <= 4
#: convex hull of T u S (area 10)
HULL = [(0, 0), (2, 0), (4, 3), (4, 4), (3, 4), (0, 2)]
#: projection of the big-M relaxation with M_T = 6, M_1 = M_2 = 3 (area 12)
BIGM_TIGHT = [(0, 0), (2, 0), (4, 2), (4, 4), (2, 4), (0, 2)]
#: with M = 10 the projection is the whole box (area 16)
BOX = [(0, 0), (4, 0), (4, 4), (0, 4)]
X_STAR = (3, 4)  # integer optimum of max -x1 + 2 x2, value 5
X_RELAX = (2, 4)  # optimum of the tight relaxation, value 6, delta = 1/3


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


def _setup_axes(ax):
    """Common axis limits, ticks and labels for both panels."""
    ax.set_xlim(-0.3, 4.6)
    ax.set_ylim(-0.3, 4.6)
    ax.set_aspect("equal")
    ax.set_xticks(range(5))
    ax.set_yticks(range(5))
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")


def either_or(outdir, lang="en"):
    """Union of a triangle and a square, its hull and two big-M relaxations.

    Used in Section sec-07-01-either-or (Example ex-07-01-either-or). Left panel:
    the feasible set T u S and its convex hull. Right panel: the projected
    LP relaxations with tight constants (hexagon) and with M = 10 (box).
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.2, 4.0))

    # --- left panel: the union T u S and its convex hull -----------------
    for poly in (TRIANGLE, SQUARE):
        ax1.add_patch(
            Polygon(
                poly,
                closed=True,
                facecolor=BLUE,
                alpha=0.55,
                edgecolor=BLUE,
                linewidth=1.5,
            )
        )
    ax1.add_patch(
        Polygon(
            HULL,
            closed=True,
            fill=False,
            edgecolor=INK,
            linestyle="--",
            linewidth=1.3,
        )
    )
    ax1.text(
        0.55, 0.55, "$T$", ha="center", va="center", color=INK, fontsize=13
    )
    ax1.text(3.5, 3.5, "$S$", ha="center", va="center", color=INK, fontsize=13)
    ax1.text(
        2.45,
        1.05,
        r"$\operatorname{conv}(T\cup S)$",
        ha="left",
        va="center",
        color=INK,
        fontsize=10.5,
    )
    # integer optimum of max -x1 + 2 x2
    ax1.plot(*X_STAR, marker="o", markersize=7, color=ORANGE, zorder=5)
    ax1.annotate(
        r"$\mathbf{x}^*=(3,4)$",
        xy=X_STAR,
        xytext=(-8, 8),
        textcoords="offset points",
        ha="right",
        va="bottom",
        fontsize=10.5,
    )
    ax1.set_title(
        _t(lang, "$T\\cup S$ and its convex hull", "$T\\cup S$ 及其凸包"),
        fontsize=11,
    )
    _setup_axes(ax1)

    # --- right panel: projections of the LP relaxations -----------------
    # M = 10: the whole box (light), tight constants: hexagon (darker)
    ax2.add_patch(
        Polygon(
            BOX,
            closed=True,
            facecolor=MUTED,
            alpha=0.12,
            edgecolor=MUTED,
            linestyle=":",
            linewidth=1.4,
        )
    )
    ax2.add_patch(
        Polygon(
            BIGM_TIGHT,
            closed=True,
            facecolor=GREEN,
            alpha=0.35,
            edgecolor=GREEN,
            linewidth=1.5,
        )
    )
    # the hull again for comparison (dashed) and the union (outline only)
    ax2.add_patch(
        Polygon(
            HULL,
            closed=True,
            fill=False,
            edgecolor=INK,
            linestyle="--",
            linewidth=1.1,
        )
    )
    for poly in (TRIANGLE, SQUARE):
        ax2.add_patch(
            Polygon(
                poly,
                closed=True,
                facecolor=BLUE,
                alpha=0.55,
                edgecolor=BLUE,
                linewidth=1.2,
            )
        )
    ax2.text(
        3.15,
        0.55,
        "$M=10$",
        ha="center",
        va="center",
        color=MUTED,
        fontsize=10.5,
    )
    ax2.text(
        2.2,
        2.1,
        "$M_T=6$\n$M_1=M_2=3$",
        ha="center",
        va="center",
        color=INK,
        fontsize=10,
    )
    # relaxed optimum (2, 4) with delta = 1/3 and the integer optimum
    ax2.plot(*X_RELAX, marker="s", markersize=6.5, color=GREEN, zorder=5)
    ax2.annotate(
        r"$(2,4),\ \delta=\frac{1}{3}$",
        xy=X_RELAX,
        xytext=(-6, -4),
        textcoords="offset points",
        ha="right",
        va="top",
        fontsize=10,
    )
    ax2.plot(*X_STAR, marker="o", markersize=7, color=ORANGE, zorder=5)
    ax2.set_title(_t(lang, "big-$M$ relaxations", "大 $M$ 松弛"), fontsize=11)
    _setup_axes(ax2)

    fig.tight_layout()
    return save_mpl(fig, "07-01-either-or", outdir)


def fixed_charge(outdir, lang="en"):
    """Fixed-charge cost and its LP-relaxation lines.

    Used in Section sec-07-01-fixed-charge. C(x) = 0 at x = 0 and
    F + c x on (0, u] with F = 500 $, c = 20 $/t, u = 100 t. The LP
    relaxation with x <= M delta charges (F / M + c) x: 25 x for M = u and
    20.5 x for M = 1000. Values at x = 60 t: 1700, 1500, 1230 $.
    """
    F, c, u = 500.0, 20.0, 100.0
    fig, ax = plt.subplots(figsize=(5.6, 4.0))

    # true cost: the point C(0) = 0 and the segment F + c x on (0, u]
    ax.plot([0, u], [F, F + c * u], color=BLUE, linewidth=2.2, label="$C(x)$")
    ax.plot(
        [0],
        [F],
        marker="o",
        markersize=6,
        markerfacecolor="white",
        markeredgecolor=BLUE,
        zorder=5,
    )  # open end of the segment at x = 0
    ax.plot(
        [0], [0], marker="o", markersize=6, color=BLUE, zorder=5
    )  # C(0) = 0

    # LP-relaxation costs (F / M + c) x for M = u and M = 1000
    ax.plot(
        [0, u],
        [0, (F / u + c) * u],
        color=ORANGE,
        linestyle="--",
        linewidth=1.8,
        label=_t(lang, "relaxation, $M=u$: $25x$", "松弛，$M=u$：$25x$"),
    )
    ax.plot(
        [0, u],
        [0, (F / 1000 + c) * u],
        color=GREEN,
        linestyle=":",
        linewidth=2.0,
        label=_t(
            lang,
            "relaxation, $M=1000$: $20.5x$",
            "松弛，$M=1000$：$20.5x$",
        ),
    )

    # values at x = 60 t
    x0 = 60
    # C(60) above-left of its point, 20.5 x below-right of its point (both
    # clear of the other lines); the middle one sits right of x0, in the
    # widening gap between 25x and 20.5x
    for val, col, off, ha, va in (
        (F + c * x0, BLUE, (-5, 3), "right", "bottom"),
        ((F / 1000 + c) * x0, GREEN, (5, -3), "left", "top"),
    ):
        ax.plot([x0], [val], marker="o", markersize=5, color=col, zorder=6)
        ax.annotate(
            f"{val:.0f}",
            xy=(x0, val),
            xytext=off,
            textcoords="offset points",
            ha=ha,
            va=va,
            fontsize=10,
            color=col,
        )
    val = (F / u + c) * x0
    ax.plot([x0], [val], marker="o", markersize=5, color=ORANGE, zorder=6)
    ax.annotate(
        f"{val:.0f}",
        xy=(x0, val),
        xytext=(64.5, 1550),
        textcoords="data",
        ha="left",
        va="center",
        fontsize=10,
        color=ORANGE,
    )
    ax.axvline(x0, color=MUTED, linewidth=0.8, linestyle=(0, (2, 3)), zorder=0)

    ax.set_xlim(-3, 108)
    ax.set_ylim(-80, 2700)
    ax.set_xlabel(_t(lang, "$x$ (t)", "$x$（吨）"))
    ax.set_ylabel(_t(lang, "cost ($)", "费用（美元）"))
    ax.legend(loc="upper left", frameon=False, fontsize=10)
    return save_mpl(fig, "07-01-fixed-charge", outdir)


FIGURES = {
    "07-01-either-or": either_or,
    "07-01-fixed-charge": fixed_charge,
}
