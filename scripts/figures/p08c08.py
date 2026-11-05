"""Figures of chapter 8.8, Optimality Conditions and KKT (Part VIII, Examples
and Problems).

* ``08-08-difficulties``: (left) f(x) = 1 - 1/x on x >= 0.5 with its
  asymptote 1; (right) the disconnected set {x > 0 : g1(x) >= g2(x)} for
  g1 = (x - 2)^2 / x and g2 = 2.5 - 0.3 x (example ex-08-08-difficulties);
* ``08-08-descent-cone``: the triangle x1 + x2 <= 5, x1 - x2 <= 1, x1 >= 0
  with the cone of feasible descent directions of
  f = (x1 - 6)^2 + (x2 - 7)^2 at the vertex (3, 2)
  (example ex-08-08-descent-cone).
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

from .common import BLUE, GREEN, INK, MUTED, ORANGE, save_mpl


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


def difficulties(outdir, lang="en"):
    """Unattained supremum (left) and disconnected feasible set (right)."""
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6))
    # --- left: f = 1 - 1/x on [0.5, inf) ------------------------------------
    ax = axes[0]
    x = np.linspace(0.5, 8, 400)
    ax.plot(x, 1 - 1 / x, color=BLUE, lw=1.8)
    ax.axhline(1, color=MUTED, ls="--", lw=1.0)  # asymptote, never reached
    sup = _t(lang, "sup $=1$", "上确界 $=1$")
    ax.text(6.2, 1.07, sup, color=MUTED, fontsize=10)
    ax.plot(0.5, -1, "o", ms=6, color=ORANGE)  # minimum at the end point
    ax.annotate(
        _t(lang, "min at $x=0.5$", "极小点 $x=0.5$"),
        (0.5, -1),
        xytext=(10, -2),
        textcoords="offset points",
        fontsize=9,
    )
    ax.axvspan(0, 0.5, color=MUTED, alpha=0.12, lw=0)  # infeasible part
    ax.set_xlim(0, 8)
    ax.set_ylim(-1.25, 1.35)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$f(x)=1-1/x$")
    # --- right: g1 >= g2 ----------------------------------------------------
    ax = axes[1]
    x = np.linspace(0.3, 6.5, 500)
    g1 = (x - 2) ** 2 / x
    g2 = 2.5 - 0.3 * x
    ax.plot(x, g1, color=BLUE, lw=1.8, label="$g_1(x)=(x-2)^2/x$")
    ax.plot(x, g2, color=ORANGE, lw=1.8, label="$g_2(x)=2.5-0.3x$")
    # roots of 1.3 x^2 - 6.5 x + 4 = 0 (crossings)
    r1 = (6.5 - np.sqrt(21.45)) / 2.6
    r2 = (6.5 + np.sqrt(21.45)) / 2.6
    # dotted drop lines from the crossings to the x-axis (not through the
    # legend above or the root labels below)
    for r in (r1, r2):
        ax.plot([r, r], [0, 2.5 - 0.3 * r], color=MUTED, ls=":", lw=0.9)
    # the two pieces of Omega on the x-axis (thick green); Omega = (0, r1]
    # and [r2, inf)
    ax.plot([0, r1], [0, 0], color=GREEN, lw=5, solid_capstyle="butt")
    ax.plot([r2, 6.5], [0, 0], color=GREEN, lw=5, solid_capstyle="butt")
    ax.text(
        (r2 + 6.5) / 2,
        0.22,
        r"$\Omega$",
        color=GREEN,
        ha="center",
        fontsize=12,
    )
    ax.text(0.5, 0.22, r"$\Omega$", color=GREEN, ha="center", fontsize=12)
    ax.text(r1, -0.45, f"{r1:.3f}", ha="center", fontsize=9, color=MUTED)
    ax.text(r2, -0.45, f"{r2:.3f}", ha="center", fontsize=9, color=MUTED)
    ax.set_xlim(0, 6.6)
    ax.set_ylim(-0.7, 5)
    ax.set_xlabel("$x$")
    ax.legend(loc="upper center", frameon=False, fontsize=9)
    for a in axes:
        a.grid(False)
    fig.tight_layout()
    return save_mpl(fig, "08-08-difficulties", outdir)


def descent_cone(outdir, lang="en"):
    """Feasible descent cone at the vertex (3, 2).

    Feasible directions: d1 <= -|d2| (between the edges towards (0, 5) and
    (0, -1)). Descent: -6 d1 - 10 d2 < 0. Intersection: directions (-1, s),
    0.6 < s <= 1, shaded. Also the gradient (-6, -10) (scaled), the line
    orthogonal to it, the interior point (1, 1) and the optimum (2, 3).
    """
    fig, ax = plt.subplots(figsize=(5.4, 4.8))
    # feasible triangle (0,5), (3,2), (0,-1)
    ax.add_patch(
        Polygon(
            [(0, 5), (3, 2), (0, -1)],
            closed=True,
            color=BLUE,
            alpha=0.12,
            ec=BLUE,
            lw=1.3,
        )
    )
    ax.text(0.5, 2.0, r"$\Omega$", color=BLUE, fontsize=14)
    v = np.array([3.0, 2.0])
    # outline of the triangle (the fill alone is faint)
    ax.plot([0, 3, 0, 0], [5, 2, -1, 5], color=BLUE, lw=1.3)
    # cone of feasible descent directions: between (-1, 1) and (-1, 0.6)
    L = 2.9
    tip1 = v + L * np.array([-1, 1]) / np.sqrt(2)
    tip2 = v + L * np.array([-1, 0.6]) / np.hypot(1, 0.6)
    ax.add_patch(
        Polygon([v, tip1, tip2], closed=True, color=GREEN, alpha=0.3, lw=0)
    )
    # the three directions of the solution table, at different lengths so
    # that their labels (left of the tips) do not collide; the label of
    # (-3, 2) sits below its tip, clear of the dotted line
    for d, lab, length, dy, va in [
        ((-3, 3), "$(-3,3)$", 2.6, 0.02, "center"),
        ((-5, 4), "$(-5,4)$", 2.0, 0.02, "center"),
        ((-3, 2), "$(-3,2)$", 1.4, -0.08, "top"),
    ]:
        e = v + length * np.asarray(d, float) / np.hypot(*d)
        ax.annotate(
            "",
            xy=e,
            xytext=v,
            arrowprops=dict(
                arrowstyle="-|>",
                color=GREEN,
                lw=1.4,
                shrinkA=0,
                shrinkB=0,
                mutation_scale=12,
            ),
        )
        ax.text(
            e[0] - 0.1,
            e[1] + dy,
            lab,
            fontsize=9,
            color=GREEN,
            ha="right",
            va=va,
        )
    # line through (3, 2) orthogonal to the gradient: direction (-1, 0.6) / (1, -0.6)
    t = np.array([-2.6, 1.6])
    ax.plot(3 - t, 2 + 0.6 * t, color=MUTED, ls=":", lw=1.0)
    # gradient (-6, -10) at (3, 2), drawn at 1/6 length, and its negative
    g = np.array([-6.0, -10.0]) / 6
    ax.annotate(
        "",
        xy=v + g,
        xytext=v,
        arrowprops=dict(
            arrowstyle="-|>",
            color=ORANGE,
            lw=1.6,
            shrinkA=0,
            shrinkB=0,
            mutation_scale=13,
        ),
    )
    # label below the x1-axis, which would otherwise cross it
    ax.text(
        v[0] + g[0] + 0.1,
        -0.1,
        r"$\nabla f(3,2)$",
        color=ORANGE,
        fontsize=10,
        va="top",
    )
    ax.annotate(
        "",
        xy=v - g,
        xytext=v,
        arrowprops=dict(
            arrowstyle="-|>",
            color=MUTED,
            lw=1.2,
            ls="--",
            shrinkA=0,
            shrinkB=0,
            mutation_scale=12,
        ),
    )
    ax.text(
        v[0] - g[0] + 0.08,
        v[1] - g[1],
        r"$-\nabla f$",
        color=MUTED,
        fontsize=10,
    )
    # points: vertex, interior point, optimum
    ax.plot(*v, "o", ms=6, color=INK, zorder=6)
    ax.annotate(
        "$(3,2)$", v, xytext=(8, -4), textcoords="offset points", fontsize=10
    )
    ax.plot(1, 1, "o", ms=5, color=INK)
    ax.annotate(
        "$(1,1)$",
        (1, 1),
        xytext=(6, 2),  # above the edge x1 - x2 = 1
        textcoords="offset points",
        fontsize=10,
    )
    ax.plot(2, 3, "o", ms=6, color=ORANGE, zorder=6)
    ax.annotate(
        r"$\mathbf{x}^*=(2,3)$",
        (2, 3),
        xytext=(8, 6),
        textcoords="offset points",
        fontsize=10,
        color=ORANGE,
    )
    ax.set_xlim(-0.6, 5.2)
    ax.set_ylim(-1.5, 5.6)
    ax.set_aspect("equal")
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.axvline(0, color=MUTED, lw=0.6)
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.grid(False)
    return save_mpl(fig, "08-08-descent-cone", outdir)


FIGURES = {
    "08-08-difficulties": difficulties,
    "08-08-descent-cone": descent_cone,
}
