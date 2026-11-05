"""Figures of chapter 6.3, Karush–Kuhn–Tucker Conditions.

* ``06-03-kkt-geometry``: lens-shaped feasible set bounded by two curved
  constraints active at x* = (0, 0), an inactive linear constraint, contours
  of f(x) = (x1 - 1.5)^2 + 2 x2^2, and -grad f(x*) inside the cone of the
  active constraint gradients (section "Geometry");
* ``06-03-cubic-cq``: two panels on constraint qualifications: the cubic
  description (x1 + x2 - 1)^3 <= 0 of a triangle with the infeasible
  direction d = (2, 1) at (0.5, 0.5) (example ex-06-03-cubic), and the unit
  disc at (1, 0) with the tangent direction (0, 1).
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
from matplotlib.patches import Circle, Polygon

from .common import BLUE, GREEN, INK, MUTED, ORANGE, save_mpl


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


def _arrow(ax, start, vec, color=INK, lw=1.7, scale=1.0, ls="-"):
    """Arrow from ``start`` along ``scale * vec`` (data coordinates)."""
    end = (start[0] + scale * vec[0], start[1] + scale * vec[1])
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            lw=lw,
            linestyle=ls,
            shrinkA=0,
            shrinkB=0,
            mutation_scale=14,
        ),
        zorder=6,
    )
    return end


def kkt_geometry(outdir, lang="en"):
    """KKT geometry at x* = (0, 0).

    Constraints: g1 = x2 + 0.7 x1 + 0.3 x1^2 <= 0 (upper curve),
    g2 = -x2 + 0.7 x1 + 0.3 x1^2 <= 0 (lower curve), both active at 0;
    g3 = -(x1 + 1.5 x2 + 1.6) <= 0 (line, inactive). Objective
    f = (x1 - 1.5)^2 + 2 x2^2 with grad f(0) = (-3, 0). Multipliers
    lambda1 = lambda2 = 15/7. The constraint gradients are drawn at 0.8
    times, -grad f at 1/3 of their length.
    """
    fig, ax = plt.subplots(figsize=(6.0, 3.9))
    xlim, ylim = (-2.7, 2.6), (-1.35, 1.35)
    # feasible lens: |x2| <= -(0.7 x1 + 0.3 x1^2), x1 in [-7/3, 0], and g3 <= 0
    t = np.linspace(-7 / 3, 0, 400)
    up = -(0.7 * t + 0.3 * t**2)  # upper boundary g1 = 0
    lo = -up  # lower boundary g2 = 0
    # clip by g3: x2 >= -(x1 + 1.6)/1.5
    lo_c = np.maximum(lo, -(t + 1.6) / 1.5)
    keep = lo_c <= up
    poly = np.r_[
        np.c_[t[keep], up[keep]], np.c_[t[keep][::-1], lo_c[keep][::-1]]
    ]
    ax.add_patch(Polygon(poly, closed=True, color=BLUE, alpha=0.15, lw=0))
    # boundary curves g1 = 0, g2 = 0 (extended past x*), and line g3 = 0
    s = np.linspace(-2.5, 0.85, 400)
    ax.plot(s, -(0.7 * s + 0.3 * s**2), color=BLUE, lw=1.4)
    ax.plot(s, 0.7 * s + 0.3 * s**2, color=BLUE, lw=1.4)
    ax.plot(s, -(s + 1.6) / 1.5, color=BLUE, lw=1.0, ls="--")
    # curve labels at the right ends: g1 = 0 descends, g2 = 0 rises past x*
    ax.text(0.9, -0.8, "$g_1=0$", color=BLUE, fontsize=10, va="center")
    ax.text(0.9, 0.8, "$g_2=0$", color=BLUE, fontsize=10, va="center")
    # label above the upper end of the line g3 = 0, clear of the dashes
    ax.text(-2.3, 0.62, "$g_3=0$", color=BLUE, fontsize=10)
    ax.text(-1.05, 0.05, r"$\Omega$", color=BLUE, fontsize=14)
    # contours of f (ellipses about (1.5, 0)); the level 2.25 passes x*
    g1 = np.linspace(*xlim, 300)
    g2 = np.linspace(*ylim, 300)
    X1, X2 = np.meshgrid(g1, g2)
    F = (X1 - 1.5) ** 2 + 2 * X2**2
    ax.contour(
        X1,
        X2,
        F,
        levels=[0.25, 1.0, 2.25, 4.0],
        colors=MUTED,
        linewidths=0.8,
        linestyles=":",
    )
    ax.plot(1.5, 0, "+", ms=8, color=MUTED)
    # gradients at x* (half length): grad g1 = (0.7, 1), grad g2 = (0.7, -1)
    x0 = (0.0, 0.0)
    e1 = _arrow(ax, x0, (0.7, 1.0), color=GREEN, scale=0.8)
    e2 = _arrow(ax, x0, (0.7, -1.0), color=GREEN, scale=0.8)
    # gradient labels to the left of the arrow heads (clear of the curves)
    ax.text(
        e1[0] - 0.12,
        e1[1] + 0.1,
        r"$\nabla g_1(\mathbf{x}^*)$",
        color=GREEN,
        fontsize=10,
        ha="right",
    )
    ax.text(
        e2[0] - 0.12,
        e2[1] - 0.1,
        r"$\nabla g_2(\mathbf{x}^*)$",
        color=GREEN,
        fontsize=10,
        ha="right",
    )
    # -grad f(x*) = (3, 0) (drawn at 1/3 length) inside the cone
    e3 = _arrow(ax, x0, (3.0, 0.0), color=ORANGE, scale=0.33, lw=2.2)
    ax.text(
        e3[0] + 0.06,
        e3[1] + 0.12,  # clear of the center mark + at (1.5, 0)
        r"$-\nabla f(\mathbf{x}^*)$",
        color=ORANGE,
        fontsize=10,
    )
    ax.plot(*x0, "o", ms=6, color=INK, zorder=7)
    ax.annotate(
        r"$\mathbf{x}^*$",
        x0,
        xytext=(-18, 6),
        textcoords="offset points",
        fontsize=11,
    )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.grid(False)
    return save_mpl(fig, "06-03-kkt-geometry", outdir)


def cubic_cq(outdir, lang="en"):
    """Left: cubic description of the triangle; right: unit disc.

    Left: Omega = triangle (0,0), (1,0), (0,1); x* = (0.5, 0.5) on the edge
    x1 + x2 = 1 where grad g1 = 3 (x1 + x2 - 1)^2 (1, 1) = 0; d = (2, 1)
    (drawn at 0.35 length) leaves Omega since g1(x* + a d) = 27 a^3.
    Right: disc x1^2 + x2^2 <= 1 at (1, 0): grad g = (2, 0) (half length),
    d = (0, 1) in L(x*) leaves the disc; the circle (a feasible curve) has
    velocity (0, 1).
    """
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0))
    # --- left panel ---------------------------------------------------------
    ax = axes[0]
    ax.add_patch(
        Polygon(
            [(0, 0), (1, 0), (0, 1)],
            closed=True,
            fc=to_rgba(BLUE, 0.15),  # light fill, opaque edge
            ec=BLUE,
            lw=1.3,
        )
    )
    ax.text(0.2, 0.25, r"$\Omega$", color=BLUE, fontsize=13)
    # edge label inside the triangle, parallel to the edge
    ax.text(
        0.6,
        0.22,
        r"$(x_1+x_2-1)^3=0$",
        color=BLUE,
        fontsize=9,
        rotation=-45,
        ha="center",
        va="center",
    )
    xs = (0.5, 0.5)
    end = _arrow(ax, xs, (2, 1), color=ORANGE, scale=0.35)
    ax.text(
        end[0] + 0.03,
        end[1] + 0.03,
        r"$\mathbf{d}=(2,1)$",
        color=ORANGE,
        fontsize=10,
    )
    ax.plot(*xs, "o", ms=6, color=INK, zorder=7)
    ax.annotate(
        r"$\mathbf{x}^*=(0.5,0.5)$",
        xs,
        # right of x*, between d and the edge (clear of the opaque edge)
        xytext=(14, -2),
        textcoords="offset points",
        ha="left",
        va="top",
        fontsize=10,
    )
    ax.text(
        0.55,
        -0.27,
        r"$\nabla g_1(\mathbf{x}^*)=\mathbf{0}$",
        fontsize=10,
        color=INK,
        ha="center",
    )
    ax.set_xlim(-0.25, 1.45)
    ax.set_ylim(-0.4, 1.2)
    ax.set_title(
        _t(
            lang,
            "cubic constraint: $L=\\mathbb{R}^2$",
            "三次约束：$L=\\mathbb{R}^2$",
        ),
        fontsize=11,
    )
    # --- right panel --------------------------------------------------------
    ax = axes[1]
    # light fill, opaque boundary circle
    ax.add_patch(Circle((0, 0), 1, fc=to_rgba(BLUE, 0.15), ec=BLUE, lw=1.3))
    ax.text(-0.3, -0.15, r"$\Omega$", color=BLUE, fontsize=13)
    xs = (1.0, 0.0)
    eg = _arrow(ax, xs, (2, 0), color=GREEN, scale=0.25)
    ax.text(
        eg[0] + 0.07,
        eg[1] + 0.12,
        r"$\nabla g=(2,0)$",
        color=GREEN,
        fontsize=10,
        va="center",
    )
    ed = _arrow(ax, xs, (0, 1), color=ORANGE, scale=0.55)
    ax.text(
        ed[0] + 0.05,
        ed[1],
        r"$\mathbf{d}=(0,1)$",
        color=ORANGE,
        fontsize=10,
        va="center",
    )
    # tangent line x1 = 1 (dotted)
    ax.plot([1, 1], [-0.75, 0.75], color=MUTED, lw=0.8, ls=":")
    ax.plot(*xs, "o", ms=6, color=INK, zorder=7)
    ax.annotate(
        "$(1,0)$",
        xs,
        xytext=(-34, -4),
        textcoords="offset points",
        fontsize=10,
    )
    ax.set_xlim(-1.25, 1.85)
    ax.set_ylim(-1.25, 1.25)
    ax.set_title(
        _t(
            lang,
            "disc: $\\mathbf{d}\\in L$, but not feasible",
            "圆盘：$\\mathbf{d}\\in L$，但不是可行方向",
        ),
        fontsize=11,
    )
    for ax in axes:
        ax.set_aspect("equal")
        ax.axhline(0, color=MUTED, lw=0.6)
        ax.axvline(0, color=MUTED, lw=0.6)
        ax.set_xlabel("$x_1$")
        ax.grid(False)
    axes[0].set_ylabel("$x_2$")
    fig.tight_layout()
    return save_mpl(fig, "06-03-cubic-cq", outdir)


FIGURES = {
    "06-03-kkt-geometry": kkt_geometry,
    "06-03-cubic-cq": cubic_cq,
}
