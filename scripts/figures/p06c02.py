"""Figures of chapter 6.2, Optimality Conditions.

* ``06-02-directions``: feasible and descent directions on the triangle
  x1 + x2 <= 2, x >= 0 for f(x) = (x1-1)^2 + (x2-2)^2, at the interior
  point (1, 0.5) and at the vertex (2, 0) (example ex-06-02-directions);
* ``06-02-saddle``: contours of f(x) = x1^2 + x1 x2 - x1 - x2
  = (x1 - 1)(x1 + x2) with the saddle point (1, -1) and the feasible
  quadrant x1 >= 1, x2 >= 0 (example ex-06-02-saddle).
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
from matplotlib.patches import Circle, Polygon

from .common import BLUE, GREEN, INK, MUTED, ORANGE, save_mpl


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


def _arrow(ax, start, vec, color=INK, lw=1.6, ls="-", scale=1.0):
    """Draw the vector ``vec`` (times ``scale``) as an arrow from ``start``."""
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
            mutation_scale=13,
        ),
        zorder=6,
    )
    return end


def directions(outdir, lang="en"):
    """Two panels: directions at (1, 0.5) (left) and at (2, 0) (right).

    Triangle Omega with vertices (0,0), (2,0), (0,2); contours of
    f = (x1-1)^2 + (x2-2)^2 (circles about (1, 2)) through the two points
    (levels 2.25 and 5) plus one inner level. Left: d1..d4 and the gradient
    (0, -3) (drawn at 1/3 length). Right: d5..d7 and the gradient (2, -4)
    (drawn at 1/4 length), with the tangent of the contour f = 5 at (2, 0).
    Used in ex-06-02-directions of 06-02-optimality-conditions.md.
    """
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4))
    center = (1.0, 2.0)  # unconstrained minimizer of f
    # (point, contour levels, directions {label: vector}, gradient, gradient scale)
    panels = [
        (
            (1.0, 0.5),
            [0.5625, 2.25],
            {
                r"$\mathbf{d}_1$": (-0.4, 1),
                r"$\mathbf{d}_2$": (0, 1),
                r"$\mathbf{d}_3$": (0.45, 1),
                r"$\mathbf{d}_4$": (-1, -0.8),
            },
            (0.0, -3.0),
            1 / 3,
        ),
        (
            (2.0, 0.0),
            [1.25, 5.0],
            {
                r"$\mathbf{d}_5$": (-1, 0.5),
                r"$\mathbf{d}_6$": (0.6, 1),
                r"$\mathbf{d}_7$": (1, 0.13),
            },
            (2.0, -4.0),
            1 / 4,
        ),
    ]
    # label offsets (in data units) of the direction tips
    tip_off = {
        r"$\mathbf{d}_1$": (-0.12, 0.08),
        r"$\mathbf{d}_2$": (0.0, 0.1),
        r"$\mathbf{d}_3$": (0.1, 0.06),
        r"$\mathbf{d}_4$": (-0.12, -0.06),
        r"$\mathbf{d}_5$": (-0.1, 0.1),
        r"$\mathbf{d}_6$": (0.08, 0.08),
        r"$\mathbf{d}_7$": (0.12, 0.0),
    }
    for k, (ax, (pt, levels, dirs, grad, gs)) in enumerate(zip(axes, panels)):
        # feasible triangle
        ax.add_patch(
            Polygon(
                [(0, 0), (2, 0), (0, 2)],
                closed=True,
                fc=to_rgba(BLUE, 0.13),  # light fill, opaque edge
                ec=BLUE,
                lw=1.3,
            )
        )
        ax.text(
            0.15,
            1.05 if k == 0 else 0.35,
            r"$\Omega$",
            color=BLUE,
            fontsize=13,
        )
        # contours of f: circles about (1, 2) with radius sqrt(level)
        for lev in levels:
            ax.add_patch(
                Circle(center, np.sqrt(lev), fill=False, ec=MUTED, lw=0.9)
            )
        ax.plot(*center, "o", ms=4, color=MUTED)
        ax.annotate(
            "$(1,2)$",
            center,
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=9,
            color=MUTED,
        )
        # contour labels f = level (placed on the right of each circle)
        lev_lab = levels[-1]
        r = np.sqrt(lev_lab)
        if k == 0:
            # label on the right of the circle through (1, 0.5)
            lab_xy = (
                center[0] + r * np.cos(0.35) + 0.04,
                center[1] + r * np.sin(0.35),
            )
        else:
            # the circle through (2, 0) leaves the window; label below its bottom
            lab_xy = (0.25, -0.5)
        ax.text(*lab_xy, f"$f={lev_lab:g}$", fontsize=9, color=MUTED)
        # directions
        for lab, v in dirs.items():
            col = GREEN if (lab in (r"$\mathbf{d}_5$",)) else INK
            end = _arrow(ax, pt, v, color=col, scale=0.8)
            dx, dy = tip_off[lab]
            ax.text(
                end[0] + dx,
                end[1] + dy,
                lab,
                ha="center",
                va="center",
                fontsize=11,
                color=col,
            )
        # gradient (scaled down to fit), orange
        gend = _arrow(ax, pt, grad, color=ORANGE, scale=gs)
        glab = r"$\nabla f=(0,-3)$" if k == 0 else r"$\nabla f=(2,-4)$"
        ax.text(
            gend[0] + 0.08,
            gend[1] - 0.05,
            glab,
            color=ORANGE,
            fontsize=10,
            ha="left",
            va="top",
        )
        if k == 1:
            # tangent of the contour f = 5 at (2, 0): normal (2, -4), direction (2, 1)
            t = np.array([-0.9, 0.6])
            ax.plot(2 + 2 * t, 0 + 1 * t, color=MUTED, lw=0.8, ls=":")
        ax.plot(*pt, "o", ms=6, color=INK, zorder=7)
        ptlab = "$(1,0.5)$" if k == 0 else "$(2,0)$"
        # left: below the edge x1 + x2 = 2; right: between d7 and the gradient
        ax.annotate(
            ptlab,
            pt,
            xytext=(6, -24) if k == 0 else (16, -15),
            textcoords="offset points",
            fontsize=10,
        )
        ax.set_xlim(-0.6, 3.1)
        ax.set_ylim(-1.3, 3.0)
        ax.set_aspect("equal")
        ax.axhline(0, color=MUTED, lw=0.7)
        ax.axvline(0, color=MUTED, lw=0.7)
        ax.set_xlabel("$x_1$")
        if k == 0:
            ax.set_ylabel("$x_2$")
        ax.set_title(
            (
                _t(lang, "interior point $(1,0.5)$", "内点 $(1,0.5)$")
                if k == 0
                else _t(lang, "vertex $(2,0)$", "顶点 $(2,0)$")
            ),
            fontsize=11,
        )
        ax.grid(False)
    fig.tight_layout()
    return save_mpl(fig, "06-02-directions", outdir)


def saddle(outdir, lang="en"):
    """Contours of f = (x1 - 1)(x1 + x2) around the saddle point (1, -1).

    Zero level set: lines x1 = 1 and x1 + x2 = 0 (dashed). Shaded: the
    feasible set x1 >= 1, x2 >= 0 of part (b); its edge x1 = 1, x2 >= 0
    (thick green) is a ray of global minima with f = 0. Used in
    ex-06-02-saddle of 06-02-optimality-conditions.md.
    """
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    xlim, ylim = (-1.5, 3.5), (-3.5, 2.5)
    g1 = np.linspace(*xlim, 400)
    g2 = np.linspace(*ylim, 400)
    X1, X2 = np.meshgrid(g1, g2)
    F = (X1 - 1) * (X2 + X1)
    # feasible quadrant x1 >= 1, x2 >= 0
    ax.add_patch(
        Polygon(
            [(1, 0), (xlim[1], 0), (xlim[1], ylim[1]), (1, ylim[1])],
            closed=True,
            color=BLUE,
            alpha=0.13,
            lw=0,
        )
    )
    # contours (positive solid blue, negative solid orange), zero dashed
    pos = [0.5, 1.5, 3, 5]
    neg = [-5, -3, -1.5, -0.5]
    ax.contour(X1, X2, F, levels=pos, colors=BLUE, linewidths=0.9)
    ax.contour(
        X1,
        X2,
        F,
        levels=neg,
        colors=ORANGE,
        linewidths=0.9,
        linestyles="solid",
    )
    ax.plot([1, 1], ylim, color=INK, lw=1.0, ls="--")
    ax.plot(g1, -g1, color=INK, lw=1.0, ls="--")
    # ray of minima x1 = 1, x2 >= 0
    ax.plot([1, 1], [0, ylim[1]], color=GREEN, lw=3.2, solid_capstyle="butt")
    # saddle point and the point (1, 0)
    ax.plot(1, -1, "o", ms=7, color=INK, zorder=5)
    ax.annotate(
        _t(lang, "saddle $(1,-1)$", "鞍点 $(1,-1)$"),
        (1, -1),
        xytext=(10, 4),
        textcoords="offset points",
        fontsize=10,
    )
    ax.plot(1, 0, "o", ms=7, color=GREEN, zorder=5)
    ax.annotate(
        "$(1,0)$",
        (1, 0),
        xytext=(6, 5),
        textcoords="offset points",
        fontsize=10,
        color=GREEN,
    )
    # sign labels of the four sectors
    ax.text(2.6, 1.6, "$f>0$", color=BLUE, fontsize=11)
    ax.text(-0.75, 2.1, "$f<0$", color=ORANGE, fontsize=11)
    ax.text(1.6, -3.25, "$f<0$", color=ORANGE, fontsize=11)
    # below the outermost contour f = 5, so that no level curve crosses it
    ax.text(-1.35, -2.7, "$f>0$", color=BLUE, fontsize=11)
    ax.text(1.85, 0.6, r"$\Omega$", color=BLUE, fontsize=13)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.grid(False)
    return save_mpl(fig, "06-02-saddle", outdir)


FIGURES = {
    "06-02-directions": directions,
    "06-02-saddle": saddle,
}
