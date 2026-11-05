"""Figures of chapter 2.2, Convex Functions.

Two matplotlib figures from closed-form functions (no randomness):

* ``02-02-chord-tangent``: chord above and tangent below the graph of a convex
  function of one variable (sections sec-02-02-definition, sec-02-02-first-order);
* ``02-02-epigraph-sublevel``: epigraph and a sublevel interval of a function
  of one variable, and sublevel sets of a convex quadratic in R^2
  (sections sec-02-02-epigraph, sec-02-02-sublevel).
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba

from .common import BLUE, INK, MUTED, ORANGE, save_mpl

FACE = to_rgba(BLUE, 0.18)  # fill of sets, edges stay opaque


def f1(x):
    """Convex function of one variable used in both figures."""
    return 0.5 * (x - 1.2) ** 2 + 0.5


def df1(x):
    """Derivative of ``f1``."""
    return x - 1.2


def _clean(ax, xlim, ylim):
    """1-D function panel: no ticks (the points are labeled instead)."""
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)


def chord_tangent(outdir, lang="en"):
    """Chord above the graph (a) and tangent below the graph (b).

    f(x) = (x - 1.2)^2 / 2 + 0.5 on [-0.6, 3.2].
    (a) x = 0.2, y = 2.8, theta = 0.4: the point theta x + (1 - theta) y =
    1.76; the chord value 0.4 f(x) + 0.6 f(y) = 1.468 lies above
    f(1.76) = 0.657. (b) tangent at x = 1.0; at y = 2.6 the graph value
    f(y) = 1.48 lies above the tangent value f(x) + f'(x)(y - x) = 0.20.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.6))
    xs = np.linspace(-0.6, 3.2, 300)
    xlim, ylim = (-0.7, 3.4), (-0.35, 2.6)

    # (a) chord between (x, f(x)) and (y, f(y))
    x, y, th = 0.2, 2.8, 0.4
    z = th * x + (1 - th) * y  # 1.76
    chord_z = th * f1(x) + (1 - th) * f1(y)  # 1.468
    ax1.plot(xs, f1(xs), color=INK, lw=2.0)
    ax1.plot([x, y], [f1(x), f1(y)], color=ORANGE, lw=1.8)
    for p, lab in ((x, "$x$"), (z, r"$\theta x+(1-\theta)y$"), (y, "$y$")):
        top = chord_z if p == z else f1(p)
        ax1.plot([p, p], [ylim[0], top], color=MUTED, lw=0.9, ls="--")
        ax1.text(p, ylim[0] - 0.08, lab, ha="center", va="top")
    ax1.plot([x, y], [f1(x), f1(y)], "o", color=ORANGE, ms=5, zorder=5)
    ax1.plot([z, z], [f1(z), chord_z], "o", color=INK, ms=5, zorder=5)
    ax1.annotate(
        r"$\theta f(x)+(1-\theta)f(y)$", xy=(z, chord_z), xytext=(z - 0.25, chord_z + 0.75),
        ha="right", arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9),
    )  # fmt: skip
    ax1.annotate(
        r"$f(\theta x+(1-\theta)y)$", xy=(z, f1(z)), xytext=(z + 0.3, f1(z) - 0.45),
        ha="left", arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9),
    )  # fmt: skip
    ax1.text(3.05, f1(3.05) + 0.05, "$f$", ha="right", va="bottom")
    _clean(ax1, xlim, ylim)
    ax1.set_title("(a)")

    # (b) tangent at x
    x, y = 1.0, 2.6
    tan = lambda s: f1(x) + df1(x) * (s - x)  # noqa: E731
    ax2.plot(xs, f1(xs), color=INK, lw=2.0)
    ts = np.linspace(-0.6, 3.2, 2)
    ax2.plot(ts, tan(ts), color=ORANGE, lw=1.8)
    ax2.plot(x, f1(x), "o", color=ORANGE, ms=5, zorder=5)
    for p, lab in ((x, "$x$"), (y, "$y$")):
        ax2.plot([p, p], [ylim[0], f1(p)], color=MUTED, lw=0.9, ls="--")
        ax2.text(p, ylim[0] - 0.08, lab, ha="center", va="top")
    ax2.plot([y, y], [tan(y), f1(y)], "o", color=INK, ms=5, zorder=5)
    ax2.text(y + 0.1, f1(y), "$f(y)$", va="center")
    ax2.text(y + 0.1, tan(y) + 0.02, "$f(x)+f′(x)(y-x)$", va="center")
    ax2.text(x - 0.05, f1(x) + 0.12, "$f(x)$", ha="right")
    _clean(ax2, xlim, ylim)
    ax2.set_title("(b)")
    return save_mpl(fig, "02-02-chord-tangent", outdir)


def fq(x1, x2):
    """Convex quadratic of panel (b): Hessian [[2, -1/2], [-1/2, 4]] > 0."""
    return (x1 - 1) ** 2 + 2 * (x2 - 1) ** 2 - 0.5 * x1 * x2


def epigraph_sublevel(outdir, lang="en"):
    """Epigraph and sublevel interval (a); sublevel sets in R^2 (b).

    (a) f(x) = (x - 1.2)^2 / 2 + 0.5; the epigraph is shaded up to the top
    of the window; at level gamma = 1.4 the sublevel set S_gamma is the
    interval [1.2 - sqrt(1.8), 1.2 + sqrt(1.8)] = [-0.142, 2.542], drawn
    on the x axis. (b) contours of fq at levels -0.25, 0.25, 1 (bold), 2, 3.5
    (minimum -0.613 at (1.290, 1.161)); the sublevel set S_1 is shaded.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.8))

    # (a) epigraph: region between the graph and the top of the window
    xs = np.linspace(-0.6, 3.2, 300)
    top = 2.6
    ax1.fill_between(xs, f1(xs), top, color=FACE, lw=0)
    ax1.plot(xs, f1(xs), color=INK, lw=2.0)
    ax1.text(1.2, 1.95, r"$\operatorname{epi} f$", ha="center", fontsize=13)
    gamma = 1.4
    lo, hi = 1.2 - np.sqrt(2 * (gamma - 0.5)), 1.2 + np.sqrt(2 * (gamma - 0.5))
    ax1.axhline(gamma, color=ORANGE, lw=1.4, ls="--")
    # level label above the dashed line at its right end (clear of the dashes)
    ax1.text(
        3.35, gamma + 0.04, r"$\gamma$", color=ORANGE, ha="right", va="bottom"
    )
    for p in (lo, hi):  # drop lines from the crossing points to the axis
        ax1.plot([p, p], [-0.3, gamma], color=MUTED, lw=0.9, ls="--")
    ax1.plot([lo, hi], [-0.3, -0.3], color=ORANGE, lw=4, solid_capstyle="butt")
    ax1.text(
        (lo + hi) / 2,
        -0.15,
        r"$S_\gamma$",
        color=ORANGE,
        ha="center",
        va="bottom",
    )
    ax1.text(3.05, f1(3.05) - 0.08, "$f$", ha="left", va="top")
    _clean(ax1, (-0.7, 3.4), (-0.45, 2.6))
    ax1.set_xlabel("$x$")
    ax1.set_title("(a)")

    # (b) contour lines and the shaded sublevel set S_1
    g1, g2 = np.meshgrid(
        np.linspace(-0.5, 3.0, 400), np.linspace(-0.3, 2.6, 400)
    )
    F = fq(g1, g2)
    ax2.contourf(g1, g2, F, levels=[F.min() - 1, 1.0], colors=[FACE])
    levels = [-0.25, 0.25, 2.0, 3.5]  # thin contours (all solid lines)
    cs = ax2.contour(g1, g2, F, levels=levels, colors=[MUTED], linewidths=1.0,
                     linestyles="solid")  # fmt: skip
    ax2.clabel(
        cs,
        fmt={v: f"{v:g}".replace("-", "\u2212") for v in levels},
        fontsize=9,
        inline_spacing=3,
    )
    # the boundary f = 1 of the shaded sublevel set, bold
    ax2.contour(g1, g2, F, levels=[1.0], colors=[INK], linewidths=1.8)
    # minimiser of fq (solves the 2 x 2 system grad fq = 0)
    xm = np.linalg.solve([[2, -0.5], [-0.5, 4]], [2, 4])
    ax2.plot(*xm, "o", color=INK, ms=4)
    ax2.text(0.45, 0.62, "$S_1$", ha="center", va="center", fontsize=13)
    ax2.set_xlim(-0.5, 3.0)
    ax2.set_ylim(-0.3, 2.6)
    ax2.set_aspect("equal")
    ax2.set_xlabel("$x_1$")
    ax2.set_ylabel("$x_2$")
    ax2.set_title("(b)")
    return save_mpl(fig, "02-02-epigraph-sublevel", outdir)


FIGURES = {
    "02-02-chord-tangent": chord_tangent,
    "02-02-epigraph-sublevel": epigraph_sublevel,
}
