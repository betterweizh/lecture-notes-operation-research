"""Figures of chapter 1.2, Vector and Matrix Calculus.

Three matplotlib figures, all from closed-form functions:

* ``01-02-linear-level-lines``: level lines of f(x) = x1 + x2 with the gradient
  and the negative gradient at (2, 1) (example ex-01-02-linear-level-lines);
* ``01-02-paraboloid-contours``: graph and contours of
  f(x) = (x1 - 2)^2 + (x2 - 2)^2 with the gradient at (1.4, 1.4) and the
  tangent of the level circle (example ex-01-02-paraboloid-contours);
* ``01-02-taylor-models``: f(x) = e^x - 2x with its first- and second-order
  Taylor models at 0 (example ex-01-02-taylor-models).
"""

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np

from .common import BLUE, GREEN, INK, MUTED, ORANGE, save_mpl


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


def _arrow(ax, start, end, color=INK, lw=1.8, ls="-", zorder=5):
    """Draw a vector as an arrow from ``start`` to ``end`` (data coords)."""
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
        zorder=zorder,
    )


def linear_level_lines(outdir, lang="en"):
    """Level lines f = c of f(x) = x1 + x2 and the gradient at (2, 1).

    The gradient (1, 1) is the same at every point; it is drawn at the
    point x = (2, 1) (on the line f = 3), together with -grad f, and is
    perpendicular to every level line. Used in example
    ex-01-02-linear-level-lines of 01-02-matrix-calculus.md.
    """
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    lim = (-1.5, 4.0)
    g = np.array(lim)
    # level lines x2 = c - x1 for c = -1, ..., 5 (labeled c = 0..3 as in
    # the example, the others muted)
    for c in range(-1, 6):
        strong = 0 <= c <= 3
        ax.plot(
            g,
            c - g,
            color=BLUE if strong else MUTED,
            lw=1.5 if strong else 0.9,
        )
        if strong:
            # label near the left/top end of each line, inside the window
            x_lab = max(lim[0] + 0.1, c - lim[1] + 0.35)
            ax.text(
                x_lab,
                c - x_lab + 0.12,
                f"$f={c}$",
                color=BLUE,
                fontsize=11,
                rotation=-45,
                rotation_mode="anchor",
            )
    x = np.array([2.0, 1.0])
    grad = np.array([1.0, 1.0])
    _arrow(ax, x, x + grad, color=ORANGE, lw=2.2)
    _arrow(ax, x, x - grad, color=GREEN, lw=2.2)
    ax.plot(*x, "o", color=INK, ms=6, zorder=6)
    ax.text(*(x + [0.12, -0.3]), r"$(2,1)$", fontsize=11)
    ax.text(
        *(x + grad + [-0.7, 0.18]),
        r"$\nabla f(2,1)$",
        color=ORANGE,
        fontsize=12,
    )
    ax.text(
        *(x - grad + [-1.05, -0.35]),
        r"$-\nabla f(2,1)$",
        color=GREEN,
        fontsize=12,
    )
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect("equal")
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    return save_mpl(fig, "01-02-linear-level-lines", outdir)


def paraboloid_contours(outdir, lang="en"):
    """Graph and contours of f(x) = (x1 - 2)^2 + (x2 - 2)^2.

    (a) the graph z = f(x) over [1, 3]^2 with the level f = 0.72 drawn on the
        surface and its shadow on the floor;
    (b) contours f = 0.1, 0.25, 0.5, 0.72, 1, 1.5 (circles about (2, 2)); at
        x = (1.4, 1.4), on the level f = 0.72, the gradient (-1.2, -1.2)
        (drawn at 1/3 of its length) points away from the minimiser and is
        perpendicular to the dotted tangent line, direction (1, -1).
    Used in example ex-01-02-paraboloid-contours of 01-02-matrix-calculus.md.
    """
    fig = plt.figure(figsize=(9.6, 4.3))
    f = lambda a, b: (a - 2) ** 2 + (b - 2) ** 2
    x0 = np.array([1.4, 1.4])
    level = f(*x0)  # 0.72
    r = np.sqrt(level)  # radius 0.8485 of the level circle
    s = np.linspace(0, 2 * np.pi, 200)

    # --- (a) surface ---------------------------------------------------
    ax3 = fig.add_subplot(1, 2, 1, projection="3d")
    g = np.linspace(1, 3, 25)
    X1, X2 = np.meshgrid(g, g)
    ax3.plot_surface(
        X1,
        X2,
        f(X1, X2),
        color=BLUE,
        alpha=0.25,
        edgecolor=BLUE,
        linewidth=0.3,
        rstride=1,
        cstride=1,
        shade=False,
    )
    # the level f = 0.72 on the surface and on the floor z = 0
    cx, cy = 2 + r * np.cos(s), 2 + r * np.sin(s)
    ax3.plot(cx, cy, np.full_like(s, level), color=ORANGE, lw=1.8)
    ax3.plot(cx, cy, np.zeros_like(s), color=ORANGE, lw=1.2, ls="--")
    ax3.plot([x0[0]], [x0[1]], [level], "o", color=INK, ms=5)
    ax3.plot([x0[0]], [x0[1]], [0], "o", color=INK, ms=4)
    ax3.set_xlabel("$x_1$")
    ax3.set_ylabel("$x_2$")
    ax3.set_zlabel("$f$")
    ax3.set_xticks([1, 2, 3])
    ax3.set_yticks([1, 2, 3])
    ax3.set_zticks([0, 1, 2])
    ax3.view_init(elev=24, azim=-58)
    for axis in (ax3.xaxis, ax3.yaxis, ax3.zaxis):
        axis.pane.fill = False
        axis.pane.set_edgecolor(MUTED)
    ax3.set_title("(a)", loc="left")

    # --- (b) contours --------------------------------------------------
    ax = fig.add_subplot(1, 2, 2)
    for c in (0.1, 0.25, 0.5, level, 1.0, 1.5):
        rad = np.sqrt(c)
        hl = np.isclose(c, level)
        ax.plot(
            2 + rad * np.cos(s),
            2 + rad * np.sin(s),
            color=ORANGE if hl else BLUE,
            lw=1.8 if hl else 1.0,
        )
    ax.plot(2, 2, "+", color=INK, ms=9)
    # tangent line through x0, direction (1, -1)
    t = np.array([-0.55, 0.55])
    ax.plot(x0[0] + t, x0[1] - t, color=INK, ls=":", lw=1.5)
    grad = np.array([2 * (x0[0] - 2), 2 * (x0[1] - 2)])  # (-1.2, -1.2)
    _arrow(ax, x0, x0 + grad / 3, color=ORANGE, lw=2.2)
    _arrow(ax, x0, x0 - grad / 3, color=GREEN, lw=2.2)
    ax.plot(*x0, "o", color=INK, ms=5, zorder=6)
    # white halo keeps labels legible where they cross contour lines
    halo = [pe.withStroke(linewidth=3, foreground="white")]
    ax.text(
        *(x0 + grad / 3 + [-0.02, -0.16]),
        r"$\nabla f(\bar{\mathbf{x}})$",
        color=ORANGE,
        fontsize=12,
        path_effects=halo,
    )
    ax.text(
        *(x0 - grad / 3 + [0.04, -0.13]),
        r"$-\nabla f(\bar{\mathbf{x}})$",
        color=GREEN,
        fontsize=12,
        path_effects=halo,
    )
    ax.annotate(
        r"$\bar{\mathbf{x}}=(1.4,1.4)$",
        xy=x0,
        xytext=(0.78, 2.5),
        fontsize=11,
        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8, shrinkB=4),
        path_effects=halo,
    )
    ax.text(
        2.62, 2.62, "$f=0.72$", color=ORANGE, fontsize=11, path_effects=halo
    )
    ax.set_xlim(0.7, 3.2)
    ax.set_ylim(0.7, 3.2)
    ax.set_aspect("equal")
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.set_title("(b)", loc="left")
    fig.tight_layout()
    return save_mpl(fig, "01-02-paraboloid-contours", outdir)


def taylor_models(outdir, lang="en"):
    """f(x) = e^x - 2x and its Taylor models at xbar = 0.

    First order: 1 - x (tangent line). Second order: 1 - x + x^2/2
    (osculating parabola), whose minimiser x = 1 is the Newton step from 0;
    the minimiser of f is ln 2 = 0.6931. Used in example
    ex-01-02-taylor-models of 01-02-matrix-calculus.md.
    """
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    x = np.linspace(-1.2, 2.2, 300)
    f = np.exp(x) - 2 * x
    lin = 1 - x
    quad = 1 - x + x**2 / 2
    ax.plot(x, f, color=INK, lw=2.0, label=r"$f(x)=e^x-2x$")
    ax.plot(
        x,
        lin,
        color=BLUE,
        lw=1.5,
        ls="--",
        label=_t(lang, "first order", "一阶模型"),
    )
    ax.plot(
        x,
        quad,
        color=ORANGE,
        lw=1.5,
        label=_t(lang, "second order", "二阶模型"),
    )
    ax.plot(0, 1, "o", color=INK, ms=6, zorder=5)
    ax.annotate(
        r"$\bar{x}=0$",
        xy=(0, 1),
        xytext=(-0.75, 0.42),
        fontsize=11,
        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8),
    )
    # minimisers of f and of the quadratic model
    xm = np.log(2)
    ax.plot(xm, np.exp(xm) - 2 * xm, "o", color=GREEN, ms=6, zorder=5)
    ax.plot(1, 0.5, "o", color=ORANGE, ms=6, zorder=5)
    ax.annotate(
        r"$\ln 2$",
        xy=(xm, np.exp(xm) - 2 * xm),
        xytext=(0.45, 1.25),
        fontsize=11,
        color=GREEN,
        arrowprops=dict(arrowstyle="-", color=GREEN, lw=0.8),
    )
    ax.annotate(
        r"$x=1$",
        xy=(1, 0.5),
        xytext=(1.25, 0.05),
        fontsize=11,
        color=ORANGE,
        arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.8),
    )
    ax.set_xlim(-1.2, 2.2)
    ax.set_ylim(-0.6, 3.4)
    ax.set_xlabel("$x$")
    ax.legend(loc="upper center", fontsize=10)
    return save_mpl(fig, "01-02-taylor-models", outdir)


FIGURES = {
    "01-02-linear-level-lines": linear_level_lines,
    "01-02-paraboloid-contours": paraboloid_contours,
    "01-02-taylor-models": taylor_models,
}
