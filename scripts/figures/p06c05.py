"""Figures of chapter 6.5, Algorithms for Unconstrained Optimization.

* ``06-05-zigzag``: contours of f(x) = (x1 - 2)^2 + 10 (x2 - 2)^2 with
  (a) the steepest-descent path from (-4, -3) (example ex-06-05-steepest-descent),
  the Newton step (ex-06-05-newton-quadratic) and the conjugate-gradient path,
  and (b) the steepest-descent path from the worst-case start (-8, 1)
  (theorem thm-06-05-sd-rate).
"""

import matplotlib.pyplot as plt
import numpy as np

from .common import BLUE, GREEN, INK, MUTED, ORANGE, save_mpl

# f(x) = 1/2 x^T Q x - b^T x + const with Q = diag(2, 20), b = (4, 40)
Q = np.diag([2.0, 20.0])
B = np.array([4.0, 40.0])
XSTAR = np.array([2.0, 2.0])


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


def _f(x1, x2):
    return (x1 - 2) ** 2 + 10 * (x2 - 2) ** 2


def _steepest_descent(x0, n):
    """n steepest-descent steps with the exact step g^T g / g^T Q g."""
    xs = [np.array(x0, float)]
    for _ in range(n):
        g = Q @ xs[-1] - B
        alpha = (g @ g) / (g @ Q @ g)
        xs.append(xs[-1] - alpha * g)
    return np.array(xs)


def _conjugate_gradient(x0):
    """Fletcher-Reeves CG with exact steps (two steps reach x* here)."""
    x = np.array(x0, float)
    g = Q @ x - B
    d = -g
    xs = [x.copy()]
    for _ in range(2):
        alpha = -(g @ d) / (d @ Q @ d)
        x = x + alpha * d
        gn = Q @ x - B
        d = -gn + (gn @ gn) / (g @ g) * d
        g = gn
        xs.append(x.copy())
    return np.array(xs)


def _contours(ax, xlim, ylim, levels):
    """Draw muted contour lines of f in the window."""
    X1, X2 = np.meshgrid(np.linspace(*xlim, 400), np.linspace(*ylim, 400))
    ax.contour(X1, X2, _f(X1, X2), levels=levels, colors=MUTED,
               linewidths=0.7, alpha=0.7)  # fmt: skip


def _path(ax, xs, color, label, lw=1.6, ms=4, zorder=4):
    """Draw a polyline of iterates with arrow heads on every segment."""
    ax.plot(xs[:, 0], xs[:, 1], "o", color=color, ms=ms, zorder=zorder)
    for p, q in zip(xs[:-1], xs[1:]):
        ax.annotate("", xy=q, xytext=p, zorder=zorder,
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                    shrinkA=0, shrinkB=2, mutation_scale=11))  # fmt: skip
    ax.plot([], [], color=color, lw=lw, label=label)  # legend entry


def zigzag(outdir, lang="en"):
    """Steepest descent, Newton and CG on an ill-conditioned quadratic.

    Used in chapter 6.5 (fig-06-05-zigzag). Panel (a): start (-4, -3); steepest
    descent (9 iterates), one Newton step straight to (2, 2), CG / BFGS in two
    steps. Panel (b): worst-case start (-8, 1), where every step reduces f by
    the factor ((kappa - 1)/(kappa + 1))^2 = 0.6694.
    """
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0),
                             gridspec_kw={"width_ratios": [1.15, 1]})  # fmt: skip

    # (a) start (-4, -3)
    ax = axes[0]
    xlim, ylim = (-5, 3.5), (-3.6, 3.0)
    _contours(ax, xlim, ylim, [0.3, 3, 29.12, 100, 286])
    sd = _t(lang, "steepest descent", "最速下降法")
    _path(ax, _steepest_descent((-4, -3), 8), BLUE, sd)
    _path(ax, np.array([[-4, -3], XSTAR]), ORANGE,
          _t(lang, "Newton (one step)", "牛顿法（一步）"), lw=1.8)  # fmt: skip
    _path(
        ax,
        _conjugate_gradient((-4, -3)),
        GREEN,
        _t(lang, "CG / BFGS (two steps)", "CG / BFGS（两步）"),
        lw=1.4,
    )
    ax.plot(*XSTAR, "*", color=INK, ms=11, zorder=6)
    ax.annotate(
        r"$\mathbf{x}^*$", XSTAR, xytext=(6, 6), textcoords="offset points"
    )
    ax.annotate(r"$\mathbf{x}^{(0)}$", (-4, -3), xytext=(6, -4),
                textcoords="offset points")  # fmt: skip
    ax.annotate(r"$\mathbf{x}^{(1)}$", (-3.3922, 2.0647), xytext=(-6, 8),
                textcoords="offset points", ha="center")  # fmt: skip
    ax.annotate(r"$\mathbf{x}^{(2)}$", (1.3891, 1.4909), xytext=(4, -18),
                textcoords="offset points", ha="left")  # fmt: skip
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.set_title(_t(lang, r"(a) start $(-4,-3)$", r"(a) 初始点 $(-4,-3)$"),
                 fontsize=11, loc="left")  # fmt: skip
    ax.legend(loc="lower right", fontsize=8.5, frameon=True, framealpha=0.85,
              edgecolor="none")  # fmt: skip

    # (b) worst-case start (-8, 1)
    ax = axes[1]
    xlim, ylim = (-9, 3.5), (0.4, 3.6)
    _contours(ax, xlim, ylim, [1, 4, 9.8, 22.09, 49.29, 110])
    _path(ax, _steepest_descent((-8, 1), 14), BLUE, sd, ms=3)
    ax.plot(*XSTAR, "*", color=INK, ms=11, zorder=6)
    ax.annotate(r"$(-8,1)$", (-8, 1), xytext=(4, -14), textcoords="offset points",
                ha="left")  # fmt: skip
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    # unequal axis scales here (x2 stretched), so right angles look oblique
    ax.set_ylabel("$x_2$")
    ax.set_xlabel("$x_1$")
    ax.set_title(
        _t(lang, r"(b) worst-case start $(-8,1)$", r"(b) 最坏初始点 $(-8,1)$"),
        fontsize=11,
        loc="left",
    )
    fig.tight_layout()
    return save_mpl(fig, "06-05-zigzag", outdir)


FIGURES = {"06-05-zigzag": zigzag}
