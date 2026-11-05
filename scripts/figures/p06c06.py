"""Figures of chapter 6.6, Algorithms for Constrained Optimization.

* ``06-06-frank-wolfe``: the triangle 3x1 + 2x2 <= 6, x >= 0, contours of
  f(x) = -5x1 + x1^2 - 8x2 + 2x2^2 and the Frank–Wolfe iterates
  x(0), ..., x(9) from (0, 0) (example ex-06-06-frank-wolfe);
* ``06-06-penalty-family``: the penalized objectives
  x + rho max(0, 2 - x)^2 for several rho and their minimizers
  x(rho) = 2 - 1/(2 rho) (example ex-06-06-penalty-1d).
"""

import matplotlib.pyplot as plt
import numpy as np

from .common import BLUE, GREEN, INK, MUTED, ORANGE, PALETTE, save_mpl

# vertices of the feasible triangle 3x1 + 2x2 <= 6, x >= 0
VERTS = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 3.0]])
HESS = np.diag([2.0, 4.0])


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


def _f(x1, x2):
    return -5 * x1 + x1**2 - 8 * x2 + 2 * x2**2


def _grad(x):
    return np.array([2 * x[0] - 5, 4 * x[1] - 8])


def _frank_wolfe(x0, n):
    """n Frank–Wolfe iterations with exact line search on [0, 1].

    Returns the iterates x(0..n) and the LP vertices y(0..n-1).
    """
    xs, ys = [np.array(x0, float)], []
    for _ in range(n):
        x = xs[-1]
        g = _grad(x)
        y = VERTS[np.argmin(VERTS @ g)]  # LP over the triangle: best vertex
        d = y - x
        alpha = min(1.0, -(g @ d) / (d @ HESS @ d))  # exact step, clipped
        ys.append(y)
        xs.append(x + alpha * d)
    return np.array(xs), np.array(ys)


def frank_wolfe(outdir, lang="en"):
    """Frank–Wolfe iterates on the triangle with contours of f.

    Used in chapter 6.6 (fig-06-06-frank-wolfe). The unconstrained minimizer
    (2.5, 2) lies outside; the optimum x* = (1, 1.5) is on the hypotenuse.
    """
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    # feasible triangle
    ax.fill(VERTS[:, 0], VERTS[:, 1], color=BLUE, alpha=0.12, lw=0)
    ax.plot(*np.vstack([VERTS, VERTS[:1]]).T, color=BLUE, lw=1.4)
    ax.text(1.78, 0.78, r"$3x_1+2x_2=6$", color=BLUE, rotation=-56, fontsize=10,
            ha="center", va="center")  # fmt: skip
    ax.text(0.45, 0.35, r"$\Omega$", color=BLUE, fontsize=13)
    # contours of f around its unconstrained minimizer (2.5, 2)
    X1, X2 = np.meshgrid(
        np.linspace(-0.4, 3.4, 300), np.linspace(-0.4, 3.4, 300)
    )
    cs = ax.contour(X1, X2, _f(X1, X2), levels=[-13.5, -12.5, -11.5, -10.5, -8, -4],
                    colors=MUTED, linewidths=0.7, alpha=0.8,
                    linestyles="solid")  # fmt: skip
    # contour labels with a true minus sign, placed by hand on the right
    # (manual positions: nearest point of the labeled contour is used)
    ax.clabel(
        cs,
        levels=[-13.5, -10.5, -8, -4],
        fontsize=7.5,
        fmt=lambda v: f"{v:g}".replace("-", "\u2212"),
        inline_spacing=2,
        manual=[(1.75, 2.6), (2.95, 0.68), (2.85, 0.26), (2.9, -0.2)],
    )
    # Frank–Wolfe path
    xs, ys = _frank_wolfe((0.0, 0.0), 9)
    for k in range(len(ys)):
        # dotted line towards the LP vertex y(k) for the first steps
        if k < 4:
            ax.plot([xs[k, 0], ys[k, 0]], [xs[k, 1], ys[k, 1]], color=MUTED,
                    ls=":", lw=1)  # fmt: skip
        ax.annotate("", xy=xs[k + 1], xytext=xs[k],
                    arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.4,
                                    shrinkA=0, shrinkB=1, mutation_scale=10))  # fmt: skip
    ax.plot(xs[:, 0], xs[:, 1], "o", color=ORANGE, ms=3.5, zorder=5)
    for k, off in ((0, (6, -12)), (1, (6, 4)), (2, (8, -6)), (3, (-30, 4))):
        ax.annotate(rf"$\mathbf{{x}}^{{({k})}}$", xs[k], xytext=off,
                    textcoords="offset points", fontsize=10)  # fmt: skip
    # optimum and unconstrained minimizer
    ax.plot(1, 1.5, "*", color=INK, ms=11, zorder=6)
    ax.annotate(r"$\mathbf{x}^*=(1,1.5)$", (1, 1.5), xytext=(10, 6),
                textcoords="offset points", fontsize=10)  # fmt: skip
    ax.plot(2.5, 2, "+", color=MUTED, ms=9, mew=1.5)
    ax.annotate("(2.5, 2)", (2.5, 2), xytext=(6, 4), textcoords="offset points",
                fontsize=9, color=MUTED)  # fmt: skip
    ax.set_xlim(-0.3, 3.3)
    ax.set_ylim(-0.3, 3.3)
    ax.set_aspect("equal")
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    fig.tight_layout()
    return save_mpl(fig, "06-06-frank-wolfe", outdir)


def penalty_family(outdir, lang="en"):
    """Penalized objectives x + rho max(0, 2 - x)^2 for several rho.

    Used in chapter 6.6 (fig-06-06-penalty-family). Minimizers
    x(rho) = 2 - 1/(2 rho): 1, 1.5, 1.75, 1.9 for rho = 0.5, 1, 2, 5.
    """
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    x = np.linspace(0, 3, 600)
    ax.axvspan(2, 3, color=GREEN, alpha=0.1, lw=0)  # feasible side x >= 2
    ax.text(
        2.5,
        0.25,
        _t(lang, "feasible", "可行域") + "\n$x\\geq2$",
        ha="center",
        color=GREEN,
        fontsize=10,
    )
    ax.plot(x, x, color=INK, lw=1.2, ls="--", label="$f(x)=x$")
    for rho, col in zip(
        (0.5, 1, 2, 5), PALETTE[:1] + PALETTE[1:2] + PALETTE[3:5]
    ):
        phi = x + rho * np.maximum(0, 2 - x) ** 2
        ax.plot(x, phi, color=col, lw=1.8, label=rf"$\rho={rho:g}$")
        xr = 2 - 1 / (2 * rho)  # minimizer x(rho)
        ax.plot(xr, xr + rho * (2 - xr) ** 2, "o", color=col, ms=6, zorder=5)
    ax.axvline(2, color=MUTED, lw=1)
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 4)
    ax.set_xlabel("$x$")
    ax.set_ylabel(r"$\Phi_\rho(x)$")
    # legend in the empty part of the feasible side, clear of the curves
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.2), fontsize=9)
    fig.tight_layout()
    return save_mpl(fig, "06-06-penalty-family", outdir)


FIGURES = {
    "06-06-frank-wolfe": frank_wolfe,
    "06-06-penalty-family": penalty_family,
}
