"""Figures of chapter 2.3, Convex Optimization and Lagrange Duality.

Two matplotlib figures from closed-form data (numbers checked in the
verification script of chapter 2.3):

* ``02-03-dispatch-dual``: the dual function q(mu) of the economic-dispatch
  example (ex-02-03-dispatch-dual);
* ``02-03-duality-geometry``: the set G = {(g(x), f(x))} with supporting
  lines, for the capped dispatch problem (convex, no gap) and for the
  on/off pump (nonconvex, gap 0.5) (section sec-02-03-geometry).
"""

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np

from .common import BLUE, INK, MUTED, ORANGE, save_mpl


HALO = [pe.withStroke(linewidth=3, foreground="white")]


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


# Dispatch data: C1 = 8 x1 + 0.01 x1^2, C2 = 10 x2 + 0.02 x2^2 ($/h), D = 475 MW.
A1, B1, A2, B2, DEMAND = 8.0, 0.01, 10.0, 0.02, 475.0


def q_dispatch(mu):
    """Dual function of the uncapped dispatch problem, eq-02-03-dispatch-q."""
    return DEMAND * mu - (mu - A1) ** 2 / (4 * B1) - (mu - A2) ** 2 / (4 * B2)


def cost_reduced(x1):
    """Total cost with x2 = 475 - x1 eliminated."""
    x2 = DEMAND - x1
    return A1 * x1 + B1 * x1**2 + A2 * x2 + B2 * x2**2


def dispatch_dual(outdir, lang="en"):
    """q(mu) for mu in [4, 26] with p* = 5587.5 and the maximiser mu* = 15.

    Values: q(10) = q(20) = 4650, q(15) = 5587.5 = p*.
    """
    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    mu = np.linspace(4, 26, 300)
    ax.plot(mu, q_dispatch(mu), color=BLUE, lw=2.0)
    p_star = 5587.5
    ax.axhline(p_star, color=MUTED, lw=1.2, ls="--")
    ax.text(
        25.8,
        p_star + 60,
        r"$p^*=5587.5$",
        ha="right",
        va="bottom",
        color=MUTED,
    )
    # the maximiser
    ax.plot(15, q_dispatch(15), "o", color=ORANGE, ms=6, zorder=5)
    ax.plot([15, 15], [1000, q_dispatch(15)], color=ORANGE, lw=1.0, ls=":")
    ax.annotate(
        r"$(\mu^*,\,d^*)=(15,\ 5587.5)$", xy=(15, q_dispatch(15)), xytext=(16.0, 5950),
        ha="left", va="bottom", arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9),
    )  # fmt: skip
    # two sample lower bounds
    for m in (10, 20):
        ax.plot(m, q_dispatch(m), "o", color=BLUE, ms=4, zorder=5)
    ax.text(9.6, q_dispatch(10) + 40, r"$q(10)=4650$", va="bottom", ha="right")
    ax.text(20.4, q_dispatch(20) + 40, r"$q(20)=4650$", va="bottom", ha="left")
    ax.text(23.6, q_dispatch(23.6) + 150, r"$q(\mu)$", color=BLUE, ha="left")
    ax.set_xlim(4, 26)
    ax.set_ylim(1000, 6400)
    ax.set_xlabel(_t(lang, r"$\mu$ (\$/MWh)", r"$\mu$（美元/MWh）"))
    ax.set_ylabel(_t(lang, r"\$/h", "美元/小时"))
    return save_mpl(fig, "02-03-dispatch-dual", outdir)


def duality_geometry(outdir, lang="en"):
    """Set G and supporting lines t + lambda u = q(lambda).

    (a) capped dispatch reduced to x1: G = {(x1 - 300, C(x1))}, a parabola
    in (u, t) with vertex (50, 5587.5). Supporting line with lambda* = 3
    touches G at (0, 5662.5) = (0, p*); the line with lambda = 1 touches G
    at u = 33.33 and has intercept q(1) = 5629.17 < p*.
    (b) on/off pump: G = {(0.5, 0), (-0.5, 1)}; the best line (lambda = 1)
    passes through both points and meets u = 0 at d* = 0.5 < p* = 1.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.9))

    # (a) parabola G over x1 in [150, 450]
    x1 = np.linspace(150, 450, 300)
    u = x1 - 300
    t = cost_reduced(x1)
    ax1.plot(u, t, color=INK, lw=2.0)  # the curve G (named in the caption)
    ax1.axvspan(-160, 0, color=MUTED, alpha=0.06, lw=0)  # feasible side u <= 0
    ax1.axvline(0, color=MUTED, lw=1.0)
    us = np.array([-160, 160])
    # lambda* = 3: t = 5662.5 - 3u, tangent at u = 0
    ax1.plot(us, 5662.5 - 3 * us, color=ORANGE, lw=1.8)
    # the curve G crosses this label: white halo keeps it legible
    ax1.text(-155, 5662.5 + 3 * 155 + 15, _t(lang, r"slope $-\lambda^*=-3$", r"斜率 $-\lambda^*=-3$"), color=ORANGE,
             ha="left", va="bottom", path_effects=HALO)  # fmt: skip
    # lambda = 1: intercept q(1) = 5629.17, tangent at u = 33.33
    q1 = cost_reduced(300 + 100 / 3) + 100 / 3
    ax1.plot(us, q1 - us, color=BLUE, lw=1.4, ls="--")
    ax1.text(
        -155,
        q1 + 155 + 12,
        _t(lang, r"slope $-1$", r"斜率 $-1$"),
        color=BLUE,
        ha="left",
        va="bottom",
    )
    ax1.plot(0, 5662.5, "o", color=ORANGE, ms=6, zorder=5)
    ax1.plot(0, q1, "o", color=BLUE, ms=5, zorder=5)
    ax1.annotate(
        r"$p^*=d^*=q(3)$", xy=(0, 5662.5), xytext=(12, 5990),
        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9),
    )  # fmt: skip
    ax1.annotate(
        r"$q(1)$", xy=(0, q1), xytext=(-100, 5520),
        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9),
    )  # fmt: skip
    ax1.set_xlim(-160, 160)
    ax1.set_ylim(5450, 6300)
    ax1.set_xlabel(r"$u=g(\mathbf{x})$ (MW)")
    ax1.set_ylabel(
        _t(
            lang,
            r"$t=f(\mathbf{x})$ (\$/h)",
            r"$t=f(\mathbf{x})$（美元/小时）",
        )
    )
    ax1.set_title("(a)")

    # (b) two points; feasible side u <= 0 shaded
    pts = np.array([[0.5, 0.0], [-0.5, 1.0]])
    ax2.axvspan(-1.0, 0, color=MUTED, alpha=0.06, lw=0)
    ax2.axvline(0, color=MUTED, lw=1.0)
    us = np.array([-1.0, 1.0])
    ax2.plot(
        us, 0.5 - us, color=ORANGE, lw=1.8
    )  # lambda = 1 through both points
    ax2.plot(*pts.T, "o", color=INK, ms=7, zorder=5)
    ax2.text(0.55, 0.05, "$x=0$", va="bottom")
    ax2.text(-0.45, 1.05, "$x=1$", va="bottom")
    ax2.plot(0, 1, "s", mfc="white", mec=INK, ms=6, zorder=6)
    ax2.plot(0, 0.5, "o", color=ORANGE, ms=6, zorder=6)
    ax2.text(0.06, 1.0, r"$p^*=1$", va="center")
    ax2.text(0.06, 0.5, r"$d^*=0.5$", va="center", color=ORANGE)
    ax2.set_xlim(-1.0, 1.0)
    ax2.set_ylim(-0.3, 1.4)
    ax2.set_xlabel(r"$u=g(x)$")
    ax2.set_ylabel(r"$t=f(x)$")
    ax2.set_title("(b)")
    return save_mpl(fig, "02-03-duality-geometry", outdir)


FIGURES = {
    "02-03-dispatch-dual": dispatch_dual,
    "02-03-duality-geometry": duality_geometry,
}
