"""Figures of chapter 7.2, Linearization Techniques.

* ``07-02-mccormick``: the product w = xy on the box [1, 4] x [1, 3] (Example
  ex-07-02-mccormick) along the two diagonals of the box, with the band between
  the McCormick under- and over-estimators;
* ``07-02-abs-max``: three piecewise-linear distance functions of the corridor
  example ex-07-02-corridor (weighted sum, maximum and minimum of |x - a_i|).
"""

import matplotlib.pyplot as plt
import numpy as np

from .common import BLUE, GREEN, INK, MUTED, ORANGE, save_mpl

# Box of the McCormick example
XL, XU, YL, YU = 1.0, 4.0, 1.0, 3.0

# Corridor example: town positions (km) and populations (thousand)
TOWNS = np.array([0, 12, 20, 34, 50], float)
POP = np.array([8, 5, 12, 4, 6], float)
LENGTH = 50.0


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


def _envelope(x, y):
    """McCormick under-estimator (max of two planes) and over-estimator
    (min of two planes) of x*y on the box [XL, XU] x [YL, YU]."""
    under = np.maximum(YL * x + XL * y - XL * YL, YU * x + XU * y - XU * YU)
    over = np.minimum(YL * x + XU * y - XU * YL, YU * x + XL * y - XL * YU)
    return under, over


def mccormick(outdir, lang="en"):
    """xy and its McCormick envelope along the two diagonals of the box.

    Used in Section sec-07-02-mccormick. Parameter s in [0, 1] runs along
    (left) the main diagonal (1,1) -> (4,3), where xy is convex, and
    (right) the other diagonal (1,3) -> (4,1), where xy is concave.
    The band touches xy at the corners and is widest (1.5) at the center.
    """
    s = np.linspace(0, 1, 301)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6), sharey=True)
    diagonals = [
        # (x(s), y(s), title)
        (
            XL + s * (XU - XL),
            YL + s * (YU - YL),
            _t(lang, r"main diagonal $(1,1)\to(4,3)$",
               r"主对角线 $(1,1)\to(4,3)$"),  # fmt: skip
        ),
        (
            XL + s * (XU - XL),
            YU - s * (YU - YL),
            _t(lang, r"other diagonal $(1,3)\to(4,1)$",
               r"另一条对角线 $(1,3)\to(4,1)$"),  # fmt: skip
        ),
    ]
    for ax, (x, y, title) in zip(axes, diagonals):
        under, over = _envelope(x, y)
        # shaded band between the envelopes, then the envelopes and xy
        ax.fill_between(s, under, over, color=GREEN, alpha=0.18, linewidth=0)
        ax.plot(
            s,
            under,
            color=GREEN,
            linestyle="--",
            linewidth=1.6,
            label=_t(lang, "under-estimator", "下估计函数"),
        )
        ax.plot(
            s,
            over,
            color=ORANGE,
            linestyle="--",
            linewidth=1.6,
            label=_t(lang, "over-estimator", "上估计函数"),
        )
        ax.plot(s, x * y, color=BLUE, linewidth=2.2, label=r"$\omega = xy$")
        # the center of the box, where the gap is largest
        ax.axvline(
            0.5, color=MUTED, linewidth=0.8, linestyle=(0, (2, 3)), zorder=0
        )
        ax.set_title(title, fontsize=10.5)
        ax.set_xlabel(
            _t(lang, "position $s$ along the diagonal", "沿对角线的位置 $s$")
        )
        ax.set_xlim(0, 1)
    axes[0].set_ylabel(r"$\omega$")
    axes[0].legend(loc="upper left", frameon=False, fontsize=9.5)
    # annotate the gap of 1.5 at the center of the main diagonal
    x0, y0 = 2.5, 2.0
    u0, o0 = _envelope(np.array(x0), np.array(y0))
    axes[0].annotate(
        "",
        xy=(0.5, float(o0)),
        xytext=(0.5, x0 * y0),
        arrowprops=dict(arrowstyle="<->", color=INK, lw=1.0),
    )
    # label in the upper part of the gap: xy rises through its middle
    axes[0].text(
        0.53,
        6.4,
        "1.5",
        ha="left",
        va="center",
        fontsize=10,
        color=INK,
    )
    axes[1].annotate(
        "",
        xy=(0.5, float(u0)),
        xytext=(0.5, x0 * y0),
        arrowprops=dict(arrowstyle="<->", color=INK, lw=1.0),
    )
    axes[1].text(
        0.53,
        (float(u0) + x0 * y0) / 2,
        "1.5",
        ha="left",
        va="center",
        fontsize=10,
        color=INK,
    )
    fig.tight_layout()
    return save_mpl(fig, "07-02-mccormick", outdir)


def abs_max(outdir, lang="en"):
    """Weighted sum, maximum and minimum of the distances |x - a_i|.

    Used in Section sec-07-02-abs (Example ex-07-02-corridor). Optima: weighted
    median x = 20 (value 436), center x = 25 (value 25), max-min point
    x = 42 (value 8).
    """
    x = np.linspace(0, LENGTH, 1001)
    dist = np.abs(x[:, None] - TOWNS[None, :])
    km = _t(lang, "km", "千米")
    panels = [
        # (values, optimum x, optimum value, title, y label)
        (
            dist @ POP,
            20,
            436,
            r"$\sum_i w_i\,|x-a_i|$",
            _t(lang, "thousand resident-km", "千人·千米"),
        ),
        (dist.max(axis=1), 25, 25, r"$\max_i\,|x-a_i|$", km),
        (dist.min(axis=1), 42, 8, r"$\min_i\,|x-a_i|$", km),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.3))
    for k, (ax, (vals, xo, vo, title, ylab)) in enumerate(zip(axes, panels)):
        ax.plot(x, vals, color=BLUE, linewidth=2.0)
        # optimum: minimum for the first two panels, maximum for the third
        ax.plot([xo], [vo], marker="o", markersize=6.5, color=ORANGE, zorder=5)
        ax.annotate(
            f"$x={xo}$",
            xy=(xo, vo),
            xytext=(0, -14) if k < 2 else (0, 8),
            textcoords="offset points",
            ha="center",
            va="top" if k < 2 else "bottom",
            fontsize=10,
        )
        # town positions as ticks on the horizontal axis
        ax.plot(
            TOWNS,
            np.zeros_like(TOWNS),
            linestyle="none",
            marker="|",
            markersize=10,
            color=INK,
            clip_on=False,
            zorder=4,
        )
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(_t(lang, "$x$ (km)", "$x$（千米）"))
        ax.set_ylabel(ylab, fontsize=9.5)
        ax.set_xlim(0, LENGTH)
        ax.set_ylim(0, vals.max() * 1.12)
        ax.set_xticks([0, 10, 20, 30, 40, 50])
    fig.tight_layout()
    return save_mpl(fig, "07-02-abs-max", outdir)


FIGURES = {
    "07-02-mccormick": mccormick,
    "07-02-abs-max": abs_max,
}
