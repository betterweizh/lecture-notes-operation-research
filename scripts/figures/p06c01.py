"""Figures of chapter 6.1, Nonlinear Programming Models.

Three matplotlib figures, all computed from closed forms or from the data
of the chapter's examples:

* ``06-01-three-feasible-sets``: contours of f(x) = (x1-2)^2 + (x2-2)^2 over
  the three feasible sets of example ex-06-01-three-feasible-sets;
* ``06-01-regression-fit``: observed against fitted spending at the
  least-squares estimate of example ex-06-01-consumer-regression;
* ``06-01-x-power-2x``: f(x) = (x^2)^x on [-5, 1.7] with its local and global
  extrema and its convex/concave parts (section sec-06-01-local-global).
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
from matplotlib.patches import Polygon

from .common import BLUE, GREEN, GRID, INK, MUTED, ORANGE, PALETTE, save_mpl


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


def three_feasible_sets(outdir, lang="en"):
    """Contours of the squared distance to (2, 2) over three feasible sets.

    Panel 1: Omega = R^2, optimum (2, 2). Panel 2: triangle x1 + x2 <= 5,
    x >= 0, optimum (2, 2) inside. Panel 3: disc (x1-1)^2 + (x2-1)^2 <= 1
    cut by x1 + x2 >= 1, optimum (1 + 1/sqrt2, 1 + 1/sqrt2) on the circle.
    Used in ex-06-01-three-feasible-sets of 06-01-nlp-models.md.
    """
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.5))
    # common contour grid of f(x) = (x1-2)^2 + (x2-2)^2
    g = np.linspace(-0.6, 5.6, 400)
    X1, X2 = np.meshgrid(g, g)
    F = (X1 - 2) ** 2 + (X2 - 2) ** 2
    s = 1 / np.sqrt(2)
    # (window, optimum, contour levels) per panel
    panels = [
        ((-0.5, 4.5), (2.0, 2.0), [0.25, 1, 2.25, 4, 6.25]),
        ((-0.5, 5.6), (2.0, 2.0), [0.25, 1, 2.25, 4, 6.25, 9]),
        ((-0.3, 2.6), (1 + s, 1 + s), [3 - 2 * np.sqrt(2), 0.5, 1, 2, 3]),
    ]
    for k, (ax, (lim, opt, levels)) in enumerate(zip(axes, panels)):
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect("equal")
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.axvline(0, color=MUTED, lw=0.8)
        # feasible set of each panel
        if k == 0:
            # Omega = R^2: shade the whole window lightly
            ax.add_patch(
                Polygon(
                    [
                        (lim[0], lim[0]),
                        (lim[1], lim[0]),
                        (lim[1], lim[1]),
                        (lim[0], lim[1]),
                    ],
                    closed=True,
                    color=BLUE,
                    alpha=0.08,
                    lw=0,
                )
            )
            ax.text(
                3.3, 3.9, r"$\Omega=\mathbb{R}^2$", color=BLUE, ha="center"
            )
        elif k == 1:
            # triangle (0,0), (5,0), (0,5)
            ax.add_patch(
                Polygon(
                    [(0, 0), (5, 0), (0, 5)],
                    closed=True,
                    fc=to_rgba(BLUE, 0.15),  # light fill, opaque edge
                    lw=1.4,
                    ec=BLUE,
                )
            )
            ax.text(3.0, 0.6, r"$\Omega$", color=BLUE, fontsize=13)
            # line label just outside the triangle, below the contours' center
            ax.text(
                4.55,
                1.0,
                r"$x_1+x_2=5$",
                color=BLUE,
                fontsize=10,
                rotation=-45,
                ha="center",
                va="center",
            )
        else:
            # disc centered at (1,1), radius 1, minus the corner below x1+x2=1
            t = np.linspace(0, 2 * np.pi, 400)
            cx, cy = 1 + np.cos(t), 1 + np.sin(t)
            keep = cx + cy >= 1
            # boundary of the feasible part: arc where x1+x2 >= 1, closed by the chord
            ang = np.linspace(
                -np.pi / 2, np.pi, 400
            )  # from (1,0) round to (0,1)
            arc = np.c_[1 + np.cos(ang), 1 + np.sin(ang)]
            ax.add_patch(
                Polygon(
                    arc, closed=True, fc=to_rgba(BLUE, 0.15), lw=1.4, ec=BLUE
                )
            )
            # the cut-off corner (dotted circle part) and the cutting line
            ax.plot(cx[~keep], cy[~keep], color=BLUE, lw=1.0, ls=":")
            ax.plot([-0.3, 1.3], [1.3, -0.3], color=MUTED, lw=1.0, ls="--")
            ax.text(0.9, 0.95, r"$\Omega$", color=BLUE, fontsize=13)
            # label at the lower end of the cutting line, clear of the dotted arc
            ax.text(
                1.36,
                -0.2,
                r"$x_1+x_2=1$",
                color=MUTED,
                fontsize=9,
                ha="left",
                va="center",
            )
        # contours of f
        ax.contour(X1, X2, F, levels=levels, colors=MUTED, linewidths=0.8)
        # point (2, 2) and the optimum
        ax.plot(2, 2, "o", ms=4, color=INK, zorder=4)
        ax.plot(*opt, "o", ms=7, color=ORANGE, zorder=5)
        if k < 2:
            ax.annotate(
                r"$\mathbf{x}^*=(2,2)$",
                xy=opt,
                xytext=(-8, 8),
                textcoords="offset points",
                color=ORANGE,
                ha="right",
            )
        else:
            ax.annotate(
                r"$(2,2)$",
                xy=(2, 2),
                xytext=(-8, 6),
                textcoords="offset points",
                ha="right",
            )
            ax.annotate(
                r"$\mathbf{x}^*\approx(1.707,1.707)$",
                xy=opt,
                xytext=(-8, -16),
                textcoords="offset points",
                ha="right",
                color=ORANGE,
                fontsize=10,
            )
        ax.set_title(
            _t(lang, f"Problem {k + 1}", f"问题 {k + 1}"), fontsize=11
        )
        ax.set_xlabel("$x_1$")
        if k == 0:
            ax.set_ylabel("$x_2$")
        ax.grid(False)
    fig.tight_layout()
    return save_mpl(fig, "06-01-three-feasible-sets", outdir)


# --- data of the regression example (prices, quantities, budgets) ---------
E = np.array([1500.0, 2000.0, 1500.0, 1200.0])  # budget E_t ($)
P = np.array(
    [
        [2.5, 4.0, 5.0],  # prices p_ti ($/unit), rows = months
        [2.7, 3.5, 4.9],
        [2.8, 4.0, 5.6],
        [2.6, 3.6, 4.5],
    ]
)
Q = np.array(
    [
        [198, 124, 100],  # purchased quantities Q_ti (units)
        [200, 200, 110],
        [100, 80, 60],
        [152, 115, 65],
    ],
    dtype=float,
)
# least-squares estimate (multistart SLSQP, see the chapter)
BETA = np.array([0.24737933, 0.46354975, 0.28907091])
QSTAR = np.array([89.85027829, 16.67488124, 37.3199961])


def regression_fit(outdir, lang="en"):
    """Observed against fitted spending for 12 product-month pairs.

    Fitted spending p_ti q_i + beta_i (E_t - sum_j p_tj q_j) at the estimate
    (BETA, QSTAR); observed spending p_ti Q_ti. One color per product,
    month numbers written next to the points. Used in
    ex-06-01-consumer-regression of 06-01-nlp-models.md.
    """
    fitted = P * QSTAR + np.outer(E - P @ QSTAR, BETA)
    observed = P * Q
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    lo, hi = 250, 800
    # 45-degree line: perfect fit
    ax.plot([lo, hi], [lo, hi], color=MUTED, ls="--", lw=1.0)
    markers = ["o", "s", "^"]
    for i in range(3):
        ax.scatter(
            observed[:, i],
            fitted[:, i],
            s=36,
            color=PALETTE[i],
            marker=markers[i],
            label=_t(lang, f"product {i + 1}", f"产品 {i + 1}"),
            zorder=3,
        )
        # month number next to each point (offset per product so that
        # labels of neighboring points do not collide)
        off, ha = [((0, -13), "center"), ((7, -4), "left"), ((7, 4), "left")][
            i
        ]
        for t in range(4):
            o, h = off, ha
            if (i, t) == (0, 3):
                # product 1, month 4: left of the point, off the 45-degree line
                o, h = (-8, -3), "right"
            ax.annotate(
                str(t + 1),
                (observed[t, i], fitted[t, i]),
                xytext=o,
                textcoords="offset points",
                ha=h,
                fontsize=9,
                color=PALETTE[i],
            )
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel(
        _t(
            lang,
            r"observed spending $p_{ti}Q_{ti}$ (\$)",
            r"观测支出 $p_{ti}Q_{ti}$（美元）",
        )
    )
    ax.set_ylabel(_t(lang, r"fitted spending (\$)", "拟合支出（美元）"))
    ax.legend(loc="lower right", frameon=False, fontsize=10)
    # the three month-3 points form the cluster at the upper left
    ax.text(262, 555, _t(lang, "month 3", "第 3 个月"), color=INK, fontsize=10)
    return save_mpl(fig, "06-01-regression-fit", outdir)


def x_power_2x(outdir, lang="en"):
    """f(x) = (x^2)^x on [-5, 1.7] with extrema and curvature regions.

    Local max at -1/e (value e^{2/e}), local min at 1/e (value e^{-2/e}),
    global min at x = -5 (value 25^-5 ~ 1e-7), global max at 1.7. f is
    concave on (-0.808, 0) (shaded) and convex elsewhere. Used in
    sec-06-01-local-global of 06-01-nlp-models.md.
    """
    x = np.linspace(-5, 1.7, 2001)
    x = x[np.abs(x) > 1e-9]
    f = np.exp(2 * x * np.log(np.abs(x)))
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    # concave part (-0.808, 0) shaded
    ax.axvspan(-0.8079, 0.0, color=GREEN, alpha=0.12, lw=0)
    concave = _t(lang, "concave", "凹区间")
    convex = _t(lang, "convex", "凸区间")
    ax.text(-0.40, 5.6, concave, ha="center", fontsize=9, color=GREEN)
    ax.text(-2.9, 5.6, convex, ha="center", fontsize=9, color=MUTED)
    ax.text(0.85, 5.6, convex, ha="center", fontsize=9, color=MUTED)
    ax.plot(x, f, color=BLUE, lw=1.8)
    # the value f(0) = 1 (limit), drawn as a filled point
    ax.plot(0, 1, "o", ms=3.5, color=BLUE)
    # extrema: (x, f, label, text offset)
    e = np.e
    pts = [
        (
            -5,
            25.0**-5,
            _t(lang, "global min", "全局极小点") + "\n$x=-5$",
            (6, 14),
            "left",
        ),
        (
            -1 / e,
            e ** (2 / e),
            _t(lang, "local max", "局部极大点") + "\n$x=-1/e$",
            (-8, 8),
            "right",
        ),
        (
            1 / e,
            e ** (-2 / e),
            _t(lang, "local min", "局部极小点") + "\n$x=1/e$",
            (14, -14),
            "left",
        ),
        (
            1.7,
            1.7**3.4,
            _t(lang, "global max", "全局极大点") + "\n$x=1.7$",
            (-8, -4),
            "right",
        ),
    ]
    for px, py, lab, off, ha in pts:
        ax.plot(px, py, "o", ms=6, color=ORANGE, zorder=5)
        ax.annotate(
            lab,
            (px, py),
            xytext=off,
            textcoords="offset points",
            ha=ha,
            fontsize=9,
            color=INK,
        )
    ax.set_xlim(-5.2, 1.9)
    ax.set_ylim(-0.4, 6.3)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$f(x)=(x^2)^x$")
    ax.axhline(0, color=GRID, lw=0.8)
    return save_mpl(fig, "06-01-x-power-2x", outdir)


FIGURES = {
    "06-01-three-feasible-sets": three_feasible_sets,
    "06-01-regression-fit": regression_fit,
    "06-01-x-power-2x": x_power_2x,
}
