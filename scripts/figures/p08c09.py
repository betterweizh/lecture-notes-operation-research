"""Figures of chapter 8.9, Nonlinear Programming Algorithms (Part VIII,
Examples and Problems).

* ``08-09-fw-graphical``: feasible triangle of min 4(x1-10)^2 + (x2-4)^2
  s.t. x1 - x2 <= 10, x1/5 - x2 >= 3, x1 >= 0; (a) level lines of the first
  linear subproblem, (b) contours of f and the Frank-Wolfe step
  (example ex-08-09-fw-graphical);
* ``08-09-flow-network``: the 4-node network of the convex-cost flow example
  with marginal arc costs (example ex-08-09-fw-network-flow);
* ``08-09-penalty-path``: contours of (x1-2)^4 + (x1-2x2)^2, the parabola
  x2 = x1^2 and the penalty minimizers x(rho) (example ex-08-09-penalty-parabola).
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from .common import (
    BLUE,
    GREEN,
    INK,
    MUTED,
    ORANGE,
    network_svg,
    save_mpl,
    write_svg,
)

# vertices of the feasible triangle (exact values 35/4 and -5/4)
TRI = np.array([[0.0, -3.0], [0.0, -10.0], [8.75, -1.25]])


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


def _triangle_panel(ax):
    """Shade the feasible triangle and draw its three boundary lines."""
    x = np.linspace(-2, 16, 10)
    ax.fill(TRI[:, 0], TRI[:, 1], color=BLUE, alpha=0.15, lw=0)
    ax.plot(x, x - 10, color=BLUE, lw=1.3)  # x1 - x2 = 10
    ax.plot(x, x / 5 - 3, color=ORANGE, lw=1.3)  # x1/5 - x2 = 3
    ax.axvline(0, color=GREEN, lw=1.3)  # x1 = 0
    ax.text(13.3, 2.0, r"$x_1-x_2=10$", color=BLUE, fontsize=9.5, rotation=0)
    ax.text(11.2, -2.9, r"$\frac{1}{5}x_1-x_2=3$", color=ORANGE, fontsize=9.5)
    ax.text(0.3, 3.6, r"$x_1=0$", color=GREEN, fontsize=9.5)
    # vertex labels
    ax.annotate("$(0,-3)$", (0, -3), xytext=(-6, 6), textcoords="offset points",
                ha="right", fontsize=9)  # fmt: skip
    ax.annotate("$(0,-10)$", (0, -10), xytext=(-6, -2), textcoords="offset points",
                ha="right", fontsize=9)  # fmt: skip
    ax.annotate(r"$(\frac{35}{4},-\frac{5}{4})$", (8.75, -1.25), xytext=(-6, 9),
                textcoords="offset points", ha="right", fontsize=9)  # fmt: skip
    ax.plot(TRI[:, 0], TRI[:, 1], "o", color=INK, ms=4, zorder=5)
    ax.set_xlim(-4, 17)
    ax.set_ylim(-12, 5)
    ax.set_aspect("equal")
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")


def fw_graphical(outdir, lang="en"):
    """Graphical Frank-Wolfe: LP by inspection and the single step.

    Used in chapter 8.9 (fig-08-09-fw-graphical). Panel (a): level lines
    -80 y1 - 14 y2 = c of the first linear subproblem; the minimum over the
    triangle is at (35/4, -5/4). Panel (b): contours of f, the step from
    (0, -3) to (35/4, -5/4) (alpha = 1, line minimizer 1.171 clipped) and the
    unconstrained minimizer (10, 4).
    """
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.4))

    # (a) linear subproblem: level lines of g^T y with g = (-80, -14)
    ax = axes[0]
    _triangle_panel(ax)
    y1 = np.linspace(-4, 17, 10)
    for c in (140, 42, -300, -682.5):
        # -80 y1 - 14 y2 = c  ->  y2 = (-c - 80 y1) / 14
        ax.plot(y1, (-c - 80 * y1) / 14, color=MUTED, ls="--", lw=0.9)
    # direction (80, 14) of decreasing c, right of the last level line so
    # that no dashed line crosses the label
    ax.annotate("", xy=(12.3, -10.92), xytext=(10.7, -11.2),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.2))  # fmt: skip
    ax.text(10.7, -9.5, _t(lang, "decreasing", "递减方向") + "\n$-80y_1-14y_2$",
            fontsize=8.5, color=MUTED, va="center")  # fmt: skip
    ax.set_title(
        _t(lang, "(a) linear subproblem at $(0,-3)$",
           "(a) $(0,-3)$ 处的线性化子问题"),
        fontsize=11, loc="left",
    )  # fmt: skip

    # (b) contours of f and the Frank-Wolfe step
    ax = axes[1]
    X1, X2 = np.meshgrid(np.linspace(-4, 17, 300), np.linspace(-12, 5, 300))
    F = 4 * (X1 - 10) ** 2 + (X2 - 4) ** 2
    ax.contour(X1, X2, F, levels=[4, 33.8125, 100, 250, 449], colors=MUTED,
               linewidths=0.7, alpha=0.8)  # fmt: skip
    _triangle_panel(ax)
    ax.annotate("", xy=(8.75, -1.25), xytext=(0, -3),
                arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2.4,
                                shrinkA=0, shrinkB=0, mutation_scale=14))  # fmt: skip
    ax.plot(10, 4, "+", color=INK, ms=9, mew=1.5)
    ax.annotate("$(10,4)$", (10, 4), xytext=(6, -4), textcoords="offset points",
                fontsize=9)  # fmt: skip
    ax.set_title(
        _t(lang, "(b) step", "(b) 迭代步")
        + r" $\mathbf{x}^{(0)}\to\mathbf{x}^{(1)}=\mathbf{x}^*$",
        fontsize=11, loc="left",
    )  # fmt: skip
    fig.tight_layout()
    return save_mpl(fig, "08-09-fw-graphical", outdir)


def flow_network(outdir, lang="en"):
    """Four-node network of the convex-cost flow example.

    Used in chapter 8.9 (fig-08-09-flow-network). Arc labels are the marginal
    costs A_a + B_a x (dollars per unit); arcs with positive optimal flow
    (all except (3, 2)) are highlighted.
    """
    nodes = {"1": (95, 140), "2": (255, 50), "3": (255, 230), "4": (415, 140)}
    arcs = [
        ("1", "2", "4+0.02x"),
        ("1", "3", "6+0.01x"),
        ("2", "3", "1+0.04x"),
        ("3", "2", "1+0.04x"),
        ("2", "4", "7+0.01x"),
        ("3", "4", "3+0.03x"),
    ]
    body = network_svg(
        nodes,
        arcs,
        highlight_arcs=[
            ("1", "2"),
            ("1", "3"),
            ("2", "3"),
            ("2", "4"),
            ("3", "4"),
        ],
        node_notes={"1": "+210.5", "4": "−210.5"},
        note_pos={"1": "left", "4": "right"},
    )
    return write_svg("08-09-flow-network", body, 510, 280, outdir)


def penalty_path(outdir, lang="en"):
    """Penalty minimizers x(rho) approaching the parabola x2 = x1^2.

    Used in chapter 8.9 (fig-08-09-penalty-path). Each x(rho) minimizes
    f + rho (x1^2 - x2)^2 (BFGS, warm-started from the previous rho).
    """
    f = lambda x: (x[0] - 2) ** 4 + (x[0] - 2 * x[1]) ** 2
    h = lambda x: x[0] ** 2 - x[1]
    pts, x0 = [], np.array([2.0, 1.0])
    for rho in (0.1, 1, 10, 100, 1000):
        res = minimize(lambda x, r=rho: f(x) + r * h(x) ** 2, x0, method="BFGS",
                       options={"gtol": 1e-10})  # fmt: skip
        x0 = res.x
        pts.append((rho, res.x))

    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    X1, X2 = np.meshgrid(
        np.linspace(0.4, 2.4, 300), np.linspace(0.3, 1.6, 300)
    )
    F = (X1 - 2) ** 4 + (X1 - 2 * X2) ** 2
    ax.contour(X1, X2, F, levels=[0.02, 0.1, 0.3, 0.6, 1.0, 1.946, 3],
               colors=MUTED, linewidths=0.7, alpha=0.8)  # fmt: skip
    t = np.linspace(0.4, 1.3, 200)
    ax.plot(t, t**2, color=BLUE, lw=2, label=r"$x_2=x_1^2$")
    xs = np.array([p for _, p in pts])
    ax.plot(xs[:, 0], xs[:, 1], "-", color=ORANGE, lw=1, alpha=0.7)
    ax.plot(xs[:, 0], xs[:, 1], "o", color=ORANGE, ms=5, zorder=5,
            label=r"$\mathbf{x}(\rho)$")  # fmt: skip
    for (rho, p), off in zip(pts[:4], [(6, -4), (6, -10), (8, -4), (8, 4)]):
        ax.annotate(rf"$\rho={rho:g}$", p, xytext=off, textcoords="offset points",
                    fontsize=9)  # fmt: skip
    ax.plot(0.945583, 0.894127, "*", color=INK, ms=11, zorder=6)
    ax.annotate(r"$\mathbf{x}^*$", (0.945583, 0.894127), xytext=(-22, 2),
                textcoords="offset points")  # fmt: skip
    ax.plot(2, 1, "+", color=INK, ms=9, mew=1.5)
    ax.annotate("$(2,1)$", (2, 1), xytext=(6, -12), textcoords="offset points",
                fontsize=9)  # fmt: skip
    ax.set_xlim(0.4, 2.4)
    ax.set_ylim(0.3, 1.6)
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    return save_mpl(fig, "08-09-penalty-path", outdir)


FIGURES = {
    "08-09-fw-graphical": fw_graphical,
    "08-09-flow-network": flow_network,
    "08-09-penalty-path": penalty_path,
}
