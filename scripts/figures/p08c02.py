"""Figures of chapter 8.2, The Simplex Method (Part VIII).

08-02-free-variable-region   feasible region with a free variable x1 (ex-08-02-a3-q3)
08-02-cut-cube               unit cube with the corner (1,1,1) cut off (ex-08-02-a3-q4/q5)
08-02-degenerate-vertex      pentagon with a degenerate corner C (ex-08-02-a3-q12)
08-02-temporary-degeneracy   region with the degenerate corner (1,0) (ex-08-02-a3-q13)
08-02-optimal-ray            unbounded region with a ray of optima (ex-08-02-a3-q14)

The 3-D cube and the pentagon redraw undimensioned textbook sketches with
data chosen to match them (stated in the chapter).
"""

import matplotlib.pyplot as plt

from .common import (
    BLUE,
    GREEN,
    INK,
    MUTED,
    ORANGE,
    PALETTE,
    draw_line,
    lp_region,
    save_mpl,
)


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


def _halo(fig):
    """White outline behind every text in the axes so that labels stay legible
    where they cross constraint lines (inverted with the figure in dark mode).
    """
    import matplotlib.patheffects as pe

    for ax in fig.axes:
        for t in ax.texts:
            t.set_path_effects(
                [pe.withStroke(linewidth=3, foreground="white")]
            )


PURPLE = PALETTE[3]


def _point(
    ax,
    p,
    text,
    offset,
    ha="left",
    va="bottom",
    color=INK,
    hollow=False,
    marker="o",
    size=6,
):
    """Mark point p and write text at an offset (in points)."""
    ax.plot(
        *p,
        marker=marker,
        markersize=size,
        zorder=5,
        clip_on=False,
        markerfacecolor="none" if hollow else color,
        markeredgecolor=color,
        markeredgewidth=1.3,
    )
    if text:
        ax.annotate(
            text,
            xy=p,
            xytext=offset,
            textcoords="offset points",
            ha=ha,
            va=va,
            fontsize=10.5,
            color=INK,
            zorder=6,
        )


def free_variable_region(outdir, lang="en"):
    """Region of max x1 + 3 x2, x1 + 2 x2 <= 2, -x1 + x2 <= 4, x2 >= 0, x1 free.

    Used in ex-08-02-a3-q3 (fig-08-02-free-variable-region). The triangle has
    corners (-4, 0), (2, 0), (-2, 2); the BFSs (0, 0) and (0, 1) with
    x1+ = x1- = 0 are not corners (squares). Dashed: contour x1 + 3 x2 = 4.
    """
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    xlim, ylim = (-5, 3), (-0.4, 3.4)
    # x1 is free, so nonneg=False and x2 >= 0 is a constraint of its own.
    constraints = [
        {"a": (1, 2), "op": "<=", "b": 2, "label_at": 0.14},
        {"a": (-1, 1), "op": "<=", "b": 4, "label_at": 0.82},
        {"a": (0, 1), "op": ">=", "b": 0, "label": None},
    ]
    lp_region(
        ax, constraints, xlim, ylim, nonneg=False, colors=[BLUE, GREEN, INK]
    )
    ax.axvline(0, color=MUTED, linewidth=0.8, zorder=1)  # the x2-axis (x1 = 0)
    # Optimal contour x1 + 3 x2 = 4.
    draw_line(
        ax, 1, 3, 4, xlim, ylim, color=MUTED, linestyle="--", linewidth=1.1
    )
    # Corner points (circles) and the two non-corner BFSs (squares).
    _point(ax, (-4, 0), r"$(-4, 0)$", (-4, 6), ha="right")
    _point(ax, (2, 0), r"$(2, 0)$", (4, 6))
    _point(ax, (-2, 2), None, (0, 0), color=ORANGE)
    # Optimum label placed in free space with a leader arrow.
    ax.annotate(
        r"$\mathbf{x}^* = (-2, 2)$",
        xy=(-2, 2),
        xytext=(0.2, 2.7),
        fontsize=10.5,
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8, shrinkB=5),
    )
    _point(ax, (0, 0), r"$(0, 0)$", (5, 5), marker="s", hollow=True)
    # (0, 1): label inside the triangle, clear of the dashed contour
    _point(ax, (0, 1), r"$(0, 1)$", (-6, -4), ha="right", va="top",
           marker="s", hollow=True)  # fmt: skip
    _halo(fig)
    return save_mpl(fig, "08-02-free-variable-region", outdir)


def cut_cube(outdir, lang="en"):
    """Unit cube with the corner (1,1,1) cut off by x1 + x2 + x3 = 5/2.

    Used in ex-08-02-a3-q4 and ex-08-02-a3-q5 (fig-08-02-cut-cube). Extreme points
    A..G are cube corners, H = (1, 1/2, 1), I = (1, 1, 1/2), J = (1/2, 1, 1).
    Viewed from the positive octant, so the three edges at A are hidden
    (dashed).
    """
    fig = plt.figure(figsize=(5.0, 4.6))
    ax = fig.add_subplot(projection="3d")
    V = {
        "A": (0, 0, 0), "B": (1, 0, 0), "C": (0, 1, 0), "D": (0, 0, 1),
        "E": (1, 1, 0), "F": (0, 1, 1), "G": (1, 0, 1),
        "H": (1, 0.5, 1), "I": (1, 1, 0.5), "J": (0.5, 1, 1),
    }  # fmt: skip
    edges = ["AB", "AC", "AD", "BE", "BG", "CE", "CF", "DF", "DG",
             "EI", "FJ", "GH", "HI", "IJ", "JH"]  # fmt: skip
    hidden = {"AB", "AC", "AD"}
    for e in edges:
        p, q = V[e[0]], V[e[1]]
        ax.plot(
            *zip(p, q),
            color=MUTED if e in hidden else INK,
            linestyle="--" if e in hidden else "-",
            linewidth=1.4,
        )
    # The cutting triangle H-I-J, lightly filled.
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    tri = Poly3DCollection(
        [[V["H"], V["I"], V["J"]]],
        facecolor=ORANGE,
        alpha=0.25,
        edgecolor=ORANGE,
    )
    ax.add_collection3d(tri)
    # Labels of the crowded points A, H, I, J are placed in screen space:
    # offset (points) from the projected point, chosen to clear the edges.
    screen_labels = {
        "A": (-14, 4, "right"),
        "H": (-6, 13, "right"),
        "I": (10, 4, "left"),
        "J": (8, -16, "left"),
    }
    # Vertex labels, nudged outward from the cube center.
    for name, p in V.items():
        d = [c - 0.5 for c in p]
        lab = [c + 0.13 * di / max(0.5, abs(di) + 1e-9) for c, di in zip(p, d)]
        ax.scatter(*p, color=INK, s=14, depthshade=False)
        if name in screen_labels:
            continue  # placed in screen coordinates below
        ax.text(*lab, name, fontsize=12, ha="center", va="center", color=INK)
    # No box axes: short axis arrows from the cube corners B, C, D instead.
    ax.set_axis_off()
    for end, name in [
        ((1.45, 0, 0), "$x_1$"),
        ((0, 1.45, 0), "$x_2$"),
        ((0, 0, 1.45), "$x_3$"),
    ]:
        start = tuple(min(c, 1) for c in end)
        ax.plot(*zip(start, end), color=MUTED, linewidth=1.0)
        ax.text(
            *[c * 1.08 for c in end],
            name,
            fontsize=12,
            color=MUTED,
            ha="center",
            va="center",
        )
    ax.view_init(elev=22, azim=35)
    ax.set_box_aspect((1, 1, 1))
    from mpl_toolkits.mplot3d import proj3d

    fig.canvas.draw()  # fixes the projection used below
    for name, (dx, dy, ha) in screen_labels.items():
        x2, y2, _ = proj3d.proj_transform(*V[name], ax.get_proj())
        ax.annotate(
            name,
            xy=(x2, y2),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=12,
            ha=ha,
            va="center",
            color=INK,
        )
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis.pane.set_edgecolor("none")
    _halo(fig)
    return save_mpl(fig, "08-02-cut-cube", outdir)


def degenerate_vertex(outdir, lang="en"):
    """Pentagon A B C D E with four lines through the corner C = (8, 2).

    Used in ex-08-02-a3-q12 (fig-08-02-degenerate-vertex). Constraints:
    x1 - 2 x2 <= 4 (edge BC), 4 x1 + x2 <= 34 (edge CD), -2 x1 + 7 x2 <= 28
    (edge ED), and the redundant 2 x1 - x2 <= 14, 4 x1 - 3 x2 <= 26, which
    also pass through C. Dashed: contour 2 x1 + x2 = 20 through D = (7, 6).
    """
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    xlim, ylim = (0, 11), (0, 8)
    constraints = [
        {"a": (1, -2), "op": "<=", "b": 4, "label_at": 0.9},
        {"a": (4, 1), "op": "<=", "b": 34, "label_at": 0.1},
        {"a": (-2, 7), "op": "<=", "b": 28, "label_at": 0.25},
        {"a": (2, -1), "op": "<=", "b": 14, "label_at": 0.58},
        {"a": (4, -3), "op": "<=", "b": 26, "label_at": 0.86},
    ]
    lp_region(
        ax,
        constraints,
        xlim,
        ylim,
        colors=[BLUE, GREEN, PURPLE, PALETTE[4], PALETTE[5]],
        label_size=9.5,
    )
    draw_line(
        ax, 2, 1, 20, xlim, ylim, color=MUTED, linestyle="--", linewidth=1.1
    )
    for name, p, off, ha in [
        ("A", (0, 0), (6, 4), "left"),
        ("B", (4, 0), (2, 5), "left"),
        ("C", (8, 2), (12, -14), "left"),
        ("D", (7, 6), (-8, -20), "right"),  # below the edge ED, left of CD
        ("E", (0, 4), (6, 4), "left"),
    ]:
        _point(
            ax,
            p,
            f"${name}$",
            off,
            ha=ha,
            color=ORANGE if name == "C" else INK,
        )
    _halo(fig)
    return save_mpl(fig, "08-02-degenerate-vertex", outdir)


def temporary_degeneracy(outdir, lang="en"):
    """Region of max 3 x1 + 2 x2 with the degenerate corner (1, 0).

    Used in ex-08-02-a3-q13 (fig-08-02-temporary-degeneracy). Constraints
    4 x1 - x2 <= 4, 4 x1 + 3 x2 <= 6, 4 x1 + x2 <= 4, x >= 0; corners
    (0,0), (1,0), (3/4,1), (0,2). The simplex path (0,0) -> (1,0) -> (3/4,1)
    is drawn with arrows; dashed: contour 3 x1 + 2 x2 = 17/4.
    """
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    xlim, ylim = (0, 2), (0, 2.6)
    constraints = [
        {"a": (4, -1), "op": "<=", "b": 4, "label_at": 0.7},
        {"a": (4, 3), "op": "<=", "b": 6, "label_at": 0.3},
        {"a": (4, 1), "op": "<=", "b": 4, "label_at": 0.12},
    ]
    lp_region(ax, constraints, xlim, ylim, colors=[BLUE, GREEN, PURPLE])
    draw_line(
        ax,
        3,
        2,
        17 / 4,
        xlim,
        ylim,
        color=MUTED,
        linestyle="--",
        linewidth=1.1,
    )
    path = [(0, 0), (1, 0), (0.75, 1)]
    for p, q in zip(path, path[1:]):
        ax.annotate(
            "",
            xy=q,
            xytext=p,
            arrowprops=dict(
                arrowstyle="-|>",
                color=ORANGE,
                lw=2.0,
                shrinkA=5,
                shrinkB=6,
                mutation_scale=13,
            ),
            zorder=4,
        )
    _point(ax, (0, 0), r"$(0, 0)$", (6, 6))
    # The degenerate corner: three lines leave it, so its label sits in the
    # free space to the right, with a leader line.
    _point(ax, (1, 0), None, (0, 0))
    ax.annotate(
        _t(lang, "$(1, 0)$\ndegenerate", "$(1, 0)$\n退化"),
        xy=(1, 0),
        xytext=(1.55, 0.3),
        fontsize=10.5,
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8, shrinkB=5),
    )
    _point(
        ax,
        (0.75, 1),
        r"$\mathbf{x}^* = (\frac{3}{4}, 1)$",
        (8, 4),
        color=ORANGE,
    )
    _point(ax, (0, 2), r"$(0, 2)$", (6, 6))
    _halo(fig)
    return save_mpl(fig, "08-02-temporary-degeneracy", outdir)


def optimal_ray(outdir, lang="en"):
    """Unbounded region with a ray of optimal solutions.

    Used in ex-08-02-a3-q14 (b) (fig-08-02-optimal-ray): max z = -3 x1 + x2,
    2 x1 - x2 <= 2, -3 x1 + x2 <= 3, x >= 0. The objective is constant
    (z = 3) on the boundary ray (t, 3 + 3t), t >= 0, drawn thick.
    """
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    xlim, ylim = (0, 3.5), (0, 10)
    constraints = [
        {"a": (2, -1), "op": "<=", "b": 2, "label_at": 0.75},
        {"a": (-3, 1), "op": "<=", "b": 3, "label": None},
    ]
    lp_region(ax, constraints, xlim, ylim, colors=[BLUE, GREEN])
    # The optimal ray from (0, 3) in direction (1, 3), clipped to the window.
    ax.plot(
        [0, 7 / 3],
        [3, 10],
        color=ORANGE,
        linewidth=4.5,
        zorder=3,
        solid_capstyle="butt",
    )
    # A lower contour z = -3 x1 + x2 = 0 for comparison (dashed).
    draw_line(
        ax,
        -3,
        1,
        0,
        xlim,
        ylim,
        color=MUTED,
        linestyle="--",
        linewidth=1.1,
        label=r"$z = 0$",
        label_at=0.75,
    )
    _point(ax, (0, 3), r"$(0, 3)$", (8, -4), va="top", color=ORANGE)
    _point(ax, (1, 0), r"$(1, 0)$", (6, 6))
    # The ray lies on -3 x1 + x2 = 3, so that line gets this label instead.
    ax.text(
        0.06,
        8.9,
        "$-3x_1 + x_2 = 3$"
        + _t(lang, ":\noptimal ray, $z^* = 3$", "：\n最优射线，$z^* = 3$"),
        fontsize=10.5,
        color=ORANGE,
        ha="left",
        va="center",
    )
    _halo(fig)
    return save_mpl(fig, "08-02-optimal-ray", outdir)


FIGURES = {
    "08-02-free-variable-region": free_variable_region,
    "08-02-cut-cube": cut_cube,
    "08-02-degenerate-vertex": degenerate_vertex,
    "08-02-temporary-degeneracy": temporary_degeneracy,
    "08-02-optimal-ray": optimal_ray,
}
