"""Figures of chapter 3.4 (Artificial Variables: Two-Phase and Big-M Methods).

03-04-phase-one-path: the iterates of the big-M / two-phase example in the
(x1, x2)-plane; the first iterates lie outside the feasible set.
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


def phase_one_path(outdir, lang="en"):
    """Path of the simplex method for min 4 x1 + x2 (big-M and two-phase).

    Used in 03-04-two-phase-big-m (fig-03-04-phase-one-path). Constraints:
    3 x1 + x2 = 3 (equation), 4 x1 + 3 x2 >= 6, x1 + 2 x2 <= 4, x >= 0.
    The feasible set is the segment of 3 x1 + x2 = 3 between (2/5, 9/5)
    and (3/5, 6/5). Iterates: (0, 0) and (1, 0) (artificial variables
    positive, outside the feasible set), (3/5, 6/5) (end of Phase I) and
    the optimum (2/5, 9/5).
    """
    fig, ax = plt.subplots(figsize=(5.4, 4.3))
    xlim, ylim = (0, 2.2), (0, 2.4)

    # Region of the two inequality constraints (with x >= 0), lightly shaded:
    # 4 x1 + 3 x2 >= 6 and x1 + 2 x2 <= 4. It is the triangle with vertices
    # (3/2, 0), (4, 0) and (0, 2) (the two lines meet at (0, 2)); the axes
    # limits clip it to the window.
    xs = [1.5, 4.0, 0.0]
    ys = [0.0, 0.0, 2.0]
    ax.fill(xs, ys, color=BLUE, alpha=0.10, lw=0, zorder=0)

    # The three constraint lines.
    draw_line(ax, 3, 1, 3, xlim, ylim, color=BLUE, linewidth=1.4)
    # Label of the equality line set apart (top left) to keep the optimum clear.
    ax.text(
        0.03, 2.22, r"$3x_1 + x_2 = 3$", color=BLUE, ha="left", va="center"
    )
    draw_line(
        ax,
        4,
        3,
        6,
        xlim,
        ylim,
        color=GREEN,
        linewidth=1.4,
        label=r"$4x_1 + 3x_2 = 6$",
        label_at=0.86,
    )
    draw_line(
        ax,
        1,
        2,
        4,
        xlim,
        ylim,
        color=PALETTE[3],
        linewidth=1.4,
        label=r"$x_1 + 2x_2 = 4$",
        label_at=0.8,
    )

    # Feasible set: thick segment of 3 x1 + x2 = 3.
    ax.plot(
        [0.4, 0.6],
        [1.8, 1.2],
        color=INK,
        linewidth=5,
        solid_capstyle="round",
        zorder=3,
    )

    # Optimal contour 4 x1 + x2 = 17/5 (dashed).
    draw_line(
        ax,
        4,
        1,
        17 / 5,
        xlim,
        ylim,
        color=MUTED,
        linestyle="--",
        linewidth=1.1,
    )

    # Iterates in order; the first two are infeasible for the original LP.
    path = [(0, 0), (1, 0), (0.6, 1.2), (0.4, 1.8)]
    for p, q in zip(path, path[1:]):
        ax.annotate(
            "", xy=q, xytext=p,
            arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2.0,
                            shrinkA=5, shrinkB=6, mutation_scale=13),
            zorder=4,
        )  # fmt: skip
    for k, p in enumerate(path):
        hollow = k < 2
        ax.plot(
            *p,
            marker="o",
            markersize=7,
            markerfacecolor="none" if hollow else (ORANGE if k == 3 else INK),
            markeredgecolor=INK if k < 3 else ORANGE,
            markeredgewidth=1.4,
            zorder=5,
            clip_on=False
        )

    # Labels of the iterates.
    ax.annotate(
        r"$(0, 0)$",
        xy=(0, 0),
        xytext=(6, 6),
        textcoords="offset points",
        fontsize=10.5,
    )
    ax.annotate(
        r"$(1, 0)$",
        xy=(1, 0),
        xytext=(4, 6),
        textcoords="offset points",
        fontsize=10.5,
    )
    ax.annotate(
        _t(
            lang,
            r"$(\frac{3}{5}, \frac{6}{5})$ end of Phase I",
            r"$(\frac{3}{5}, \frac{6}{5})$ 第一阶段结束",
        ),
        xy=(0.6, 1.2),
        xytext=(10, -4),
        textcoords="offset points",
        fontsize=10.5,
        va="top",
    )
    ax.annotate(
        r"$\mathbf{x}^* = (\frac{2}{5}, \frac{9}{5})$",
        xy=(0.4, 1.8),
        xytext=(14, 10),
        textcoords="offset points",
        fontsize=10.5,
    )

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    _halo(fig)
    return save_mpl(fig, "03-04-phase-one-path", outdir)


FIGURES = {
    "03-04-phase-one-path": phase_one_path,
}
