"""Figures of chapter 3.7 (Sensitivity Analysis).

* ``03-07-rhs-ranging``      right-hand-side ranging of the resistor constraint
  of the product-mix LP (max 3x1 + 4x2), section "A geometric view";
* ``03-07-cost-ranging``     objective contours rotating about the optimal corner
  of the same LP, section "A geometric view";
* ``03-07-parametric-value`` optimal value z*(b1) of the running example,
  section "Parametric right-hand side".

All numbers were recomputed in exact arithmetic (scratch scripts of the
chapter) and checked with scipy's HiGHS solver.
"""

from fractions import Fraction as F

import matplotlib.pyplot as plt

from .common import (
    BLUE,
    GREEN,
    GRID,
    INK,
    MUTED,
    ORANGE,
    PALETTE,
    draw_line,
    lp_region,
    save_mpl,
)

# Product-mix LP used by both geometric figures:
#   max 3 x1 + 4 x2
#   s.t. 2 x1 + 3 x2 <= 1200 (resistors), 2 x1 + x2 <= 1000 (capacitors),
#        4 x2 <= 800 (chips), x >= 0.
XLIM = (0, 640)
YLIM = (0, 420)


def _region(ax, resistor_label_at=0.25, muted=False):
    """Shade the region with b1 = 1200 and label the three constraint lines.

    With ``muted=True`` the constraint lines are drawn thin and gray, so that
    objective contours drawn on top stand out.
    """
    cols = (MUTED, MUTED, MUTED) if muted else (BLUE, ORANGE, GREEN)
    # In the muted version the two binding lines carry no label: the extreme
    # objective contours drawn on top of them are labeled instead.
    constraints = [
        {"a": (2, 3), "op": "<=", "b": 1200,
         "label": None if muted else r"$2x_1+3x_2=1200$",
         "label_at": resistor_label_at, "color": cols[0]},
        {"a": (2, 1), "op": "<=", "b": 1000,
         "label": None if muted else r"$2x_1+x_2=1000$",
         "label_at": 0.16, "color": cols[1]},
        {"a": (0, 4), "op": "<=", "b": 800, "label": r"$4x_2=800$",
         "label_at": 0.12, "color": cols[2]},  # clear of the x2 axis
    ]  # fmt: skip
    return lp_region(ax, constraints, XLIM, YLIM, exact=True)


def rhs_ranging(outdir, lang="en"):
    """Resistor line 2x1 + 3x2 = b1 for b1 = 1000, 1200, 1400.

    The optimal corner (intersection with the capacitor line) slides from
    (500, 0) to (400, 200); the basis stays the same over [1000, 1400].
    Used in chapter 3.7, section "A geometric view".
    """
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    _region(ax, resistor_label_at=0.25)

    # The two end points of the resistor range, dashed and muted.
    for b1, at in ((1000, 0.45), (1400, 0.10)):
        draw_line(
            ax, 2, 3, b1, XLIM, YLIM, color=MUTED, linestyle="--", linewidth=1.2,
            label=rf"$b_1={b1}$", label_at=at,
        )  # fmt: skip

    # Optimal corners for b1 = 1000, 1200, 1400 (solve with the capacitor line):
    # 2x1 + x2 = 1000 and 2x1 + 3x2 = b1  =>  x2 = (b1 - 1000)/2, x1 = 500 - x2/2.
    corners = [(500, 0), (450, 100), (400, 200)]
    for (x1, x2), col in zip(corners, (MUTED, INK, MUTED)):
        ax.plot(x1, x2, "o", color=col, markersize=6, zorder=5, clip_on=False)
    # Arrow along the capacitor line showing the movement of the optimum.
    # (drawn a little to the right of the line so that the line stays visible)
    ax.annotate(
        "", xy=(426, 186), xytext=(510, 18),
        arrowprops=dict(arrowstyle="->", color=INK, lw=1.3, shrinkA=6, shrinkB=6),
    )  # fmt: skip
    ax.annotate(
        r"$\mathbf{x}^*$",
        xy=(450, 100),
        xytext=(408, 78),
        fontsize=11,
        color=INK,
    )
    ax.annotate(
        r"$(500,0)$", xy=(500, 0), xytext=(508, 12), fontsize=10, color=MUTED
    )
    ax.annotate(
        r"$(400,200)$",
        xy=(400, 200),
        xytext=(410, 214),
        fontsize=10,
        color=MUTED,
    )
    return save_mpl(fig, "03-07-rhs-ranging", outdir)


def cost_ranging(outdir, lang="en"):
    """Objective contours c1 x1 + 4 x2 = z through the optimal corner (450, 100).

    c1 = 8/3 (the contour coincides with the resistor line), c1 = 3 (current)
    and c1 = 8 (it coincides with the capacitor line): the corner stays optimal for
    8/3 <= c1 <= 8. Used in chapter 3.7, section "A geometric view".
    """
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    _region(ax, resistor_label_at=0.30, muted=True)

    # Contours through (450, 100): z = c1*450 + 4*100. The label of c1 = 8/3
    # goes to the lower-right end, the others to the upper-left end, where the
    # three lines are far enough apart.
    contours = [
        (F(8, 3), r"$c_1=8/3$", BLUE, 0.86),
        (F(3), r"$c_1=3$", INK, 0.12),
        (F(8), r"$c_1=8$", ORANGE, 0.10),
    ]
    for c1, label, col, at in contours:
        z = c1 * 450 + 4 * 100
        draw_line(
            ax, c1, 4, z, XLIM, YLIM, color=col, linestyle=(0, (5, 3)), linewidth=1.8,
            label=label, label_at=at,
        )  # fmt: skip
    ax.plot(450, 100, "o", color=INK, markersize=6, zorder=5)
    ax.annotate(
        r"$\mathbf{x}^*=(450,100)$",
        xy=(450, 100),
        xytext=(285, 62),
        fontsize=10,
    )
    return save_mpl(fig, "03-07-cost-ranging", outdir)


def parametric_value(outdir, lang="en"):
    """Optimal value z*(b1) of the running example of chapter 3.7.

    Pieces (basis, interval, z*):
      (x3, x5, x6)  b1 <= 0        10 - 4/3 b1
      (x1, x3, x5)  0 <= b1 <= 7   10 - b1
      (x1, x2, x3)  7 <= b1 <= 13  2 b1 - 11
      (x1, x2, x5)  13 <= b1 <= 40/3  5 b1 - 50
      infeasible for b1 > 40/3.
    Used in chapter 3.7, section "Parametric right-hand side".
    """
    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    pieces = [
        (
            -4.0,
            0.0,
            lambda t: 10 - 4 * t / 3,
            r"$(x_3,x_5,x_6)$",
            r"$\pi_1=-4/3$",
        ),
        (0.0, 7.0, lambda t: 10 - t, r"$(x_1,x_3,x_5)$", r"$\pi_1=-1$"),
        (7.0, 13.0, lambda t: 2 * t - 11, r"$(x_1,x_2,x_3)$", r"$\pi_1=2$"),
        (13.0, 40 / 3, lambda t: 5 * t - 50, r"$(x_1,x_2,x_5)$", r"$\pi_1=5$"),
    ]
    # Label offsets (points) for each piece: basis label and slope label.
    # (the label of the falling piece 2 sits up and to the right, clear of it)
    offsets = [(0, -62), (30, 8), (-30, 8), (-62, 6)]
    for k, ((lo, hi, f, basis, slope), (dx, dy)) in enumerate(
        zip(pieces, offsets)
    ):
        ax.plot([lo, hi], [f(lo), f(hi)], color=PALETTE[k], linewidth=2.2)
        mid = (lo + hi) / 2
        ax.annotate(
            basis + "\n" + slope, xy=(mid, f(mid)), xytext=(dx, dy),
            textcoords="offset points", ha="center", va="bottom", fontsize=9,
            color=PALETTE[k],
        )  # fmt: skip
    # Breakpoints.
    for t in (0, 7, 13, 40 / 3):
        z = min(p[2](t) for p in pieces if p[0] <= t <= p[1] + 1e-12)
        ax.plot(t, z, "o", color=INK, markersize=3.5, zorder=5)
    # Current right-hand side b1 = 11, z* = 11.
    ax.plot(
        11, 11, "o", color=INK, markersize=7, zorder=6, mfc="none", mew=1.5
    )
    # Infeasible region b1 > 40/3.
    ax.axvspan(40 / 3, 16, color=MUTED, alpha=0.18, lw=0)
    ax.text(
        14.7,
        4,
        "infeasible" if lang == "en" else "不可行",
        rotation=90,
        ha="center",
        va="center",
        color=MUTED,
        fontsize=10,
    )
    ax.set_xlim(-4, 16)
    ax.set_ylim(0, 19)
    ax.set_xticks([-4, 0, 4, 7, 11, 13, 16])
    ax.set_xlabel(r"$b_1$")
    ax.set_ylabel(r"$z^*(b_1)$")
    ax.grid(True, color=GRID, linewidth=0.6)
    return save_mpl(fig, "03-07-parametric-value", outdir)


FIGURES = {
    "03-07-rhs-ranging": rhs_ranging,
    "03-07-cost-ranging": cost_ranging,
    "03-07-parametric-value": parametric_value,
}
