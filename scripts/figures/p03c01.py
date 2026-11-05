"""Figures of chapter 3.1, Linear Programming: Model Building.

* ``03-01-inventory-flow``: flow view of the inventory-balance equations
  I_t = I_{t-1} + p_t - d_t for three periods (hand-written SVG);
* ``03-01-piecewise-profit``: the concave piecewise-linear contribution of
  juice glasses (slope 500 up to x1 = 3, slope 300 beyond) used in the
  piecewise-linear product-mix example (matplotlib).

Labels are the book's symbols and numbers (language neutral).
"""

import matplotlib.pyplot as plt

from .common import (
    BLUE,
    INK,
    MUTED,
    ORANGE,
    network_svg,
    save_mpl,
    svg_line,
    svg_text,
    write_svg,
)


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


def inventory_flow(outdir, lang="en"):
    """Three period nodes linked by inventory arcs; production in, demand out.

    Used in Sec. "Multi-period models" of chapter 3.1 to show that the
    balance equation (eq-03-01-inventory-balance) is "inflow = outflow" at the
    node of each period.
    """
    # Period nodes on one horizontal line (px, y down).
    y = 135
    xs = {"1": 150, "2": 290, "3": 430}
    nodes = {k: (x, y) for k, x in xs.items()}
    # Inventory arcs between consecutive periods, labeled I_1, I_2.
    arcs = [("1", "2", "I_1"), ("2", "3", "I_2")]
    body = network_svg(nodes, arcs, node_radius=20, highlight_arcs=arcs)

    r = 20  # node radius, used to stop the extra arrows at the circles
    # Initial inventory I_0 enters period 1 from the left.
    body += svg_line(
        40, y, xs["1"], y, arrow=True, shorten_end=r, color=BLUE, width=2.6
    )
    body += svg_text(78, y - 14, "I_0")
    # Final inventory I_3 leaves period 3 to the right.
    body += svg_line(
        xs["3"], y, 540, y, arrow=True, shorten_start=r, color=BLUE, width=2.6
    )
    body += svg_text(500, y - 14, "I_3")
    for k, x in xs.items():
        # Production p_t enters the node from above.
        body += svg_line(
            x, 40, x, y, arrow=True, shorten_end=r, color=ORANGE, width=2
        )
        body += svg_text(x + 18, 64, f"p_{k}", anchor="start")
        # Demand d_t leaves the node downwards.
        body += svg_line(
            x, y, x, 232, arrow=True, shorten_start=r, color=INK, width=2
        )
        body += svg_text(x + 18, 208, f"d_{k}", anchor="start")
    return write_svg("03-01-inventory-flow", body, 560, 250, outdir)


def piecewise_profit(outdir, lang="en"):
    """Concave piecewise-linear contribution r(x1) with a breakpoint at x1 = 3.

    Used in Sec. "Piecewise-linear profits and costs" of chapter 3.1
    (ex-03-01-piecewise-profit). The dashed line is the constant-rate
    contribution 500 x1 of the original product mix, for comparison.
    """
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    x_max = 8.0  # demand limit of juice glasses (hundreds of cases/week)

    # r(x1) = 500 x1 on [0, 3], then 1500 + 300 (x1 - 3): two segments.
    ax.plot([0, 3], [0, 1500], color=BLUE, linewidth=2.4, zorder=3)
    ax.plot(
        [3, x_max],
        [1500, 1500 + 300 * (x_max - 3)],
        color=BLUE,
        linewidth=2.4,
        zorder=3,
    )
    # Constant-rate contribution 500 x1, dashed and muted.
    ax.plot(
        [0, x_max],
        [0, 500 * x_max],
        color=MUTED,
        linestyle="--",
        linewidth=1.4,
        zorder=2,
    )
    ax.text(
        6.3, 500 * 6.3 + 120, "$500x_1$", color=MUTED, ha="right", va="bottom"
    )

    # Breakpoint (3, 1500) and a guide line down to the x1 axis.
    ax.plot([3], [1500], marker="o", color=ORANGE, markersize=6, zorder=4)
    ax.plot([3, 3], [0, 1500], color=MUTED, linestyle=":", linewidth=1.2)
    ax.annotate(
        "$(3,\\ 1500)$",
        xy=(3, 1500),
        xytext=(10, -16),
        textcoords="offset points",
    )

    # Slope labels beside the two segments (the only words: translated).
    slope = {"zh": "斜率"}.get(lang, "slope")
    ax.text(0.9, 1250, f"{slope} 500", color=BLUE, ha="center")
    ax.text(6.0, 1750, f"{slope} 300", color=BLUE, ha="center")

    ax.set_xlim(0, x_max + 0.2)
    ax.set_ylim(0, 4200)
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$r(x_1)$")
    _halo(fig)
    return save_mpl(fig, "03-01-piecewise-profit", outdir)


FIGURES = {
    "03-01-inventory-flow": inventory_flow,
    "03-01-piecewise-profit": piecewise_profit,
}
