"""Example figure module: copy it to ``pPPcCC.py`` to start your own.

This file is NOT part of the book: ``make_figures.py`` loads only modules
named ``pPPcCC.py``, and the leading underscore keeps this one out. To draw
its two figures for a look::

    python3 scripts/make_figures.py --example --outdir /tmp/fig-preview/example
    python3 scripts/figures/preview.py /tmp/fig-preview/example/00-example-lp.svg \
        /tmp/fig-preview/example/00-example-network.svg

What a figure module contains:

1. imports from ``.common`` (style is applied on import);
2. one function per figure, ``func(outdir, lang="en")``, with a docstring
   saying what the figure shows and where it is used, comments on the data
   and on each drawing step; it writes ``<name>.svg`` and returns the path;
3. the ``FIGURES`` dict at the end, mapping each name to its function.
   Names start with the chapter number: ``"03-02-feasible-region"``.
"""

import matplotlib.pyplot as plt

from .common import (
    MUTED,
    ORANGE,
    draw_line,
    lp_region,
    network_svg,
    save_mpl,
    svg_text,
    write_svg,
)


def example_lp(outdir, lang="en"):
    """Feasible region of a product-mix LP with its optimal vertex.

    Example only. Shows how ``lp_region`` shades the region of
    max 500 x1 + 450 x2 s.t. 6 x1 + 5 x2 <= 60, 10 x1 + 20 x2 <= 150,
    x1 <= 8, x >= 0, and how ``draw_line`` adds an objective contour.
    """
    fig, ax = plt.subplots(figsize=(5.4, 4.2))

    # Constraints a1*x1 + a2*x2 (op) b. Integers keep the vertices exact.
    constraints = [
        # Tuple form: labelled with its equation "6x1 + 5x2 = 60";
        # the dict form below moves the label along the line.
        {"a": (6, 5), "op": "<=", "b": 60, "label_at": 0.12},
        (10, 20, "<=", 150),
        # Dict form for per-constraint options: own label, label position.
        {
            "a": (1, 0),
            "op": "<=",
            "b": 8,
            "label": "$x_1 = 8$",
            "label_at": 0.9,
        },
    ]
    vertices = lp_region(
        ax,
        constraints,
        xlim=(0, 11),
        ylim=(0, 9),
        mark_vertices=True,
        exact=True,
    )

    # The optimal vertex (45/7, 30/7) has the largest objective value.
    c1, c2 = 500, 450
    x_opt = max(vertices, key=lambda p: c1 * p[0] + c2 * p[1])
    z_opt = c1 * x_opt[0] + c2 * x_opt[1]  # exact: Fraction(36000, 7)

    # Objective contour through the optimum, dashed and muted.
    draw_line(
        ax, c1, c2, z_opt, (0, 11), (0, 9),
        color=MUTED, linestyle="--", label=f"$z = {z_opt}$", label_at=0.88,
    )  # fmt: skip

    # Mark and label the optimal vertex.
    x_star = (float(x_opt[0]), float(x_opt[1]))
    ax.plot(*x_star, marker="o", markersize=7, color=ORANGE, zorder=5)
    ax.annotate(
        rf"$\mathbf{{x}}^* = ({x_opt[0]},\ {x_opt[1]})$",
        xy=x_star,
        xytext=(-14, -22),
        textcoords="offset points",
        ha="right",
    )
    return save_mpl(fig, "00-example-lp", outdir)


def example_network(outdir, lang="en"):
    """A small directed network with arc lengths and a highlighted path.

    Example only. Shows ``network_svg`` with an antiparallel pair (a, b)
    and (b, a), which is curved apart automatically, a highlighted path
    s-a-t, node notes with labels [V, P], and a panel title.
    """
    # Node coordinates in px (origin top left, y down).
    nodes = {"s": (50, 120), "a": (180, 55), "b": (180, 185), "t": (310, 120)}
    # Arcs (tail, head, length).
    arcs = [
        ("s", "a", 4),
        ("s", "b", 2),
        ("a", "b", 1),
        ("b", "a", 3),
        ("a", "t", 5),
        ("b", "t", 8),
    ]
    body = network_svg(
        nodes,
        arcs,
        highlight_arcs=[("s", "a"), ("a", "t")],  # a shortest path s-a-t
        highlight_nodes=["s", "t"],
        node_notes={
            "s": "[0, –]",
            "a": "[4, s]",
            "b": "[2, s]",
            "t": "[9, a]",
        },
        note_pos={"s": "below", "t": "below", "b": "below"},
    )
    # Panel title in the top-left corner.
    body += svg_text(
        10, 14, "Shortest path s → t", anchor="start", size=13, math=False
    )
    return write_svg("00-example-network", body, 360, 230, outdir)


FIGURES = {
    "00-example-lp": example_lp,
    "00-example-network": example_network,
}
