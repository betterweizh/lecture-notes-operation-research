"""Figures of chapter 4.4, Minimum-Cost Flow.

Hand-written SVG built with ``network_svg``:

* ``04-04-supply-chain``: the six-node supply chain, (a) data (net supplies,
  unit costs, bounds) and (b) the optimal flows of the repaired instance
  with the spanning-tree (basic) arcs highlighted;
* ``04-04-transportation``: the bipartite transportation network with unit
  costs, net supplies and the arcs of the optimal plan highlighted.

Labels are numbers, node names and symbols only (language neutral).
"""

from .common import network_svg, svg_text, write_svg

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

#: Supply chain: node positions (px) and net supplies b_i (units).
SC_NODES = {
    "1": (40, 52),
    "2": (40, 222),
    "3": (148, 150),
    "4": (270, 52),
    "5": (270, 230),
    "6": (360, 140),
}
SC_SUPPLY = {
    "1": "[100]",
    "2": "[200]",
    "3": "[50]",
    "4": "[−150]",
    "5": "[−80]",
    "6": "[−120]",
}
#: Arcs (tail, head, unit cost $, bounds or None for [0, inf)).
SC_ARCS = [
    ("1", "2", 3, None),
    ("1", "3", 4, None),
    ("1", "4", 1, (50, 80)),
    ("2", "3", 5, None),
    ("2", "5", 6, None),
    ("3", "4", 1, (70, 120)),
    ("3", "5", 2, None),
    ("4", "6", 2, (100, 120)),
    ("5", "6", 4, None),
]
#: Optimal flows of the repaired instance (0 <= x46 <= 120), z* = $1 810
#: (verify_ch17.py), and the five arcs strictly between their bounds.
SC_FLOW = {("1", "2"): 0, ("1", "3"): 20, ("1", "4"): 80, ("2", "3"): 50, ("2", "5"): 150,
           ("3", "4"): 120, ("3", "5"): 0, ("4", "6"): 50, ("5", "6"): 70}  # fmt: skip
SC_TREE = [("1", "3"), ("2", "3"), ("2", "5"), ("4", "6"), ("5", "6")]

# Note positions: supplies of the left nodes on the left, others outside.
SC_NOTE_POS = {"1": "above", "2": "below", "3": "below", "4": "above", "5": "below",
               "6": "right"}  # fmt: skip

#: Transportation example: plants A, B (supplies in t) and markets 1-3.
TR_NODES = {
    "A": (60, 60),
    "B": (60, 190),
    "1": (300, 30),
    "2": (300, 125),
    "3": (300, 220),
}
TR_SUPPLY = {
    "A": "[40]",
    "B": "[50]",
    "1": "[−20]",
    "2": "[−25]",
    "3": "[−45]",
}
TR_COST = {
    ("A", "1"): 4,
    ("A", "2"): 6,
    ("A", "3"): 9,
    ("B", "1"): 2,
    ("B", "2"): 7,
    ("B", "3"): 5,
}
#: Optimal plan (cost $445): A1 15, A2 25, B1 5, B3 45.
TR_USED = [("A", "1"), ("A", "2"), ("B", "1"), ("B", "3")]


def supply_chain(outdir, lang="en"):
    """Supply chain of Example 4.4.1: data (a) and optimal flows (b).

    Panel (a): label "c" or "c (l, u)" on every arc, [b_i] next to nodes.
    Panel (b): optimal flow x*_ij on every arc; highlighted arcs are basic
    (strictly between their bounds) and form a spanning tree.
    """
    arcs_a = [
        (i, j, f"{c} ({b[0]}, {b[1]})" if b else c) for i, j, c, b in SC_ARCS
    ]
    pos = {
        ("1", "4"): 0.5,
        ("3", "4"): 0.5,
        ("4", "6"): 0.45,
        ("2", "3"): 0.45,
    }
    side = {("3", "4"): -1, ("1", "3"): -1, ("3", "5"): -1, ("4", "6"): 1}
    body = network_svg(
        SC_NODES, arcs_a, node_notes=SC_SUPPLY, note_pos=SC_NOTE_POS,
        label_pos=pos, label_side=side, offset=(14, 40),
    )  # fmt: skip
    body += svg_text(10, 14, "(a)", anchor="start", math=False)
    arcs_b = [(i, j, SC_FLOW[(i, j)]) for i, j, c, b in SC_ARCS]
    body += network_svg(
        SC_NODES, arcs_b, node_notes=SC_SUPPLY, note_pos=SC_NOTE_POS,
        highlight_arcs=SC_TREE, label_pos=pos, label_side=side, offset=(474, 40),
    )  # fmt: skip
    body += svg_text(470, 14, "(b)", anchor="start", math=False)
    return write_svg("04-04-supply-chain", body, 920, 310, outdir)


def transportation(outdir, lang="en"):
    """Transportation network of Example 4.4.2 with the optimal plan.

    Arc labels are unit costs ($/t); [b_i] next to the nodes; the four arcs
    that carry flow in the optimal plan are highlighted.
    """
    arcs = [(p, m, c) for (p, m), c in TR_COST.items()]
    body = network_svg(
        TR_NODES,
        arcs,
        node_notes=TR_SUPPLY,
        note_pos={"A": "left", "B": "left", "1": "right", "2": "right", "3": "right"},
        highlight_arcs=TR_USED,
        # Labels spread along the arcs so that the crossing arcs stay clear.
        label_pos={("A", "1"): 0.5, ("A", "2"): 0.75, ("A", "3"): 0.22,
                   ("B", "1"): 0.22, ("B", "2"): 0.75, ("B", "3"): 0.5},
        label_side={("A", "3"): -1, ("B", "3"): -1},
        offset=(50, 14),
    )  # fmt: skip
    return write_svg("04-04-transportation", body, 420, 266, outdir)


FIGURES = {
    "04-04-supply-chain": supply_chain,
    "04-04-transportation": transportation,
}
