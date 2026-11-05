"""Figures of chapter 5.4, Shortest Paths by Dynamic Programming.

Hand-written SVG built with ``network_svg``:

* ``05-04-acyclic-network``: a six-node acyclic network with letter names,
  the topological number of each node (small gray number), the labels
  [V(j), P(j)] of the forward recursion and the shortest-path tree from the
  source D highlighted;
* ``05-04-negative-arc``: a four-node network with a negative arc
  (3 -> 2 of length -4) and the antiparallel arc 2 -> 3, used for the
  successive-approximation procedures and for Floyd's and Dantzig's
  all-pairs recursions; the shortest-path tree from node 1 is highlighted.

Labels are node names, numbers and symbols only (language neutral).
"""

from .common import MUTED, network_svg, svg_text, write_svg

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

#: Acyclic example: positions in px and arcs (tail, head, length).
ACYCLIC_NODES = {
    "D": (45, 130),
    "C": (185, 50),
    "A": (185, 215),
    "E": (335, 130),
    "F": (480, 215),
    "B": (500, 50),
}
ACYCLIC_ARCS = [
    ("D", "A", 4), ("D", "C", 2), ("C", "A", 1), ("C", "E", 6),
    ("A", "E", 2), ("A", "F", 5), ("E", "F", 2), ("E", "B", 6),
    ("F", "B", 1),
]  # fmt: skip
#: Topological numbering found by removing source nodes one at a time.
ACYCLIC_ORDER = {"D": 1, "C": 2, "A": 3, "E": 4, "F": 5, "B": 6}
#: Forward labels [V(j), P(j)] (verified with networkx, sp_verify.py).
ACYCLIC_LABELS = {
    "D": "[0, –]", "C": "[2, D]", "A": "[3, C]",
    "E": "[5, A]", "F": "[7, E]", "B": "[8, F]",
}  # fmt: skip
ACYCLIC_TREE = [("D", "C"), ("C", "A"), ("A", "E"), ("E", "F"), ("F", "B")]

#: Four-node network with a negative arc.
NEG_NODES = {"1": (45, 120), "2": (200, 45), "3": (200, 195), "4": (355, 120)}
NEG_ARCS = [
    ("1", "2", 2), ("1", "3", 5), ("2", "3", 6), ("3", "2", "−4"),
    ("2", "4", 3), ("3", "4", 1),
]  # fmt: skip
#: Shortest-path tree from node 1: 1-3, 3-2, 2-4 (V = 0, 1, 5, 4).
NEG_TREE = [("1", "3"), ("3", "2"), ("2", "4")]
NEG_LABELS = {"1": "[0, –]", "2": "[1, 3]", "3": "[5, 1]", "4": "[4, 2]"}


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def acyclic_network(outdir, lang="en"):
    """Six-node acyclic network with topological numbers and labels.

    Used in chapter 5.4 (topological ordering and the one-pass recursion).
    The gray number at the upper left of each node is its position in the
    topological order; the label [V(j), P(j)] beside it comes from the
    forward recursion; highlighted arcs form the shortest-path tree from D.
    """
    nodes = ACYCLIC_NODES
    notes = dict(ACYCLIC_LABELS)
    note_pos = {"D": "below", "C": "above", "A": "below", "E": (0, 30),
                "F": "below", "B": "above"}  # fmt: skip
    body = network_svg(
        nodes,
        ACYCLIC_ARCS,
        highlight_arcs=ACYCLIC_TREE,
        highlight_nodes=["D"],
        node_notes=notes,
        note_pos=note_pos,
        label_side={("C", "A"): -1, ("A", "E"): -1, ("E", "F"): 1},
    )
    # Topological numbers: small gray digits next to each node, placed where
    # no arc or label is (offsets in px from the node center).
    num_pos = {"D": (-22, -20), "C": (-30, -8), "A": (-31, 4),
               "E": (0, -28), "F": (26, -16), "B": (28, 0)}  # fmt: skip
    for n, (x, y) in nodes.items():
        dx, dy = num_pos[n]
        body += svg_text(
            x + dx, y + dy, str(ACYCLIC_ORDER[n]), size=12, color=MUTED
        )
    return write_svg("05-04-acyclic-network", body, 545, 260, outdir)


def negative_arc(outdir, lang="en"):
    """Four-node network with the negative arc 3 -> 2 (length -4).

    Used in chapter 5.4 for the successive-approximation procedures and the
    all-pairs recursions, and in chapter 8.6 (procedures 2 and 3). Labels
    [V(j), P(j)] are the final shortest distances from node 1 and their
    predecessors; the shortest-path tree 1-3-2-4 is highlighted.
    """
    body = network_svg(
        NEG_NODES,
        NEG_ARCS,
        highlight_arcs=NEG_TREE,
        highlight_nodes=["1"],
        node_notes=NEG_LABELS,
        note_pos={"1": "below", "2": "above", "3": "below", "4": "below"},
        antiparallel_bend=0.22,
    )
    return write_svg("05-04-negative-arc", body, 400, 245, outdir)


FIGURES = {
    "05-04-acyclic-network": acyclic_network,
    "05-04-negative-arc": negative_arc,
}
