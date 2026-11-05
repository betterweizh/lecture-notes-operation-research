"""Figures of chapter 4.2, Shortest-Path Problems.

All figures are hand-written SVG built with ``network_svg``:

* ``04-02-dijkstra-network``: the seven-node directed network with arc lengths
  (Dijkstra example);
* ``04-02-dijkstra-tree``: the same network with the final labels [V(j), P(j)]
  and the shortest-path tree from node 1 highlighted;
* ``04-02-negative-arc``: the three-node counterexample showing that Dijkstra's
  algorithm can fail with a negative arc length;
* ``04-02-floyd-network``: the four-node directed network of the
  Floyd-Warshall example (one negative arc, no negative cycle).

Labels are numbers, node names and the book's symbols (language neutral).
"""

from .common import network_svg, write_svg

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

#: Seven-node network: positions in px (scaled from the original drawing,
#: origin top left, y down) and arcs (tail, head, length).
SP_NODES = {
    "1": (40, 140),
    "2": (168, 58),
    "3": (206, 140),
    "4": (124, 230),
    "5": (372, 80),
    "6": (350, 230),
    "7": (500, 148),
}
SP_ARCS = [
    ("1", "2", 2),
    ("1", "4", 4),
    ("2", "5", 5),
    ("3", "1", 10),
    ("3", "2", 11),
    ("3", "6", 7),
    ("3", "7", 5),
    ("4", "3", 3),
    ("4", "6", 1),
    ("5", "3", 8),
    ("5", "7", 6),
    ("6", "7", 9),
]
# Label placement: 3 -> 7 and 5 -> 3 cross near their midpoints, 3 -> 1 and
# 3 -> 2 meet at node 3.
SP_LABEL_POS = {
    ("3", "7"): 0.62,
    ("5", "3"): 0.3,
    ("3", "2"): 0.45,
    ("3", "1"): 0.5,
}
SP_LABEL_SIDE = {
    ("3", "1"): 1,
    ("3", "2"): -1,
    ("4", "3"): -1,
    ("4", "6"): -1,
    ("3", "6"): 1,
}

#: Final Dijkstra labels [V(j), P(j)] from node 1 and the tree arcs (P(j), j).
SP_LABELS = {
    "1": "[0, –]",
    "2": "[2, 1]",
    "3": "[7, 4]",
    "4": "[4, 1]",
    "5": "[7, 2]",
    "6": "[5, 4]",
    "7": "[12, 3]",
}
SP_TREE = [
    ("1", "2"),
    ("1", "4"),
    ("2", "5"),
    ("4", "3"),
    ("4", "6"),
    ("3", "7"),
]


def dijkstra_network(outdir, lang="en"):
    """The seven-node directed network with arc lengths.

    Used in Example "Shortest route from node 1 to node 7" (sec-04-02-problem).
    """
    body = network_svg(
        SP_NODES,
        SP_ARCS,
        label_pos=SP_LABEL_POS,
        label_side=SP_LABEL_SIDE,
        offset=(0, 0),
    )
    return write_svg("04-02-dijkstra-network", body, 540, 260, outdir)


def dijkstra_tree(outdir, lang="en"):
    """Final labels [V(j), P(j)] and the shortest-path tree from node 1.

    Used in sec-04-02-dijkstra after the iteration table. Tree arcs (P(j), j)
    are thick; nodes on the shortest path 1 -> 4 -> 3 -> 7 are filled. Each
    node carries its permanent label; labels sit where no arc passes.
    """
    body = network_svg(
        SP_NODES,
        SP_ARCS,
        label_pos=SP_LABEL_POS,
        label_side=SP_LABEL_SIDE,
        highlight_arcs=SP_TREE,
        highlight_nodes=["1", "4", "3", "7"],
        node_notes=SP_LABELS,
        note_pos={
            "1": "left",
            "2": "above",
            "3": (2, 34),
            "4": "below",
            "5": "above",
            "6": "below",
            "7": "below",
        },
        offset=(48, 6),
    )
    return write_svg("04-02-dijkstra-tree", body, 600, 290, outdir)


def negative_arc(outdir, lang="en"):
    """Three-node counterexample: Dijkstra with a negative arc length.

    Used in Example "A negative arc defeats Dijkstra's algorithm"
    (sec-04-02-dijkstra-correctness). Arcs (1,2) = 2, (1,3) = 3, (3,2) = -2.
    Node notes show the labels Dijkstra makes permanent ([2, 1] at node 2 is
    wrong); the true shortest path 1 -> 3 -> 2 of length 1 is highlighted.
    """
    nodes = {"1": (40, 110), "2": (230, 40), "3": (230, 180)}
    arcs = [("1", "2", 2), ("1", "3", 3), ("3", "2", "−2")]
    body = network_svg(
        nodes,
        arcs,
        highlight_arcs=[("1", "3"), ("3", "2")],
        node_notes={"1": "[0, –]", "2": "[2, 1]", "3": "[3, 1]"},
        note_pos={"1": "left", "2": "right", "3": "right"},
        label_side={("3", "2"): -1},
        offset=(60, 10),
    )
    return write_svg("04-02-negative-arc", body, 360, 230, outdir)


def floyd_network(outdir, lang="en"):
    """Four-node directed network of the Floyd-Warshall example.

    Used in Example "Floyd-Warshall on a four-node network" (sec-04-02-floyd).
    Arcs: (1,2) = 5, (1,4) = -3, (2,3) = 3, (3,1) = 3, (3,4) = 7,
    (4,3) = 2. The pair (3,4), (4,3) is curved apart automatically.
    """
    nodes = {"1": (40, 40), "2": (250, 40), "3": (250, 200), "4": (40, 200)}
    arcs = [
        ("1", "2", 5),
        ("1", "4", "−3"),
        ("2", "3", 3),
        ("3", "1", 3),
        ("3", "4", 7),
        ("4", "3", 2),
    ]
    body = network_svg(
        nodes, arcs, label_side={("1", "4"): -1}, offset=(20, 10)
    )
    return write_svg("04-02-floyd-network", body, 330, 260, outdir)


FIGURES = {
    "04-02-dijkstra-network": dijkstra_network,
    "04-02-dijkstra-tree": dijkstra_tree,
    "04-02-negative-arc": negative_arc,
    "04-02-floyd-network": floyd_network,
}
