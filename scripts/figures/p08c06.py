"""Figures of chapter 8.6, DP Fundamentals and Shortest Paths (Part VIII,
Examples and Problems).

Hand-written SVG built with ``network_svg``; the networks are redrawn from
their data (Dreyfus and Law, 1977, Ch. 4, and a staged route network):

* ``08-06-acyclic-network``: six-node acyclic network with unnumbered nodes
  (letters a-f) for the topological-ordering problems;
* ``08-06-dijkstra-network``: six-node network with a cycle; the shortest-path
  tree from node 1 is highlighted and every node carries [V(j), P(j)];
* ``08-06-all-pairs-network``: six-node network with undirected edges and two
  pairs of opposite arcs with negative lengths (all-pairs problems);
* ``08-06-turn-grid``: nine-node diamond grid for the turn-penalty problem,
  with the optimal path 1-3-5-7-9 highlighted;
* ``08-06-staged-route``: (a) the seven-city route network with the shortest
  route 1-2-3-5-7 highlighted; (b) the same network drawn in five stages,
  with a zero-length "stay" arc at city 7 (dashed).

Labels are node names, numbers and panel letters only (language neutral).
"""

from .common import GRID, MUTED, network_svg, svg_text, write_svg

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

#: Acyclic network: letters a-f in px positions (pentagon with a center).
ACYC_NODES = {
    "e": (40, 140),  # left
    "a": (190, 40),  # top
    "d": (190, 140),  # center
    "b": (340, 140),  # right
    "f": (110, 245),  # bottom left
    "c": (270, 245),  # bottom right
}
ACYC_ARCS = [
    ("e", "a", 4), ("a", "d", 2), ("a", "b", 6), ("e", "d", 3), ("d", "b", 7),
    ("e", "f", 5), ("d", "f", 1), ("d", "c", 5), ("f", "c", 2), ("c", "b", 3),
]  # fmt: skip

#: Network with a cycle for Dijkstra's method; lengths >= 0.
DIJ_NODES = {"1": (40, 130), "2": (170, 45), "3": (170, 215),
             "4": (400, 45), "5": (400, 215), "6": (520, 130)}  # fmt: skip
DIJ_ARCS = [
    ("1", "2", 5), ("1", "3", 1), ("2", "3", 4), ("2", "4", 6), ("3", "4", 7),
    ("5", "2", 1), ("3", "5", 2), ("5", "4", 3), ("4", "6", 2), ("5", "6", 8),
]  # fmt: skip
DIJ_TREE = [("1", "3"), ("3", "5"), ("5", "2"), ("5", "4"), ("4", "6")]
DIJ_LABELS = {"1": "[0, –]", "2": "[4, 5]", "3": "[1, 1]",
              "4": "[6, 5]", "5": "[3, 3]", "6": "[8, 4]"}  # fmt: skip

#: All-pairs network: undirected edges (both directions, same length) and
#: directed arcs (two opposite pairs, one arc of each pair negative).
AP_NODES = {"1": (40, 135), "2": (175, 45), "3": (175, 225),
            "4": (420, 45), "5": (420, 225), "6": (555, 135)}  # fmt: skip
AP_EDGES = [
    ("1", "2", 1), ("1", "3", 2), ("2", "3", 4), ("4", "5", 7),
    ("2", "5", 9), ("3", "4", 8), ("4", "6", 5), ("5", "6", 3),
]  # fmt: skip
AP_ARCS = [("4", "2", "−2"), ("2", "4", 3), ("3", "5", "−5"), ("5", "3", 6)]

#: Turn-penalty grid: nine nodes of a diamond, arcs directed to the right.
GRID_NODES = {"1": (40, 190), "2": (130, 110), "3": (130, 270),
              "4": (220, 30), "5": (220, 190), "6": (220, 350),
              "7": (310, 110), "8": (310, 270), "9": (400, 190)}  # fmt: skip
GRID_ARCS = [
    ("1", "2", 3), ("1", "3", 1), ("2", "4", 4), ("2", "5", 2), ("3", "5", 5),
    ("3", "6", 4), ("4", "7", 4), ("5", "7", 2), ("5", "8", 3), ("6", "8", 5),
    ("7", "9", 3), ("8", "9", 4),
]  # fmt: skip
GRID_PATH = [("1", "3"), ("3", "5"), ("5", "7"), ("7", "9")]

#: Seven-city route network (distances).
ROUTE_NODES = {"1": (35, 140), "2": (150, 45), "3": (150, 140), "4": (150, 240),
               "5": (285, 85), "6": (285, 200), "7": (390, 140)}  # fmt: skip
ROUTE_ARCS = [
    ("1", "2", 5), ("1", "3", 14), ("1", "4", 6), ("2", "3", 6), ("2", "5", 12),
    ("2", "6", 12), ("3", "5", 2), ("3", "6", 3), ("4", "3", 7), ("4", "6", 9),
    ("5", "7", 4), ("6", "7", 4),
]  # fmt: skip
ROUTE_PATH = [("1", "2"), ("2", "3"), ("3", "5"), ("5", "7")]


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def acyclic_network(outdir, lang="en"):
    """Six-node acyclic network whose nodes are to be numbered.

    Used in chapter 8.6 for the source-node, topological-order and one-pass
    recursion problems. Nodes carry letters only, so that the reader can
    number them.
    """
    body = network_svg(
        ACYC_NODES,
        ACYC_ARCS,
        label_side={("d", "f"): -1, ("d", "c"): 1, ("a", "d"): 1},
        label_pos={("d", "f"): 0.45, ("d", "c"): 0.45},
    )
    return write_svg("08-06-acyclic-network", body, 385, 280, outdir)


def dijkstra_network(outdir, lang="en"):
    """Network with a cycle (2-3-5-2) and non-negative lengths.

    Used in chapter 8.6 for Dijkstra's method in DP form. Labels
    [V(j), P(j)] and the highlighted shortest-path tree from node 1 are the
    result of the example.
    """
    body = network_svg(
        DIJ_NODES,
        DIJ_ARCS,
        highlight_arcs=DIJ_TREE,
        highlight_nodes=["1"],
        node_notes=DIJ_LABELS,
        note_pos={
            "1": "below",
            "2": "above",
            "3": "below",
            "4": "above",
            "5": "below",
            "6": "below",
        },  # fmt: skip
        # The crossing diagonals 3->4 and 5->2: labels away from the crossing.
        label_pos={("3", "4"): 0.28, ("5", "2"): 0.28},
        label_side={("2", "3"): -1},
    )
    return write_svg("08-06-dijkstra-network", body, 565, 260, outdir)


def all_pairs_network(outdir, lang="en"):
    """Mixed network for the all-pairs problems.

    Undirected edges (drawn without arrow heads) may be used in both
    directions with the same length; the two pairs of opposite arcs 2-4 and
    3-5 are directed and each contains one negative arc. Used in chapter 8.6
    (Floyd's and Dantzig's recursions).
    """
    # Undirected layer: edges with labels; nodes are redrawn by the second
    # layer on top, which is harmless.
    body = network_svg(
        AP_NODES,
        AP_EDGES,
        directed=False,
        # The two crossing diagonals 2-5 and 3-4: labels near their ends.
        label_pos={("2", "5"): 0.72, ("3", "4"): 0.28},
        label_side={("2", "3"): -1, ("4", "5"): 1},
    )
    # Directed layer: the two antiparallel pairs, curved apart.
    body += network_svg(AP_NODES, AP_ARCS, antiparallel_bend=0.12)
    return write_svg("08-06-all-pairs-network", body, 600, 270, outdir)


def turn_grid(outdir, lang="en"):
    """Nine-node diamond grid for the turn-penalty problem (penalty 2).

    Arcs are directed to the right. The optimal path 1-3-5-7-9 (length 11,
    two turns, total 15) is highlighted. Used in chapter 8.6.
    """
    body = network_svg(
        GRID_NODES,
        GRID_ARCS,
        highlight_arcs=GRID_PATH,
        highlight_nodes=["1", "9"],
    )
    return write_svg("08-06-turn-grid", body, 440, 380, outdir)


def staged_route(outdir, lang="en"):
    """Seven-city route network and its five-stage drawing.

    Panel (a): the network with distances and the shortest route
    1-2-3-5-7 (17) highlighted. Panel (b): stage k = number of arcs used
    so far, a node may appear in several stages; the dashed arc 7 -> 7 of
    length 0 lets routes with three arcs end at stage 5. Used in chapter 8.6.
    """
    body = network_svg(
        ROUTE_NODES,
        ROUTE_ARCS,
        highlight_arcs=ROUTE_PATH,
        highlight_nodes=["1", "7"],
        label_pos={("2", "6"): 0.72, ("3", "5"): 0.35, ("1", "3"): 0.5},
        label_side={("2", "3"): -1, ("4", "3"): 1},
        offset=(0, 20),
    )
    body += svg_text(10, 14, "(a)", anchor="start", math=False)

    # Panel (b): columns k = 1..5; keys "node@stage", labels = node.
    ox = 450
    layers = {
        1: ["1"],
        2: ["2", "3", "4"],
        3: ["3", "5", "6"],
        4: ["5", "6", "7"],
        5: ["7"],
    }
    ys = {
        1: [160],
        2: [70, 160, 250],
        3: [70, 160, 250],
        4: [70, 160, 250],
        5: [160],
    }
    nodes, labels = {}, {}
    for k, names in layers.items():
        for name, y in zip(names, ys[k]):
            key = f"{name}@{k}"
            nodes[key] = (ox + 40 + 95 * (k - 1), y)
            labels[key] = name
    arcs = []
    lengths = {(i, j): w for i, j, w in ROUTE_ARCS}
    for k in range(1, 5):
        for a in layers[k]:
            for b in layers[k + 1]:
                if (a, b) in lengths:
                    arcs.append((f"{a}@{k}", f"{b}@{k + 1}"))
    arcs.append(("7@4", "7@5"))  # zero-length stay arc
    # Stage lines with labels k = 1..5.
    for k in range(1, 6):
        x = ox + 40 + 95 * (k - 1)
        body += (
            f'<line x1="{x}" y1="40" x2="{x}" y2="282" stroke="{GRID}" '
            f'stroke-width="1.2" stroke-dasharray="3 4"/>'
        )
        body += svg_text(x, 298, f"k = {k}", size=13, color=MUTED)
    body += network_svg(
        nodes,
        arcs,
        node_labels=labels,
        dashed_arcs=[("7@4", "7@5")],
        highlight_arcs=[
            ("1@1", "2@2"),
            ("2@2", "3@3"),
            ("3@3", "5@4"),
            ("5@4", "7@5"),
        ],
    )
    body += svg_text(ox, 14, "(b)", anchor="start", math=False)
    return write_svg("08-06-staged-route", body, 905, 312, outdir)


FIGURES = {
    "08-06-acyclic-network": acyclic_network,
    "08-06-dijkstra-network": dijkstra_network,
    "08-06-all-pairs-network": all_pairs_network,
    "08-06-turn-grid": turn_grid,
    "08-06-staged-route": staged_route,
}
