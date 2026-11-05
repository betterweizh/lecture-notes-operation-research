"""Figures of chapter 8.5, Network Optimization (Part VIII, Examples and Problems).

Every figure redraws a problem network from its arc list (the numbers are
those of the chapter), as hand-written SVG with ``network_svg``. Solutions
are highlighted where the figure sits inside a solution.

Labels are numbers, node names and symbols only (language neutral); the
names of the children in the hide-and-seek problem are proper names.
"""

from .common import network_svg, svg_text, write_svg


def _title(x, y, text):
    """Plain panel title such as "(i)" or "(a)" (not italicized)."""
    return svg_text(x, y, text, anchor="start", size=14, math=False)


# ---------------------------------------------------------------------------
# Terminology
# ---------------------------------------------------------------------------


def directed_networks(outdir, lang="en"):
    """Example 8.5.1: the two directed networks (i) and (ii)."""
    # (i) five nodes in a pentagon; 5 -> 1 closes the cycles through node 1.
    nodes_i = {
        "1": (40, 80),
        "2": (115, 26),
        "5": (190, 80),
        "3": (68, 170),
        "4": (162, 170),
    }
    arcs_i = [("1", "2"), ("2", "5"), ("5", "1"), ("1", "3"), ("3", "4"), ("3", "5"),
              ("4", "5"), ("4", "2")]  # fmt: skip
    # (ii) four nodes in a diamond; acyclic.
    nodes_ii = {
        "1": (30, 100),
        "3": (115, 30),
        "4": (200, 100),
        "2": (115, 170),
    }
    arcs_ii = [("1", "3"), ("1", "2"), ("2", "3"), ("3", "4"), ("2", "4")]
    body = network_svg(nodes_i, arcs_i, offset=(20, 30))
    body += _title(8, 14, "(i)")
    body += network_svg(nodes_ii, arcs_ii, offset=(290, 30))
    body += _title(280, 14, "(ii)")
    return write_svg("08-05-directed-networks", body, 520, 230, outdir)


def draw_network(outdir, lang="en"):
    """Example 8.5.2: the network N = {1..5} with nine arcs (3,4) and (4,3)."""
    nodes = {
        "1": (36, 105),
        "2": (140, 36),
        "3": (264, 36),
        "5": (140, 174),
        "4": (264, 174),
    }
    arcs = [("1", "2"), ("1", "5"), ("2", "3"), ("2", "4"), ("3", "4"), ("3", "5"),
            ("4", "3"), ("4", "5"), ("5", "2")]  # fmt: skip
    return write_svg(
        "08-05-draw-network", network_svg(nodes, arcs), 300, 210, outdir
    )


# ---------------------------------------------------------------------------
# Minimum spanning trees
# ---------------------------------------------------------------------------


def mst_eight(outdir, lang="en"):
    """Example 8.5.4: eight-node network with its minimum spanning tree (14)."""
    nodes = {"s": (36, 160), "v1": (150, 132), "v7": (150, 240), "v2": (260, 40),
             "v3": (330, 130), "v6": (310, 240), "v5": (440, 214), "v4": (470, 90)}  # fmt: skip
    labels = {n: (n if n == "s" else f"v_{n[1]}") for n in nodes}
    edges = [("s", "v1", 3), ("s", "v2", 6), ("s", "v7", 1), ("v1", "v2", 2), ("v1", "v7", 2),
             ("v2", "v3", 1), ("v3", "v4", 6), ("v3", "v6", 4), ("v3", "v5", 3),
             ("v4", "v5", 3), ("v7", "v6", 10), ("v6", "v5", 2)]  # fmt: skip
    tree = [("s", "v7"), ("v7", "v1"), ("v1", "v2"), ("v2", "v3"), ("v3", "v5"),
            ("v5", "v6"), ("v5", "v4")]  # fmt: skip
    body = network_svg(
        nodes, edges, directed=False, node_labels=labels, highlight_arcs=tree,
        bend={("s", "v2"): -0.12}, label_side={("v3", "v6"): -1},
    )  # fmt: skip
    return write_svg("08-05-mst-eight", body, 510, 280, outdir)


def wellheads(outdir, lang="en"):
    """Example 8.5.5: wellhead network, constrained MST and forbidden links.

    Node P is the delivery point (joined to wellhead 1, length not given).
    Highlighted: the high-pressure tree 1-2-3-7-4 (31 miles) and the
    low-pressure tree 1-5-6, 5-9-8 (18 miles). Dashed: cross links.
    """
    # High-pressure wells (2, 3, 4, 7) on the left, low-pressure wells
    # (5, 6, 8, 9) on the right, wellhead 1 and the delivery point P on top.
    nodes = {"P": (300, 30), "1": (300, 110), "2": (150, 100), "3": (60, 220),
             "7": (80, 370), "4": (210, 390), "5": (330, 240), "6": (310, 380),
             "9": (480, 180), "8": (500, 340)}  # fmt: skip
    edges = [("1", "2", 5), ("1", "3", 9), ("1", "4", 20), ("1", "5", 4), ("1", "8", 14),
             ("1", "9", 15), ("2", "3", 6), ("3", "4", 15), ("3", "5", 10), ("4", "5", 20),
             ("4", "6", 7), ("4", "7", 12), ("5", "6", 3), ("5", "7", 5), ("5", "8", 13),
             ("5", "9", 6), ("7", "8", 7), ("8", "9", 5), ("3", "7", 8), ("P", "1", "")]  # fmt: skip
    tree = [("1", "2"), ("2", "3"), ("3", "7"), ("7", "4"),
            ("1", "5"), ("5", "6"), ("5", "9"), ("9", "8")]  # fmt: skip
    cross = [("3", "5"), ("4", "5"), ("4", "6"), ("5", "7"), ("7", "8")]
    body = network_svg(
        nodes, edges, directed=False, highlight_arcs=tree, dashed_arcs=cross,
        muted_arcs=cross,
        # The cross link 7-8 is routed below wells 4 and 6.
        bend={("7", "8"): 0.3},
        label_pos={("1", "4"): 0.35, ("1", "8"): 0.7, ("5", "8"): 0.6, ("4", "5"): 0.4,
                   ("3", "5"): 0.4, ("5", "7"): 0.8, ("1", "3"): 0.4, ("3", "4"): 0.6},
        label_side={("5", "6"): -1, ("5", "9"): -1},
    )  # fmt: skip
    return write_svg("08-05-wellheads", body, 540, 480, outdir)


def machine_cells(outdir, lang="en"):
    """Example 8.5.6: MST of the machine dissimilarities and the three cells.

    The two longest tree links 1-6 and 2-10 are dashed; deleting them gives
    the cells {1}, {2, 4, 6, 7, 9}, {3, 5, 8, 10}. Total length 5.1694.
    """
    nodes = {"1": (36, 60), "6": (150, 60), "9": (150, 170), "7": (264, 60),
             "2": (378, 60), "4": (378, 170), "10": (492, 60), "8": (492, 170),
             "3": (606, 60), "5": (606, 170)}  # fmt: skip
    edges = [("1", "6", "0.8333"), ("6", "9", "0.6000"), ("6", "7", "0.6250"),
             ("7", "2", "0.6000"), ("2", "4", "0.4444"), ("2", "10", "0.7000"),
             ("10", "8", "0.6667"), ("10", "3", "0.5000"), ("3", "5", "0.2000")]  # fmt: skip
    cut = [("1", "6"), ("2", "10")]
    kept = [(u, v) for u, v, _ in edges if (u, v) not in cut]
    body = network_svg(
        nodes,
        edges,
        directed=False,
        highlight_arcs=kept,
        dashed_arcs=cut,
        muted_arcs=cut,
        offset=(0, 10),
    )
    return write_svg("08-05-machine-cells", body, 690, 210, outdir)


# ---------------------------------------------------------------------------
# Shortest paths
# ---------------------------------------------------------------------------


def undirected_network(outdir, lang="en"):
    """Example 8.5.7: undirected network, final labels [V, P], path tree."""
    nodes = {"0": (40, 140), "2": (160, 40), "5": (180, 150), "4": (170, 270),
             "1": (330, 60), "3": (350, 230)}  # fmt: skip
    edges = [("0", "2", 13), ("0", "4", 16), ("0", "5", 8), ("1", "3", 6), ("1", "5", 10),
             ("2", "3", 14), ("2", "5", 11), ("3", "4", 5), ("3", "5", 17), ("4", "5", 7)]  # fmt: skip
    tree = [("0", "5"), ("0", "2"), ("5", "4"), ("5", "1"), ("4", "3")]
    notes = {"0": "[0, –]", "5": "[8, 0]", "2": "[13, 0]", "4": "[15, 5]",
             "1": "[18, 5]", "3": "[20, 4]"}  # fmt: skip
    body = network_svg(
        nodes, edges, directed=False, highlight_arcs=tree, node_notes=notes,
        note_pos={"0": "left", "2": "above", "5": (-44, -26), "4": "below", "1": "right",
                  "3": "right"},
        label_pos={("2", "3"): 0.25, ("3", "5"): 0.5, ("1", "5"): 0.78},
        label_side={("2", "5"): -1, ("4", "5"): -1, ("3", "5"): -1, ("0", "5"): -1},
        offset=(30, 14),
    )  # fmt: skip
    return write_svg("08-05-undirected-network", body, 470, 320, outdir)


def dijkstra_seven(outdir, lang="en"):
    """Example 8.5.8: directed network, shortest-path tree from node 1.

    Arc 2 -> 5 is bent around node 4, which lies on the straight line.
    Node notes are the final distances V(j).
    """
    nodes = {"1": (40, 165), "2": (172, 40), "3": (165, 293), "4": (300, 172),
             "6": (434, 46), "5": (425, 306), "7": (556, 176)}  # fmt: skip
    arcs = [("1", "2", 5), ("1", "3", 1), ("3", "2", 2), ("2", "6", 6), ("6", "2", 7),
            ("2", "4", 7), ("2", "5", 1), ("3", "4", 6), ("4", "3", 7), ("3", "5", 7),
            ("4", "6", 4), ("4", "7", 6), ("5", "4", 3), ("5", "6", 5), ("5", "7", 9),
            ("6", "7", 2)]  # fmt: skip
    tree = [
        ("1", "3"),
        ("3", "2"),
        ("3", "4"),
        ("2", "5"),
        ("2", "6"),
        ("6", "7"),
    ]
    notes = {
        "1": "0",
        "2": "3",
        "3": "1",
        "4": "7",
        "5": "4",
        "6": "9",
        "7": "11",
    }
    body = network_svg(
        nodes, arcs, highlight_arcs=tree, node_notes=notes,
        note_pos={"1": "left", "2": "above", "3": "below", "4": (0, -30), "5": "below",
                  "6": "above", "7": "right"},
        bend={("2", "5"): -0.25},
        label_pos={("2", "5"): 0.4, ("5", "6"): 0.4, ("2", "4"): 0.55, ("4", "7"): 0.78},
        label_side={("1", "3"): -1, ("3", "5"): -1},
        offset=(20, 20),
    )  # fmt: skip
    return write_svg("08-05-dijkstra-seven", body, 620, 370, outdir)


def purchasing(outdir, lang="en"):
    """Example 8.5.9: purchasing network (arc costs in $), shortest path.

    Arc (i, j): order at the start of month i for months i..j-1. Arcs that
    skip one node bend above the row, arcs that skip two bend below, and the
    arc 1 -> 5 bends far above.
    """
    nodes = {str(k): (40 + 115 * (k - 1), 150) for k in range(1, 6)}
    arcs = [("1", "2", "1 700"), ("2", "3", "1 880"), ("3", "4", "2 300"), ("4", "5", "2 720"),
            ("1", "3", "3 968"), ("2", "4", "4 652"), ("3", "5", "4 316"),
            ("1", "4", "7 622"), ("2", "5", "7 244"), ("1", "5", "10 970")]  # fmt: skip
    bend = {("1", "3"): -0.32, ("2", "4"): -0.32, ("3", "5"): -0.32,
            ("1", "4"): 0.28, ("2", "5"): 0.28, ("1", "5"): -0.5}  # fmt: skip
    body = network_svg(
        nodes, arcs, highlight_arcs=[("1", "2"), ("2", "3"), ("3", "5")], bend=bend,
        offset=(0, 20),
        label_pos={("2", "4"): 0.5, ("1", "3"): 0.45, ("3", "5"): 0.55},
    )  # fmt: skip
    return write_svg("08-05-purchasing", body, 540, 270, outdir)


def knapsack(outdir, lang="en"):
    """Example 8.5.10: knapsack states [i, w] and arc values; longest path 80."""
    nodes = {"15": (40, 150), "25": (170, 85), "23": (170, 220), "35": (310, 30),
             "32": (310, 112), "33": (310, 192), "30": (310, 270), "E": (440, 150)}  # fmt: skip
    labels = {k: (k if k == "E" else f"{k[0]},{k[1]}") for k in nodes}
    arcs = [("15", "25", 0), ("15", "23", 30), ("25", "35", 0), ("25", "32", 50),
            ("23", "33", 0), ("23", "30", 50), ("35", "E", 70), ("32", "E", 0),
            ("33", "E", 0), ("30", "E", 0)]  # fmt: skip
    body = network_svg(
        nodes, arcs, node_labels=labels, node_radius=21, node_size=13,
        highlight_arcs=[("15", "23"), ("23", "30"), ("30", "E")],
        label_side={("15", "23"): -1, ("23", "30"): -1, ("30", "E"): -1},
    )  # fmt: skip
    return write_svg("08-05-knapsack", body, 480, 300, outdir)


def hide_and_seek(outdir, lang="en"):
    """Example 8.5.11: "who knows where" network; longest chain Bob -> Joe."""
    nodes = {"Joe": (50, 150), "Bob": (200, 40), "Kay": (330, 150), "Jim": (200, 260),
             "Rae": (470, 260), "Kim": (470, 40)}  # fmt: skip
    arcs = [("Joe", "Bob"), ("Joe", "Jim"), ("Joe", "Rae"), ("Kay", "Bob"), ("Kay", "Jim"),
            ("Kay", "Rae"), ("Jim", "Kay"), ("Bob", "Kay"), ("Rae", "Kim"), ("Kim", "Joe"),
            ("Kim", "Bob")]  # fmt: skip
    chain = [("Bob", "Kay"), ("Kay", "Rae"), ("Rae", "Kim"), ("Kim", "Joe")]
    body = network_svg(
        nodes,
        arcs,
        node_radius=22,
        node_size=13,
        highlight_arcs=chain,
        offset=(0, 10),
    )
    return write_svg("08-05-hide-and-seek", body, 520, 310, outdir)


# ---------------------------------------------------------------------------
# Maximum flow
# ---------------------------------------------------------------------------


def recover_flow(outdir, lang="en"):
    """Example 8.5.12: (a) capacities, (b) residual network after the first
    augmentation with the second augmenting path, (c) final residual network.

    The final residual arc between 3 and 4 points 4 -> 3 (x_34 = u_34 = 5).
    """
    nodes = {"1": (30, 110), "2": (130, 30), "3": (130, 190), "4": (230, 110)}
    caps = [
        ("1", "2", 5),
        ("1", "3", 5),
        ("2", "3", 5),
        ("2", "4", 5),
        ("3", "4", 5),
    ]
    res1 = [
        ("2", "1", 5),
        ("1", "3", 5),
        ("3", "2", 5),
        ("2", "4", 5),
        ("4", "3", 5),
    ]
    res2 = [
        ("2", "1", 5),
        ("3", "1", 5),
        ("2", "3", 5),
        ("4", "2", 5),
        ("4", "3", 5),
    ]
    body = network_svg(nodes, caps, offset=(16, 34)) + _title(6, 14, "(a)")
    body += network_svg(
        nodes,
        res1,
        offset=(296, 34),
        highlight_arcs=[("1", "3"), ("3", "2"), ("2", "4")],
    )
    body += _title(286, 14, "(b)")
    body += network_svg(nodes, res2, offset=(576, 34)) + _title(566, 14, "(c)")
    return write_svg("08-05-recover-flow", body, 840, 250, outdir)


def maxflow_eight(outdir, lang="en"):
    """Example 8.5.13: maximum flow (x/u labels); min cut S = {s, a}."""
    nodes = {"b": (30, 180), "s": (150, 180), "d": (260, 30), "a": (260, 110),
             "e": (380, 180), "t": (510, 180), "c": (270, 255), "f": (270, 340)}  # fmt: skip
    flow = [("s", "a", 6, 7), ("s", "b", 2, 2), ("c", "s", 0, 3), ("a", "d", 3, 3),
            ("a", "e", 3, 3), ("b", "d", 0, 1), ("b", "f", 2, 3), ("c", "f", 1, 4),
            ("d", "t", 3, 5), ("e", "c", 1, 3), ("e", "t", 2, 3), ("f", "t", 3, 4)]  # fmt: skip
    arcs = [(i, j, f"{x}/{u}") for i, j, x, u in flow]
    body = network_svg(
        nodes, arcs, highlight_nodes=["s", "a"],
        label_pos={("b", "d"): 0.35, ("b", "f"): 0.4, ("c", "s"): 0.45},
        label_side={("s", "b"): 1, ("c", "s"): -1, ("e", "c"): 1, ("a", "e"): -1},
        offset=(10, 10),
    )  # fmt: skip
    return write_svg("08-05-maxflow-eight", body, 550, 380, outdir)


def chores(outdir, lang="en"):
    """Example 8.5.14: matching network; assignment Rif 3, Mai 1, Ben 2, Kim 5."""
    kids = ["Rif", "Mai", "Ben", "Kim", "Ken"]
    nodes = {"s": (30, 150), "t": (450, 150)}
    for k, name in enumerate(kids):
        nodes[name] = (160, 30 + 60 * k)
        nodes[str(k + 1)] = (320, 30 + 60 * k)
    pref = {
        "Rif": [3, 4, 5],
        "Mai": [1],
        "Ben": [1, 2],
        "Kim": [1, 2, 5],
        "Ken": [2],
    }
    arcs = [("s", c) for c in kids]
    arcs += [(c, str(j)) for c in kids for j in pref[c]]
    arcs += [(str(j), "t") for j in range(1, 6)]
    match = {"Rif": "3", "Mai": "1", "Ben": "2", "Kim": "5"}
    hl = (
        [("s", c) for c in match]
        + list(match.items())
        + [(j, "t") for j in match.values()]
    )
    body = network_svg(
        nodes,
        arcs,
        node_radius=20,
        node_size=12,
        highlight_arcs=hl,
        offset=(0, 10),
    )
    # height: bottom row at y = 270 + 10 (offset) plus the node radius 20
    return write_svg("08-05-chores", body, 490, 312, outdir)


# ---------------------------------------------------------------------------
# Minimum-cost flow
# ---------------------------------------------------------------------------


def transshipment(outdir, lang="en"):
    """Example 8.5.15: transshipment network with warehouse W and dummy D.

    Arc labels: unit cost, and "(0, 40)" bounds on F1 -> W; node notes are
    the net supplies [b_i]; highlighted arcs carry flow in the optimum.
    """
    nodes = {"F1": (50, 60), "F2": (50, 230), "W": (210, 145), "S1": (380, 40),
             "S2": (380, 145), "S3": (380, 250), "D": (210, 310)}  # fmt: skip
    # Node names are plain upright text (F1, S1, W, D), as in the chapter.
    arcs = [("F1", "W", "2 (0, 40)"), ("F2", "W", 3), ("F1", "S1", 7), ("W", "S1", 3),
            ("W", "S2", 4), ("W", "S3", 6), ("F2", "S3", 5), ("F1", "D", 0), ("F2", "D", 0)]  # fmt: skip
    used = [("F1", "W"), ("F2", "W"), ("F1", "S1"), ("W", "S1"), ("W", "S2"), ("F2", "S3"),
            ("F1", "D")]  # fmt: skip
    notes = {"F1": "[60]", "F2": "[50]", "W": "[0]", "S1": "[−30]", "S2": "[−40]",
             "S3": "[−25]", "D": "[−15]"}  # fmt: skip
    body = network_svg(
        nodes, arcs, highlight_arcs=used, node_notes=notes, math=False,
        note_pos={"F1": "left", "F2": "left", "W": "above", "S1": "right", "S2": "right",
                  "S3": "right", "D": "below"},
        label_pos={("F1", "D"): 0.3, ("F2", "S3"): 0.6, ("F1", "W"): 0.45, ("F2", "W"): 0.72},
        label_side={("F2", "W"): -1, ("F2", "D"): -1},
        offset=(40, 10),
    )  # fmt: skip
    return write_svg("08-05-transshipment", body, 500, 360, outdir)


FIGURES = {
    "08-05-directed-networks": directed_networks,
    "08-05-draw-network": draw_network,
    "08-05-mst-eight": mst_eight,
    "08-05-wellheads": wellheads,
    "08-05-machine-cells": machine_cells,
    "08-05-undirected-network": undirected_network,
    "08-05-dijkstra-seven": dijkstra_seven,
    "08-05-purchasing": purchasing,
    "08-05-knapsack": knapsack,
    "08-05-hide-and-seek": hide_and_seek,
    "08-05-recover-flow": recover_flow,
    "08-05-maxflow-eight": maxflow_eight,
    "08-05-chores": chores,
    "08-05-transshipment": transshipment,
}
