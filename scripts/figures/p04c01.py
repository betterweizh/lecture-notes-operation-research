"""Figures of chapter 4.1, Networks and Minimum Spanning Trees.

All figures are hand-written SVG built with ``network_svg``:

* ``04-01-road-network``: the five-node road network, undirected with edge
  lengths and directed (terminology example);
* ``04-01-terminology``: six small panels on the same five nodes showing a
  path, a cycle, a network that is not connected, a tree, a spanning tree
  and a directed path;
* ``04-01-sewage-network``: the seven-node sewage network with link costs;
* ``04-01-prim-steps``: two intermediate connected sets of Prim's algorithm on
  the sewage network, with the edges of the cut drawn dashed;
* ``04-01-sewage-mst``: the minimum spanning tree of the sewage network with
  the two tied alternative edges dashed.

Labels are numbers and node names only (language neutral).
"""

from .common import network_svg, svg_text, write_svg

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

#: Five-node road network: node positions in px (scaled from the original
#: drawing) and edges (i, j, length in miles). The directed version orients
#: every edge from i to j as listed, so the arc between 2 and 4 is (2, 4).
ROAD_NODES = {
    "3": (30, 88),
    "1": (130, 40),
    "2": (240, 40),
    "4": (185, 128),
    "5": (352, 112),
}
ROAD_EDGES = [
    ("3", "1", 1),
    ("3", "4", 2),
    ("1", "2", 0.5),
    ("1", "4", 2.5),
    ("1", "5", 5),
    ("2", "5", 1.2),
    ("2", "4", 1.7),
    ("4", "5", 3),
]

#: Sewage network: node positions in px and links (i, j, cost in $ million).
SEWAGE_NODES = {
    "1": (30, 100),
    "2": (175, 42),
    "3": (210, 132),
    "4": (140, 232),
    "5": (400, 64),
    "6": (380, 232),
    "7": (525, 112),
}
SEWAGE_EDGES = [
    ("1", "2", 1),
    ("1", "3", 5),
    ("1", "4", 7),
    ("1", "5", 9),
    ("2", "3", 6),
    ("2", "4", 4),
    ("2", "5", 3),
    ("3", "4", 5),
    ("3", "6", 10),
    ("3", "7", 8),
    ("4", "5", 8),
    ("4", "6", 3),
    ("5", "7", 3),
    ("6", "7", 4),
]

#: Edges of the minimum spanning tree found by Prim (from node 7) and by
#: Kruskal; total cost 3 + 3 + 1 + 4 + 3 + 5 = 19 ($ million).
SEWAGE_MST = [
    ("5", "7"),
    ("2", "5"),
    ("1", "2"),
    ("2", "4"),
    ("4", "6"),
    ("1", "3"),
]
#: Non-tree edges that tie with a tree edge on their cycle (6-7 with 2-4,
#: 3-4 with 1-3): swapping them gives the other three optimal trees.
SEWAGE_TIES = [("6", "7"), ("3", "4")]

# Label nudges for the sewage drawing: several links cross near their
# midpoints (1-5 crosses 2-3 and 2-4; 4-5 crosses 3-6 and 3-7).
SEWAGE_LABEL_POS = {
    ("1", "5"): 0.62,
    ("2", "3"): 0.32,
    ("2", "4"): 0.62,
    ("4", "5"): 0.22,
    ("3", "6"): 0.6,
    ("3", "7"): 0.72,
    ("3", "4"): 0.3,
    ("1", "4"): 0.45,
}
SEWAGE_LABEL_SIDE = {
    ("1", "3"): -1,
    ("3", "4"): -1,
    ("1", "4"): -1,
    ("4", "5"): -1,
}
SEWAGE_LABEL_SHIFT = {
    ("1", "5"): (0, 2),
}


def _road_panel(highlight=(), directed=False, edges=None, muted=(), dashed=(),
                offset=(0, 0), labels=True, hl_nodes=()):  # fmt: skip
    """One drawing of the five-node network (helper for two figures).

    ``edges`` defaults to all road edges with their lengths; ``labels=False``
    drops the lengths. The other arguments are passed to ``network_svg``.
    """
    edges = ROAD_EDGES if edges is None else edges
    arcs = [(u, v, c) if labels else (u, v) for u, v, c in edges]
    body = network_svg(
        ROAD_NODES,
        arcs,
        directed=directed,
        highlight_arcs=highlight,
        highlight_nodes=hl_nodes,
        muted_arcs=muted,
        dashed_arcs=dashed,
        offset=offset,
        # 1-5 crosses 2-4: put its label right of the crossing, below the line.
        label_pos={("1", "4"): 0.42, ("1", "5"): 0.72, ("2", "4"): 0.55},
        label_side={
            ("1", "4"): -1,
            ("3", "4"): -1,
            ("4", "5"): -1,
            ("1", "5"): -1,
        },
    )
    return body


def road_network(outdir, lang="en"):
    """The five-node road network, undirected (a) and directed (b).

    Used in Example "A small road network" (sec-04-01-graphs). Panel (a) shows
    the edge lengths in miles; panel (b) orients each edge as listed in
    ``ROAD_EDGES`` (the arc between nodes 2 and 4 is (2, 4)).
    """
    # Panel (a): undirected with lengths.
    body = _road_panel(offset=(10, 30))
    body += svg_text(14, 14, "(a)", anchor="start", math=False)
    # Panel (b): same layout shifted right, arrow heads, lengths.
    body += _road_panel(directed=True, offset=(410, 30))
    body += svg_text(414, 14, "(b)", anchor="start", math=False)
    return write_svg("04-01-road-network", body, 790, 185, outdir)


def terminology(outdir, lang="en"):
    """Six panels on the five-node network illustrating the definitions.

    Used in sec-04-01-paths and sec-04-01-trees. In each panel the object is drawn
    in the highlight color (edges and the nodes it contains); the remaining
    edges of the road network are dashed gray for context. (a) path 3-1-5; (b) cycle 1-2-4-1; (c) a
    sub-network that is not connected (node 5 cut off); (d) a tree that is
    not spanning (node 5 missing); (e) a spanning tree; (f) a directed path
    3 -> 1 -> 2 -> 5 in the directed network.
    """
    all_edges = [(u, v) for u, v, _ in ROAD_EDGES]

    def rest(chosen):
        # Edges of the road network not in the chosen object (for context).
        keys = {frozenset(e) for e in chosen}
        return [e for e in all_edges if frozenset(e) not in keys]

    panels = [
        # (panel title, edges of the object, nodes not in it, directed?)
        ("(a)", [("3", "1"), ("1", "5")], ["2", "4"], False),
        ("(b)", [("1", "2"), ("2", "4"), ("1", "4")], ["3", "5"], False),
        (
            "(c)",
            [("3", "1"), ("1", "2"), ("1", "4"), ("2", "4"), ("3", "4")],
            [],
            False,
        ),
        ("(d)", [("3", "1"), ("1", "2"), ("1", "4")], ["5"], False),
        ("(e)", [("3", "1"), ("1", "2"), ("2", "5"), ("3", "4")], [], False),
        ("(f)", [("3", "1"), ("1", "2"), ("2", "5")], ["4"], True),
    ]
    body = ""
    pw, ph = 395, 175  # panel size in px
    for idx, (title, chosen, out_nodes, directed) in enumerate(panels):
        col, row = idx % 2, idx // 2
        ox, oy = 10 + col * pw, 28 + row * ph
        # No lengths in these panels: only the structure matters.
        edges = [(u, v, "") for u, v in chosen] + [
            (u, v, "") for u, v in rest(chosen)
        ]
        body += _road_panel(
            highlight=chosen,
            directed=directed,
            edges=edges,
            dashed=rest(chosen),
            muted=rest(chosen),
            offset=(ox, oy),
            labels=False,
            hl_nodes=[n for n in ROAD_NODES if n not in out_nodes],
        )
        body += svg_text(ox + 4, oy - 14, title, anchor="start", math=False)
    return write_svg("04-01-terminology", body, 790, 28 + 3 * ph - 15, outdir)


def sewage_network(outdir, lang="en"):
    """The sewage network: 7 stations, 14 feasible links, costs in $ million.

    Used in Example "A sewage collection network" (sec-04-01-mst-problem).
    Node 7 is the main station.
    """
    body = network_svg(
        SEWAGE_NODES,
        SEWAGE_EDGES,
        directed=False,
        highlight_nodes=["7"],
        label_pos=SEWAGE_LABEL_POS,
        label_side=SEWAGE_LABEL_SIDE,
        label_shift=SEWAGE_LABEL_SHIFT,
        offset=(5, 5),
    )
    return write_svg("04-01-sewage-network", body, 560, 265, outdir)


def prim_steps(outdir, lang="en"):
    """Two intermediate states of Prim's algorithm on the sewage network.

    Used in Example "Prim's algorithm on the sewage network" (sec-04-01-prim).
    (a) After iteration 2: connected set C_3 = {2, 5, 7} (filled nodes),
    tree edges 5-7 and 2-5 thick, the links of the cut between C_3 and the
    other nodes dashed with their costs; the cheapest of them, 1-2 (cost 1),
    enters next. (b) After iteration 4: C_5 = {1, 2, 4, 5, 7}; the cheapest
    cut link is 4-6 (cost 3).
    """
    states = [
        # (title, connected set, tree edges so far)
        ("(a)", {"2", "5", "7"}, [("5", "7"), ("2", "5")]),
        (
            "(b)",
            {"1", "2", "4", "5", "7"},
            [("5", "7"), ("2", "5"), ("1", "2"), ("2", "4")],
        ),
    ]
    body = ""
    for idx, (title, C, tree) in enumerate(states):
        ox, oy = 5, 28 + idx * 280
        tree_keys = {frozenset(e) for e in tree}
        # Cut links: exactly one end in C.
        cut = [(u, v, c) for u, v, c in SEWAGE_EDGES if (u in C) != (v in C)]
        # Tree edges keep their costs; links not in the cut are gray and
        # unlabeled.
        tree_lab = [
            (u, v, c)
            for u, v, c in SEWAGE_EDGES
            if frozenset((u, v)) in tree_keys
        ]
        other = [(u, v, "") for u, v, c in SEWAGE_EDGES
                 if frozenset((u, v)) not in tree_keys and (u in C) == (v in C)]  # fmt: skip
        arcs = tree_lab + cut + other
        body += network_svg(
            SEWAGE_NODES,
            arcs,
            directed=False,
            highlight_arcs=tree,
            highlight_nodes=sorted(C),
            dashed_arcs=[(u, v) for u, v, _ in cut],
            muted_arcs=[(u, v) for u, v, _ in other],
            label_pos=SEWAGE_LABEL_POS,
            label_side=SEWAGE_LABEL_SIDE,
            label_shift=SEWAGE_LABEL_SHIFT,
            offset=(ox, oy),
        )
        body += svg_text(ox + 4, oy - 14, title, anchor="start", math=False)
    return write_svg("04-01-prim-steps", body, 560, 28 + 2 * 280 - 20, outdir)


def sewage_mst(outdir, lang="en"):
    """Minimum spanning tree of the sewage network (cost 19, $ million).

    Used in sec-04-01-prim and sec-04-01-kruskal. Tree edges thick in the
    highlight color with their costs; the two non-tree links that tie with
    a tree link on their cycle (6-7 ties with 2-4, 3-4 with 1-3) dashed with
    their costs; all other links gray without labels.
    """
    tree_keys = {frozenset(e) for e in SEWAGE_MST}
    tie_keys = {frozenset(e) for e in SEWAGE_TIES}
    arcs = []
    for u, v, c in SEWAGE_EDGES:
        k = frozenset((u, v))
        arcs.append((u, v, c) if (k in tree_keys or k in tie_keys) else (u, v))
    others = [(u, v) for u, v, _ in SEWAGE_EDGES
              if frozenset((u, v)) not in tree_keys | tie_keys]  # fmt: skip
    body = network_svg(
        SEWAGE_NODES,
        arcs,
        directed=False,
        highlight_arcs=SEWAGE_MST,
        highlight_nodes=["7"],
        dashed_arcs=SEWAGE_TIES,
        muted_arcs=others,
        label_pos=SEWAGE_LABEL_POS,
        label_side=SEWAGE_LABEL_SIDE,
        label_shift=SEWAGE_LABEL_SHIFT,
        offset=(5, 5),
    )
    return write_svg("04-01-sewage-mst", body, 560, 265, outdir)


FIGURES = {
    "04-01-road-network": road_network,
    "04-01-terminology": terminology,
    "04-01-sewage-network": sewage_network,
    "04-01-prim-steps": prim_steps,
    "04-01-sewage-mst": sewage_mst,
}
