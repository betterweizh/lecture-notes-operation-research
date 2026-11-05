"""Figures of chapter 4.3, Maximum Flow.

All figures are hand-written SVG built with ``network_svg``:

* ``04-03-pipeline``: the oil pipeline network (super-source 0, refineries 1-3,
  pumping stations 4-6, terminals 7-8, super-sink 9) with capacities and the
  minimum cut S = {0, 1, 2, 3} as a dashed line;
* ``04-03-residual-example``: the diamond network with a flow x^0 (labels
  x_ij/u_ij) and its residual network G(x^0);
* ``04-03-ff-iterations``: the four residual networks G(x^0)..G(x^3) of the
  Ford-Fulkerson run on the diamond network, augmenting paths highlighted;
* ``04-03-min-cut``: the maximum flow on the diamond network with its two
  minimum cuts.

Labels are numbers, node names and symbols only (language neutral).
"""

from .common import MUTED, network_svg, svg_line, svg_text, write_svg

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

#: Pipeline network: node positions in px (columns: source, refineries,
#: pumping stations, terminals, sink).
PIPE_NODES = {
    "0": (36, 150),
    "1": (130, 58),
    "2": (130, 150),
    "3": (130, 242),
    "4": (262, 92),
    "5": (262, 208),
    "6": (368, 150),
    "7": (478, 78),
    "8": (478, 222),
    "9": (566, 150),
}
#: One-way pipelines (tail, head, capacity in thousand barrels per day);
#: "\\infty" marks the artificial arcs of the super-source and super-sink.
#: The capacity 30 of 5 -> 8 is the value assumed in the chapter.
PIPE_DIRECTED = [
    ("0", "1", "\\infty"),
    ("0", "2", "\\infty"),
    ("0", "3", "\\infty"),
    ("2", "4", 10),
    ("2", "5", 20),
    ("4", "7", 10),
    ("5", "8", 30),
    ("6", "7", 50),
    ("6", "8", 20),
    ("7", "9", "\\infty"),
    ("8", "9", "\\infty"),
]
#: Two-way pipelines (each is two opposite arcs of this capacity).
PIPE_UNDIRECTED = [
    ("1", "4", 20),
    ("3", "5", 15),
    ("4", "5", 20),
    ("4", "6", 10),
    ("5", "6", 30),
]

#: Diamond network: positions and capacities u_ij (source 1, sink 4).
DIAMOND_NODES = {
    "1": (34, 110),
    "2": (134, 34),
    "3": (164, 186),
    "4": (264, 104),
}
DIAMOND_CAP = {
    ("1", "2"): 4,
    ("1", "3"): 2,
    ("3", "2"): 3,
    ("2", "4"): 5,
    ("3", "4"): 1,
}

#: Residual capacities r_ij > 0 of G(x^0), ..., G(x^3) in the Ford-Fulkerson
#: run (recomputed in verify_ch16.py), and the augmenting path found in each.
FF_RESIDUALS = [
    {("1", "2"): 4, ("1", "3"): 2, ("3", "2"): 3, ("2", "4"): 5, ("3", "4"): 1},
    {("1", "3"): 2, ("2", "1"): 4, ("2", "4"): 1, ("3", "2"): 3, ("3", "4"): 1,
     ("4", "2"): 4},
    {("1", "3"): 1, ("2", "1"): 4, ("2", "3"): 1, ("3", "1"): 1, ("3", "2"): 2,
     ("3", "4"): 1, ("4", "2"): 5},
    {("2", "1"): 4, ("2", "3"): 1, ("3", "1"): 2, ("3", "2"): 2, ("4", "2"): 5,
     ("4", "3"): 1},
]  # fmt: skip
FF_PATHS = [
    [("1", "2"), ("2", "4")],
    [("1", "3"), ("3", "2"), ("2", "4")],
    [("1", "3"), ("3", "4")],
    [],
]


def pipeline(outdir, lang="en"):
    """Pipeline network of Example 4.3.1 with capacities and the minimum cut.

    Two-way pipelines are drawn without arrow heads; the dashed vertical
    line separates S = {0, 1, 2, 3} from the rest and crosses exactly the
    forward arcs 1-4, 2->4, 2->5 and 3-5 (capacity 65).
    """
    # Directed arcs first, then the undirected links on top (node circles
    # are drawn twice at the same place, which is harmless).
    body = network_svg(
        PIPE_NODES,
        PIPE_DIRECTED,
        label_pos={("2", "4"): 0.42, ("2", "5"): 0.42, ("6", "7"): 0.45},
        label_side={
            ("2", "5"): -1,
            ("6", "8"): -1,
            ("8", "9"): -1,
            ("0", "3"): -1,
        },
    )
    body += network_svg(
        PIPE_NODES,
        PIPE_UNDIRECTED,
        directed=False,
        label_pos={("1", "4"): 0.72, ("3", "5"): 0.72},
        label_side={("3", "5"): -1, ("4", "6"): 1, ("5", "6"): -1},
    )
    # The minimum cut: a dashed vertical line between refineries and stations.
    body += svg_line(196, 18, 196, 282, color=MUTED, dashed=True)
    body += svg_text(196, 296, "65", size=13, color=MUTED)
    return write_svg("04-03-pipeline", body, 600, 306, outdir)


def _diamond(arcs, *, offset=(0, 0), highlight=(), title=None):
    """One panel with the diamond network (helper for three figures)."""
    body = network_svg(
        DIAMOND_NODES,
        arcs,
        highlight_arcs=highlight,
        offset=offset,
        label_pos={("3", "2"): 0.45, ("2", "3"): 0.4, ("3", "1"): 0.72},
    )
    if title:
        # Panel titles like "(a)" are plain text; "G(x^0)" is math.
        plain = title.startswith("(")
        body += svg_text(
            offset[0] + 4,
            offset[1] + 4,
            title,
            anchor="start",
            size=14,
            math=not plain,
        )
    return body


def residual_example(outdir, lang="en"):
    """Example 4.3.2: flow x^0 (2 units on 1-2-4) and its residual network.

    Panel (a): labels x_ij/u_ij. Panel (b): every r_ij > 0 of G(x^0); the
    backward arcs (2, 1) and (4, 2) carry the 2 units that can be pushed back.
    """
    x0 = {("1", "2"): 2, ("2", "4"): 2}
    arcs_a = [
        (i, j, f"{x0.get((i, j), 0)}/{u}") for (i, j), u in DIAMOND_CAP.items()
    ]
    # r_ij = u_ij - x_ij + x_ji for both directions of every arc.
    res = {}
    for (i, j), u in DIAMOND_CAP.items():
        res[(i, j)] = u - x0.get((i, j), 0)
        res[(j, i)] = x0.get((i, j), 0)
    arcs_b = [(i, j, r) for (i, j), r in res.items() if r > 0]
    body = _diamond(arcs_a, offset=(10, 26), title="(a)")
    body += _diamond(arcs_b, offset=(320, 26), title="(b)")
    return write_svg("04-03-residual-example", body, 620, 236, outdir)


def ff_iterations(outdir, lang="en"):
    """Example 4.3.3: residual networks G(x^0)..G(x^3) in a 2 x 2 grid.

    The augmenting path found by the labeling procedure is highlighted in
    each panel; G(x^3) has no residual arc leaving the source.
    """
    body = ""
    for k, (res, path) in enumerate(zip(FF_RESIDUALS, FF_PATHS)):
        ox, oy = (10 + 320 * (k % 2), 30 + 236 * (k // 2))
        arcs = [(i, j, r) for (i, j), r in res.items()]
        body += _diamond(
            arcs, offset=(ox, oy), highlight=path, title=f"G(\\vx^{{({k})}})"
        )
    return write_svg("04-03-ff-iterations", body, 620, 470, outdir)


def min_cut(outdir, lang="en"):
    """Maximum flow of Example 4.3.3 (labels x_ij/u_ij) and its two minimum
    cuts.

    Left dashed line: S = {1} (arcs 1->2, 1->3, capacity 6). Right dashed
    line: S = {1, 2, 3} (arcs 2->4, 3->4, capacity 6).
    """
    xopt = {
        ("1", "2"): 4,
        ("1", "3"): 2,
        ("3", "2"): 1,
        ("2", "4"): 5,
        ("3", "4"): 1,
    }
    arcs = [(i, j, f"{xopt[(i, j)]}/{u}") for (i, j), u in DIAMOND_CAP.items()]
    body = _diamond(arcs, offset=(20, 26))
    # Cut lines (coordinates include the panel offset of 20, 26).
    body += svg_line(80, 26, 80, 232, color=MUTED, dashed=True)
    body += svg_line(252, 26, 252, 232, color=MUTED, dashed=True)
    body += svg_text(80, 246, "6", size=13, color=MUTED)
    body += svg_text(252, 246, "6", size=13, color=MUTED)
    return write_svg("04-03-min-cut", body, 310, 256, outdir)


FIGURES = {
    "04-03-pipeline": pipeline,
    "04-03-residual-example": residual_example,
    "04-03-ff-iterations": ff_iterations,
    "04-03-min-cut": min_cut,
}
