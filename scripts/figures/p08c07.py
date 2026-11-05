"""Figures of chapter 8.7, DP Applications (Part VIII, Examples and Problems).

* ``08-07-allocation-network``: the layered longest-path network of the
  one-resource allocation with five units (four activities), with the two
  optimal allocations highlighted.

Labels are numbers and symbols (language neutral).
"""

import html

from .common import MUTED, network_svg, svg_math, svg_text, write_svg

#: returns r_k(y), y = 0..5, of the four activities (data of chapter 5.3)
R = [
    [0, 3, 7, 10, 12, 13],
    [0, 1, 2, 4, 8, 13],
    [0, 2, 4, 6, 8, 10],
    [0, 1, 3, 6, 9, 12],
]
X = 5
#: optimal allocations (3, 0, 2, 0) and (4, 0, 1, 0) as node sequences
#: (stage k, units left y); the last activity takes the rest (y = 0 here)
PATHS = [
    [(1, 5), (2, 2), (3, 2), (4, 0)],
    [(1, 5), (2, 1), (3, 1), (4, 0)],
]

# Layout: stage k on the x axis, units left y on the vertical axis (5 at the top).
X0, DX = 80, 150
Y0, DY = 300, 50


def _pos(k, y):
    return (X0 + (k - 1) * DX, Y0 - y * DY)


def allocation_network(outdir, lang="en"):
    """Longest-path network of the allocation of 5 units (ex-08-07-l9-longest-path).

    Node (k, y): activity k is next and y units are left. An arc
    (k, y) -> (k+1, y') allocates x = y - y' units to activity k and has
    length r_k(x). Stage 1 has the single node (1, 5). The numbers right
    of the stage-4 nodes are the terminal lengths r_4(y) (the last activity
    takes the rest). Only the arcs out of stage 1 and the arcs of the two
    optimal paths carry their lengths, as in the text.
    """
    nodes, labels = {"1,5": _pos(1, X)}, {"1,5": ""}
    for k in (2, 3, 4):
        for x in range(X + 1):
            nodes[f"{k},{x}"] = _pos(k, x)
            labels[f"{k},{x}"] = ""

    hl = set()
    for path in PATHS:
        names = [f"{k},{x}" for k, x in path]
        hl.update(zip(names[:-1], names[1:]))

    arcs, label_pos, label_side = [], {}, {}
    # stage 1 -> 2: all six arcs labeled with r_1(x), labels near the head
    for x2 in range(X + 1):
        u, v = "1,5", f"2,{x2}"
        arcs.append((u, v, R[0][X - x2]))
        label_pos[(u, v)] = 0.78
        label_side[(u, v)] = 1
    # stages 2 -> 3 and 3 -> 4: every allocation y <= x; labels on optimal arcs
    for k in (2, 3):
        for x in range(X + 1):
            for x2 in range(x + 1):
                u, v = f"{k},{x}", f"{k + 1},{x2}"
                if (u, v) in hl:
                    arcs.append((u, v, R[k - 1][x - x2]))
                    label_side[(u, v)] = -1
                else:
                    arcs.append((u, v))

    notes = {f"4,{x}": str(R[3][x]) for x in range(X + 1)}
    note_pos = {name: (24, 0) for name in notes}
    body = network_svg(
        nodes,
        arcs,
        node_labels=labels,
        node_radius=7,
        label_size=12,
        highlight_arcs=sorted(hl),
        highlight_nodes=sorted({n for a in hl for n in a}),
        muted_arcs=[
            (a[0], a[1])
            for a in arcs
            if (a[0], a[1]) not in hl and not a[0].startswith("1,")
        ],
        label_pos=label_pos,
        label_side=label_side,
        node_notes=notes,
        note_pos=note_pos,
    )
    # stage ticks below, units-left ticks on the left, header of the terminal column
    for k in range(1, 5):
        body += svg_text(
            X0 + (k - 1) * DX, Y0 + 30, str(k), size=12, color=MUTED
        )
    for x in range(X + 1):
        body += svg_text(X0 - 40, Y0 - x * DY, str(x), size=12, color=MUTED)
    body += svg_text(
        X0 + 3 * DX + 26, Y0 - X * DY - 34, "r_4(y)", size=13, color=MUTED
    )
    # axis titles: a word in upright text followed by the italic symbol
    (word_k, sym_k), (word_x, sym_x) = AXIS_TITLES.get(lang, AXIS_TITLES["en"])
    title_k = f"{html.escape(word_k)} {svg_math(sym_k)}"
    title_x = f"{html.escape(word_x)} {svg_math(sym_x)}"
    body += svg_text(
        X0 + 1.5 * DX, Y0 + 54, title_k, size=13, color=MUTED, raw=True
    )
    body += svg_text(
        X0 - 40, Y0 - X * DY - 34, title_x, size=13, color=MUTED, raw=True
    )
    return write_svg(
        "08-07-allocation-network", body, X0 + 3 * DX + 70, Y0 + 72, outdir
    )


#: axis titles per edition
AXIS_TITLES = {
    "en": (("activity", "k"), ("units left", "y")),
    "zh": (("活动", "k"), ("待分配单位", "y")),
}

FIGURES = {
    "08-07-allocation-network": allocation_network,
}
