"""Figures of chapter 5.2, Equipment Replacement.

* ``05-02-replacement-network``: the age-versus-year decision network of the
  five-year replacement example (keep and buy arcs with their costs, the
  terminal salvage values) with the two optimal plans highlighted.

Labels are numbers and the book's symbols (language neutral).
"""

from .common import GRID, MUTED, network_svg, svg_line, svg_text, write_svg

# ---------------------------------------------------------------------------
# Data of the running example (thousands of dollars)
# ---------------------------------------------------------------------------

N, I0, PRICE = 5, 2, 50
#: operating cost c(i), trade-in value t(i), salvage value s(i)
C = {0: 10, 1: 13, 2: 20, 3: 40, 4: 70, 5: 100, 6: 100}
T = {1: 32, 2: 21, 3: 11, 4: 5, 5: 0, 6: 0}
S = {1: 25, 2: 17, 3: 8, 4: 0, 5: 0, 6: 0, 7: 0}

#: the two optimal plans BBKBK and BKBBK as node sequences (year, age)
PLANS = [
    [(1, 2), (2, 1), (3, 1), (4, 2), (5, 1), (6, 2)],
    [(1, 2), (2, 1), (3, 2), (4, 1), (5, 1), (6, 2)],
]

#: axis titles per edition
AXIS_TITLES = {"en": ("year", "age"), "zh": ("年份", "役龄")}

# Layout: year k on the x axis, age i on the y axis (age 1 at the bottom).
X0, DX = 170, 100  # x of year 1 and spacing between years (px)
Y1, DY = 340, 44  # y of age 1 and spacing between ages (px)
XK, XB, XA = 30, 72, 122  # x of the key columns K, B and of the age ticks
TOP = Y1 - 6 * DY - 34  # y of the column headers


def _ages(k):
    """Possible ages at the start of year k (eq-05-02-states)."""
    if k == 1:
        return [I0]
    return sorted(set(range(1, k)) | {I0 + k - 1})


def _pos(k, i):
    return (X0 + (k - 1) * DX, Y1 - (i - 1) * DY)


def replacement_network(outdir, lang="en"):
    """Age-versus-year network of Example 5.2.1 (ex-05-02-replacement).

    Nodes are the states (k, i); a keep arc (k, i) -> (k+1, i+1) has length
    c(i), a buy arc (k, i) -> (k+1, 1) has length p - t(i) + c(0) = 60 - t(i).
    Arc lengths depend on the age only, so they are given once per age in a
    key on the left (columns K and B); only the arcs of the two optimal
    plans carry their lengths. The terminal lengths -s(i) are written right
    of the year-6 nodes. Used in sec-05-02-network (fig-05-02-replacement-network).
    """
    nodes, labels = {}, {}
    for k in range(1, N + 2):
        for i in _ages(k):
            name = f"{k},{i}"
            nodes[name] = _pos(k, i)
            labels[name] = ""  # empty circles; axes give year and age

    # The two optimal plans: highlighted arcs and nodes.
    hl_arcs, hl_nodes = set(), set()
    for plan in PLANS:
        names = [f"{k},{i}" for k, i in plan]
        hl_nodes.update(names)
        hl_arcs.update(zip(names[:-1], names[1:]))

    # Decision arcs; only highlighted arcs get a length label.
    arcs, label_side = [], {}
    for k in range(1, N + 1):
        for i in _ages(k):
            u = f"{k},{i}"
            for v, length in (
                (f"{k + 1},{i + 1}", C[i]),
                (f"{k + 1},1", PRICE - T[i] + C[0]),
            ):
                if (u, v) in hl_arcs:
                    arcs.append((u, v, length))
                    # buy arcs (descending or flat): label below the line
                    label_side[(u, v)] = (
                        1 if v.endswith(f",{i + 1}") and i + 1 > 1 else -1
                    )
                else:
                    arcs.append((u, v))

    # Terminal lengths -s(i) right of the year-6 nodes.
    notes = {
        f"{N + 1},{i}": f"\u2212{S[i]}" if S[i] else "0" for i in _ages(N + 1)
    }
    note_pos = {name: (26, 0) for name in notes}

    body = network_svg(
        nodes,
        arcs,
        node_labels=labels,
        node_radius=6,
        label_size=12,
        highlight_arcs=sorted(hl_arcs),
        highlight_nodes=sorted(hl_nodes),
        label_side=label_side,
        node_notes=notes,
        note_pos=note_pos,
    )

    # Guide lines per year (drawn first, under the network) and year ticks.
    guides = ""
    bottom = Y1 + 28
    for k in range(1, N + 2):
        x = X0 + (k - 1) * DX
        guides += svg_line(
            x, Y1 + 12, x, Y1 - 6 * DY - 10, color=GRID, width=1
        )
        body += svg_text(x, bottom, str(k), size=12, color=MUTED)
    body = guides + body

    # Key on the left: per age i, the keep length c(i) and buy length 60 - t(i);
    # then the age ticks. Header of the terminal column on the right.
    body += svg_text(XK, TOP, "K", size=13, weight="bold", math=False)
    body += svg_text(XB, TOP, "B", size=13, weight="bold", math=False)
    body += svg_line(XK - 18, TOP + 14, XB + 18, TOP + 14, color=GRID, width=1)
    for i in range(1, 8):
        y = Y1 - (i - 1) * DY
        body += svg_text(XA, y, str(i), size=12, color=MUTED)
        if i <= 6:
            body += svg_text(XK, y, str(C[i]), size=12)
            body += svg_text(XB, y, str(PRICE - T[i] + C[0]), size=12)
    body += svg_line(XB + 24, TOP - 10, XB + 24, Y1 + 10, color=GRID, width=1)
    x_end = X0 + N * DX
    body += svg_text(x_end + 30, TOP, "\u2212s(i)", size=13, color=MUTED)

    year, age = AXIS_TITLES.get(lang, AXIS_TITLES["en"])
    body += svg_text(
        X0 + 2.5 * DX, bottom + 22, year, size=13, color=MUTED, math=False
    )
    body += svg_text(XA, TOP, age, size=13, color=MUTED, math=False)
    return write_svg(
        "05-02-replacement-network", body, x_end + 70, Y1 + 64, outdir
    )


FIGURES = {
    "05-02-replacement-network": replacement_network,
}
