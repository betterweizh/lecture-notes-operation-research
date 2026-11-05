"""Figures of chapter 5.1, Dynamic Programming Fundamentals.

All figures are hand-written SVG:

* ``05-01-diamond-network``: the 16-node diamond lattice from A to B with arc
  times, the stage lines k = 0, ..., 6, the level axis y, and the optimal
  path A-C-F-J-M-P-B (length 13) of the backward recursion highlighted,
  with the value V_k(y) written beside every node;
* ``05-01-stage-state``: schematic of one stage of a DP model: a state y at
  stage k, two decisions d leading to states t_k(y, d) at stage k + 1, the
  stage costs c_k(y, d) and the values V_{k+1};
* ``05-01-turn-penalty-path``: the same lattice with the optimal path of the
  turn-penalty problem, A-D-G-K-N-Q-B (time 14 plus one turn at K), and the
  pair [V_k(y, U), V_k(y, D)] beside every node.

Labels are node names, numbers and symbols, except the stage titles of
``05-01-stage-state`` ("stage k", translated with ``_t``).
"""

from .common import (
    GRID,
    MUTED,
    ORANGE,
    network_svg,
    svg_line,
    svg_text,
    write_svg,
)

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

#: Lattice coordinates (k, y) of the 16 nodes: stage k = column 0..6, level y.
LATTICE = {
    "A": (0, 0),
    "C": (1, 1), "D": (1, -1),
    "E": (2, 2), "F": (2, 0), "G": (2, -2),
    "H": (3, 3), "I": (3, 1), "J": (3, -1), "K": (3, -3),
    "L": (4, 2), "M": (4, 0), "N": (4, -2),
    "P": (5, 1), "Q": (5, -1),
    "B": (6, 0),
}  # fmt: skip

#: Arcs (tail, head, travel time). Every arc goes one stage to the right,
#: up (y + 1) or down (y - 1).
ARCS = [
    ("A", "C", 1), ("A", "D", 0),
    ("C", "E", 5), ("C", "F", 4), ("D", "F", 7), ("D", "G", 3),
    ("E", "H", 2), ("E", "I", 1), ("F", "I", 1), ("F", "J", 2),
    ("G", "J", 5), ("G", "K", 4),
    ("H", "L", 3), ("I", "L", 3), ("I", "M", 4), ("J", "M", 2),
    ("J", "N", 2), ("K", "N", 2),
    ("L", "P", 5), ("M", "P", 2), ("M", "Q", 8), ("N", "Q", 4),
    ("P", "B", 2), ("Q", "B", 1),
]  # fmt: skip

#: Backward values V_k(y) (minimum time to B), verified by brute force over
#: all 20 paths (scratch script dp_lattice.py).
VALUES = {
    "A": 13, "C": 12, "D": 14, "E": 9, "F": 8, "G": 11,
    "H": 10, "I": 8, "J": 6, "K": 7, "L": 7, "M": 4, "N": 5,
    "P": 2, "Q": 1, "B": 0,
}  # fmt: skip

OPT_PATH = ["A", "C", "F", "J", "M", "P", "B"]

#: Turn-penalty values [V_k(y, U), V_k(y, D)] (leave up / leave down); "–"
#: marks a direction without an arc. Verified in dp_lattice.py.
TURN_VALUES = {
    "A": (20, 17), "C": (20, 16), "D": (21, 17),
    "E": (15, 15), "F": (14, 12), "G": (14, 14),
    "H": (None, 10), "I": (13, 14), "J": (9, 10), "K": (7, None),
    "L": (None, 7), "M": (7, 12), "N": (5, None),
    "P": (None, 2), "Q": (1, None),
}  # fmt: skip

TURN_PATH = ["A", "D", "G", "K", "N", "Q", "B"]

# Drawing geometry: px per stage, px per level, origin of (k, y) = (0, 0).
DX, DY = 92, 52
X0, Y0 = 100, 205


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


def _pos(k, y):
    """Pixel position of lattice point (k, y) (y up in the lattice)."""
    return (X0 + DX * k, Y0 - DY * y)


def _nodes():
    return {n: _pos(*kv) for n, kv in LATTICE.items()}


def _label_sides():
    """Put every arc label on the outer side of the diamond's cell so labels
    of the two arcs entering a node do not collide: arcs going up get their
    label above-left, arcs going down below-left (default side +1/-1)."""
    side = {}
    for u, v, _ in ARCS:
        up = LATTICE[v][1] > LATTICE[u][1]
        side[(u, v)] = 1 if up else -1
    return side


def _label_positions():
    """Labels of up arcs sit near the head (t = 0.6) and labels of down arcs
    near the tail (t = 0.42), which keeps them clear of the value labels
    written above the nodes."""
    return {
        (u, v): (0.6 if LATTICE[v][1] > LATTICE[u][1] else 0.42)
        for u, v, _ in ARCS
    }


def _stage_lines(top, bottom, labels=True):
    """Dotted vertical line through every stage with the label k below."""
    out = ""
    for k in range(7):
        x = X0 + DX * k
        out += (
            f'<line x1="{x}" y1="{top}" x2="{x}" y2="{bottom}" '
            f'stroke="{GRID}" stroke-width="1.2" stroke-dasharray="3 4"/>'
        )
        if labels:
            out += svg_text(x, bottom + 14, f"k = {k}", size=13, color=MUTED)
    return out


def _path_arcs(path):
    return list(zip(path, path[1:]))


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def diamond_network(outdir, lang="en"):
    """The diamond lattice with arc times, stages, values and the optimum.

    Used in chapter 5.1 for the running example (backward and forward
    recursion) and again in chapter 8.6 for the lattice problems. Values
    V_k(y) are written in small boxes beside the nodes; the optimal path
    A-C-F-J-M-P-B (13) is highlighted.
    """
    nodes = _nodes()
    # Stage lines first so that the network is drawn on top of them.
    body = _stage_lines(28, 382)
    # Level axis on the left: ticks y = -3..3.
    ax = X0 - 62
    body += svg_line(ax, 380, ax, 30, color=MUTED, arrow=True, width=1.1)
    for y in range(-3, 4):
        yy = Y0 - DY * y
        body += svg_line(ax - 4, yy, ax + 4, yy, color=MUTED, width=1.1)
        # true minus sign (U+2212) for the negative levels
        body += svg_text(
            ax - 8, yy, str(y).replace("-", "\u2212"), size=12, color=MUTED,
            anchor="end",
        )  # fmt: skip
    body += svg_text(ax + 8, 30, "y", size=14, color=MUTED, anchor="start")
    body += network_svg(
        nodes,
        ARCS,
        highlight_arcs=_path_arcs(OPT_PATH),
        highlight_nodes=["A", "B"],
        label_side=_label_sides(),
        label_pos=_label_positions(),
    )
    # Values V_k(y) in bold orange above every node (A and B: to the side),
    # so that they cannot be mistaken for arc times.
    for n, (x, y) in nodes.items():
        dx, dy = {"A": (-30, -22), "B": (30, -22)}.get(n, (0, -27))
        body += svg_text(
            x + dx,
            y + dy,
            str(VALUES[n]),
            size=13,
            color=ORANGE,
            weight="bold",
            halo=True,
        )
    return write_svg("05-01-diamond-network", body, 720, 410, outdir)


def stage_state(outdir, lang="en"):
    """One stage of a DP model (schematic).

    State y at stage k; each decision d moves the system to the state
    t_k(y, d) at stage k + 1 and costs c_k(y, d); the recurrence compares
    c_k(y, d) + V_{k+1}(t_k(y, d)) over the decisions. Used in chapter 5.1,
    section on the elements of a DP model.
    """
    xk, xk1 = 120, 400
    # Stage columns: dotted lines with the stage label underneath.
    body = ""
    it = '<tspan font-style="italic">k</tspan>'
    stage_k = _t(lang, f"stage {it}", f"第 {it} 阶段")
    stage_k1 = _t(lang, f"stage {it} + 1", f"第 {it} + 1 阶段")
    for x, lab in ((xk, stage_k), (xk1, stage_k1)):
        body += (
            f'<line x1="{x}" y1="22" x2="{x}" y2="232" stroke="{GRID}" '
            f'stroke-width="1.2" stroke-dasharray="3 4"/>'
        )
        body += svg_text(x, 248, lab, size=13, color=MUTED, raw=True)
    # Other states of each stage drawn as small gray circles.
    for x, ys in ((xk, (45, 205)), (xk1, (30, 125, 222))):
        for y in ys:
            body += (
                f'<circle cx="{x}" cy="{y}" r="9" fill="none" '
                f'stroke="{MUTED}" stroke-width="1.1"/>'
            )
    # The state y and its two successor states.
    y_state, y1, y2 = 125, 70, 180
    body += svg_line(
        xk, y_state, xk1, y1, arrow=True, shorten_start=18, shorten_end=18
    )
    body += svg_line(
        xk,
        y_state,
        xk1,
        y2,
        arrow=True,
        highlight=True,
        shorten_start=18,
        shorten_end=18,
    )
    for x, y, r in ((xk, y_state, 18), (xk1, y1, 18), (xk1, y2, 18)):
        body += (
            f'<circle cx="{x}" cy="{y}" r="{r}" fill="#ffffff" '
            f'stroke="#2b2b2b" stroke-width="1.4"/>'
        )
    body += svg_text(xk, y_state, "y")
    # primes are U+2032 / U+2033 (an ASCII ' would read as a quote mark)
    body += svg_text(xk1, y1, "y′", size=13)
    body += svg_text(xk1, y2, "y″", size=13)
    # Decision and cost labels on the two arcs.
    body += svg_text(255, 78, "d′, c_k(y, d′)", size=13, halo=True)
    body += svg_text(255, 174, "d″, c_k(y, d″)", size=13, halo=True)
    # Values beside the states.
    body += svg_text(xk - 26, y_state, "V_k(y)", anchor="end", size=13)
    body += svg_text(xk1 + 26, y1, "V_{k+1}(y′)", anchor="start", size=13)
    body += svg_text(xk1 + 26, y2, "V_{k+1}(y″)", anchor="start", size=13)
    # Transition equations under the successor values.
    body += svg_text(
        xk1 + 26,
        y1 + 22,
        "y′ = t_k(y, d′)",
        anchor="start",
        size=12,
        color=MUTED,
    )
    body += svg_text(
        xk1 + 26,
        y2 + 22,
        "y″ = t_k(y, d″)",
        anchor="start",
        size=12,
        color=MUTED,
    )
    return write_svg("05-01-stage-state", body, 560, 262, outdir)


def turn_penalty_path(outdir, lang="en"):
    """Optimal path of the turn-penalty problem on the diamond lattice.

    Every turn (change between an up arc and a down arc) costs 3. The
    optimum A-D-G-K-N-Q-B takes time 14 and turns once, at K, so its total
    is 17. Beside every node: [V_k(y, U), V_k(y, D)], the least total from
    that node when it is left upwards / downwards ("–" = no such arc).
    Used in chapter 5.1, state augmentation.
    """
    nodes = _nodes()
    body = _stage_lines(22, 388)

    def fmt(v):
        return "–" if v is None else str(v)

    body += network_svg(
        nodes,
        ARCS,
        highlight_arcs=_path_arcs(TURN_PATH),
        highlight_nodes=["K"],
        label_side=_label_sides(),
        label_pos=_label_positions(),
    )
    # Pairs [V(U), V(D)] in bold orange above the nodes (A to the left,
    # B to the right, the lowest node K below).
    for n, (x, y) in nodes.items():
        if n == "B":
            text, (dx, dy) = "0", (28, -22)
        else:
            a, b = TURN_VALUES[n]
            text = f"[{fmt(a)}, {fmt(b)}]"
            dx, dy = {"A": (-40, -24), "K": (0, 28)}.get(n, (0, -27))
        body += svg_text(
            x + dx,
            y + dy,
            text,
            size=12,
            color=ORANGE,
            weight="bold",
            halo=True,
            math=False,
        )
    # Mark the single turn at K with "+3" in the highlight color.
    kx, ky = _pos(3, -3)
    body += svg_text(
        kx + 34,
        ky + 4,
        "+3",
        size=13,
        color="#2a78d6",
        weight="bold",
        anchor="start",
        math=False,
    )
    return write_svg("05-01-turn-penalty-path", body, 720, 420, outdir)


FIGURES = {
    "05-01-diamond-network": diamond_network,
    "05-01-stage-state": stage_state,
    "05-01-turn-penalty-path": turn_penalty_path,
}
