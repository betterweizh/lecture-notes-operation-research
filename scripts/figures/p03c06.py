"""Figures of chapter 3.6 (Revised and Dual Simplex Methods).

* ``03-06-revised-iteration``: one iteration of the revised simplex method as a
  cycle of four steps (price, choose the entering variable, ratio test,
  update the inverse), with the two exits (optimal, unbounded). Used in the
  section "The revised simplex algorithm".
"""

from .common import BLUE, ORANGE, svg_line, svg_text, write_svg

# Step titles in each edition (the formulas are language-neutral).
WORDS = {
    "en": {
        "price": "1  Price",
        "enter": "2  Choose entering variable",
        "ratio": "3  Ratio test",
        "update": "4  Update the inverse",
        "optimal": "optimal",
        "unbounded": "unbounded",
        # words around the formulas: "all cbar_j >= 0: optimal", "-> row r"
        "all": "all ",
        "colon": ": ",
        "row": "  → row ",
        "row_after": "",
    },
    "zh": {
        "price": "1  定价",
        "enter": "2  选择进基变量",
        "ratio": "3  比值检验",
        "update": "4  更新逆矩阵",
        "optimal": "最优",
        "unbounded": "无界",
        "all": "所有 ",
        "colon": "：",
        "row": "  → 第 ",
        "row_after": " 行",
    },
}


def _f(*parts):
    """Formula markup from (text, style) parts.

    style: "b" bold upright (vectors, matrices), "i" italic (scalars),
    "n" upright, "sub"/"sup" an upright sub- or superscript and "subi" an
    italic subscript of the previous part. Baseline shifts are in px for the
    13 px font; the next part returns to the baseline.
    """
    shift = {"sub": 4.0, "subi": 4.0, "sup": -6.0}
    out, level = [], 0.0
    for text, style in parts:
        target = shift.get(style, 0.0)
        dy = target - level
        attrs = f' dy="{dy:g}"' if dy else ""
        if style in ("sub", "subi", "sup"):
            attrs += ' font-size="72%"'
        if style == "b":
            attrs += ' font-weight="bold"'
        if style in ("i", "subi"):
            attrs += ' font-style="italic"'
        out.append(f"<tspan{attrs}>{text}</tspan>")
        level = target
    if level:
        out.append(f'<tspan dy="{-level:g}"></tspan>')
    return "".join(out)


def _box(x, y, w, h, title, formula, color):
    """A rounded box with a bold title line and a formula line."""
    body = (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" ry="8" '
        f'fill="none" stroke="{color}" stroke-width="1.6"/>'
    )
    body += svg_text(
        x + w / 2,
        y + 20,
        title,
        size=13,
        weight="bold",
        math=False,
        color=color,
    )
    body += svg_text(x + w / 2, y + 46, formula, size=13, raw=True)
    return body


def revised_iteration(outdir, lang="en"):
    """Cycle of one revised simplex iteration (see the module docstring)."""
    t = WORDS.get(lang, WORDS["en"])
    W, H = 600, 300
    bw, bh = 250, 64  # box size
    # Box positions: top-left, top-right, bottom-right, bottom-left.
    x1, x2 = 30, 320
    y1, y2 = 40, 190

    # Step 1: pi^T = c_B^T B^-1, cbar_j = c_j - pi^T A_.j
    price = _f(
        ("π", "b"),
        ("T", "sup"),
        (" = ", "n"),
        ("c", "b"),
        ("B", "subi"),
        ("T", "sup"),
        ("B", "b"),
        ("−1", "sup"),
        (",   ", "n"),
        ("c̄", "i"),
        ("j", "subi"),
        (" = ", "n"),
        ("c", "i"),
        ("j", "subi"),
        (" − ", "n"),
        ("π", "b"),
        ("T", "sup"),
        ("A", "b"),
        ("·", "sub"),
        ("j", "subi"),
    )
    # Step 2: cbar_s < 0, Abar_.s = B^-1 A_.s
    enter = _f(
        ("c̄", "i"),
        ("s", "subi"),
        (" &lt; 0,   ", "n"),
        ("Ā", "b"),
        ("·", "sub"),
        ("s", "subi"),
        (" = ", "n"),
        ("B", "b"),
        ("−1", "sup"),
        ("A", "b"),
        ("·", "sub"),
        ("s", "subi"),
    )
    # Step 3: theta = min bbar_i / abar_is over abar_is > 0 -> row r
    ratio = _f(
        ("θ", "i"),
        (" = min ", "n"),
        ("b̄", "i"),
        ("i", "subi"),
        (" / ", "n"),
        ("ā", "i"),
        ("is", "subi"),
        (t["row"], "n"),
        ("r", "i"),
        # "第 r 行" needs a word after r; English has none (no empty tspan)
        *([(t["row_after"], "n")] if t["row_after"] else []),
    )
    # Step 4: B^-1 <- E B^-1
    update = _f(
        ("B", "b"),
        ("−1", "sup"),
        (" ← ", "n"),
        ("E", "b"),
        ("B", "b"),
        ("−1", "sup"),
    )

    body = ""
    body += _box(x1, y1, bw, bh, t["price"], price, BLUE)
    body += _box(x2, y1, bw, bh, t["enter"], enter, BLUE)
    body += _box(x2, y2, bw, bh, t["ratio"], ratio, BLUE)
    body += _box(x1, y2, bw, bh, t["update"], update, BLUE)

    # Cycle arrows 1 -> 2 -> 3 -> 4 -> 1.
    body += svg_line(x1 + bw, y1 + bh / 2, x2, y1 + bh / 2, arrow=True)
    body += svg_line(x2 + bw / 2, y1 + bh, x2 + bw / 2, y2, arrow=True)
    body += svg_line(x2, y2 + bh / 2, x1 + bw, y2 + bh / 2, arrow=True)
    body += svg_line(x1 + bw / 2, y2, x1 + bw / 2, y1 + bh, arrow=True)

    # Exit from step 1: all reduced costs non-negative -> optimal.
    exit1 = (
        _f(
            (t["all"], "n"),
            ("c̄", "i"),
            ("j", "subi"),
            (" ≥ 0" + t["colon"], "n"),
        )
        + t["optimal"]
    )
    body += svg_line(x1 + 40, y1, x1 + 40, 12, arrow=True)
    body += svg_text(
        x1 + 50, 14, exit1, size=12, anchor="start", raw=True, color=ORANGE
    )
    # Exit from step 3: column without positive entry -> unbounded.
    exit3 = (
        _f(
            ("Ā", "b"),
            ("·", "sub"),
            ("s", "subi"),
            (" ≤ ", "n"),
            ("0", "b"),
            (t["colon"], "n"),
        )
        + t["unbounded"]
    )
    body += svg_line(x2 + bw / 2, y2 + bh, x2 + bw / 2, H - 14, arrow=True)
    body += svg_text(
        x2 + bw / 2 + 10,
        H - 18,
        exit3,
        size=12,
        anchor="start",
        raw=True,
        color=ORANGE,
    )
    return write_svg("03-06-revised-iteration", body, W, H, outdir)


FIGURES = {
    "03-06-revised-iteration": revised_iteration,
}
