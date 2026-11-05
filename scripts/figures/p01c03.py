"""Figures of chapter 1.3, Introduction to Operations Research.

Two hand-written SVG diagrams (no data, no randomness):

* ``01-03-modelling-process``: the five phases of an OR study with feedback
  arrows from validation (section sec-01-03-modelling-process);
* ``01-03-roadmap``: the nine parts of the book and their
  dependencies (section sec-01-03-roadmap).

Both figures contain words. The ``lang`` argument selects the label table
``LABELS`` ("en" or "zh") and falls back to English for missing keys.
"""

from .common import (
    HIGHLIGHT_FILL,
    INK,
    MUTED,
    svg_line,
    svg_text,
    write_svg,
)

# Words used in the two diagrams, per edition ("en", "zh"); another edition
# adds its own table and falls back to English for missing keys.
LABELS = {
    "en": {
        # modeling process: one or two lines per box
        "phases": [
            ("Problem", "definition"),
            ("Model", "formulation"),
            ("Solution",),
            ("Validation",),
            ("Implemen-", "tation"),
        ],
        "revise": "revise the model",
        "redefine": "redefine the problem",
        # roadmap: one box per part; chapter numbers in the subtitles
        "part1": "Part I · Preliminaries",
        "part1_sub": "1.1 Linear algebra · 1.2 Matrix calculus · "
        "1.3 Introduction to operations research",
        "theory": "Parts II–VII · Theory and methods",
        "conv": (
            "II · Convex optimization",
            "2.1–2.3 sets, functions, duality",
        ),
        "lp": ("III · Linear programming", "3.1–3.7 models, simplex, duality"),
        "net": ("IV · Network optimization", "4.1–4.4 trees, paths, flows"),
        "dp": ("V · Dynamic programming", "5.1–5.4 stages and recursions"),
        "nlp": (
            "VI · Nonlinear programming",
            "6.1–6.6 conditions, KKT, algorithms",
        ),
        "mip": ("VII · Linearization", "7.1–7.3 logic, products, pieces"),
        "kkt": ("KKT via", "LP duality"),
        "sp": "shortest paths",
        "part8": "Part VIII · Examples and problems for Parts III–VI (8.1–8.9)",
        "part9": "Part IX · Implementation in Python: solvers, models, algorithms (9.1–9.6)",
    },
    "zh": {
        # modeling process: every phase fits on one line in Chinese
        "phases": [
            ("问题定义",),
            ("模型建立",),
            ("求解",),
            ("验证",),
            ("实施",),
        ],
        "revise": "修正模型",
        "redefine": "重新定义问题",
        # roadmap: part names as in the table of contents (第一部分 ...)
        "part1": "第一部分 · 预备知识",
        "part1_sub": "1.1 线性代数基础 · 1.2 向量与矩阵微积分 · 1.3 运筹学导论",
        "theory": "第二至第七部分 · 理论与方法",
        "conv": ("第二部分 · 凸优化", "2.1–2.3 凸集、凸函数、对偶"),
        "lp": ("第三部分 · 线性规划", "3.1–3.7 模型、单纯形法、对偶"),
        "net": ("第四部分 · 网络优化", "4.1–4.4 生成树、最短路、流"),
        "dp": ("第五部分 · 动态规划", "5.1–5.4 阶段与递推"),
        "nlp": ("第六部分 · 非线性规划", "6.1–6.6 最优性条件、KKT、算法"),
        "mip": ("第七部分 · 线性化", "7.1–7.3 逻辑、乘积、分段线性"),
        "kkt": ("经 LP 对偶", "导出 KKT"),
        "sp": "最短路",
        "part8": "第八部分 · 第三至第六部分的例题与习题（8.1–8.9）",
        "part9": "第九部分 · Python 实现：求解器、模型、算法（9.1–9.6）",
    },
}


def _labels(lang):
    """Label table of edition ``lang`` with English fallbacks."""
    table = dict(LABELS["en"])
    table.update(LABELS.get(lang, {}))
    return table


def _box(
    x, y, w, h, *, fill="none", stroke=INK, width=1.4, dashed=False, rx=6
):
    """A rounded rectangle with top-left corner (x, y)."""
    dash = ' stroke-dasharray="6 4"' if dashed else ""
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{width}"{dash}/>'
    )


def _path_arrow(d, *, color=MUTED, dashed=False):
    """A path with an arrow head at its end (muted marker by default)."""
    marker = "arrow-muted" if color == MUTED else "arrow"
    dash = ' stroke-dasharray="5 4"' if dashed else ""
    return (
        f'<path d="{d}" fill="none" stroke="{color}" stroke-width="1.4"{dash} '
        f'marker-end="url(#{marker})"/>'
    )


def modelling_process(outdir, lang="en"):
    """The five phases of an OR study as boxes joined by arrows.

    Forward arrows join consecutive phases; two curved feedback arrows run
    from validation back to model formulation (above the row) and back to
    problem definition (below the row). Used in sec-01-03-modelling-process.
    """
    lab = _labels(lang)
    W, H = 700, 176
    bw, bh, gap = 120, 52, 20  # box width, height, gap between boxes
    y0 = 62  # top edge of the row of boxes
    xs = [10 + i * (bw + gap) for i in range(5)]  # left edges: 10, 150, ...
    body = []

    # Boxes with their one- or two-line titles; validation is highlighted
    # because both feedback arrows start there.
    for i, (x, words) in enumerate(zip(xs, lab["phases"])):
        fill = HIGHLIGHT_FILL if i == 3 else "none"
        body.append(_box(x, y0, bw, bh, fill=fill))
        cx = x + bw / 2
        if len(words) == 1:
            body.append(svg_text(cx, y0 + bh / 2, words[0], math=False))
        else:
            body.append(svg_text(cx, y0 + bh / 2 - 9, words[0], math=False))
            body.append(svg_text(cx, y0 + bh / 2 + 9, words[1], math=False))

    # Forward arrows between consecutive boxes (tip on the next box's edge).
    ym = y0 + bh / 2
    for i in range(4):
        body.append(svg_line(xs[i] + bw, ym, xs[i + 1], ym, arrow=True))

    # Feedback 1: validation (box 4) -> formulation (box 2), arching above.
    x_from = xs[3] + bw / 2
    x_to = xs[1] + bw / 2
    body.append(
        _path_arrow(
            f"M{x_from},{y0} C{x_from},{y0 - 46} {x_to},{y0 - 46} {x_to},{y0 - 1}"
        )
    )
    body.append(
        svg_text(
            (x_from + x_to) / 2,
            y0 - 44,
            lab["revise"],
            size=12,
            color=MUTED,
            math=False,
        )
    )

    # Feedback 2: validation -> problem definition (box 1), arching below.
    x_to = xs[0] + bw / 2
    yb = y0 + bh
    body.append(
        _path_arrow(
            f"M{x_from},{yb} C{x_from},{yb + 52} {x_to},{yb + 52} {x_to},{yb + 1}"
        )
    )
    body.append(
        svg_text(
            (x_from + x_to) / 2,
            yb + 50,
            lab["redefine"],
            size=12,
            color=MUTED,
            math=False,
        )
    )
    return write_svg("01-03-modelling-process", "\n".join(body), W, H, outdir)


def roadmap(outdir, lang="en"):
    """Map of the book: the nine parts and their dependencies.

    Top band: Part I (chapters 1.1-1.3). Inside the dashed frame, the six
    topic parts II-VII in two rows (II, III, IV above; VI, VII, V below).
    Solid arrows: the later part uses the earlier one (Part I -> II and
    III; II -> III and VI; III -> IV and VII). Dashed lines: two parts
    treat the same problem from two sides (KKT conditions from LP duality;
    shortest paths by labeling and by DP). Bottom bands: Parts VIII and IX.
    Used in sec-01-03-roadmap.
    """
    lab = _labels(lang)
    W, H = 760, 500
    bw, bh = 232, 54  # part box size
    cols = [14, 264, 514]  # left edges of the three columns
    mids = [x + bw / 2 for x in cols]  # centers: 130, 380, 630
    row_b, row_c = 136, 300  # top edges of the two rows of topic parts
    body = []

    def group(x, y, words, fill="none"):
        """A part box with a bold title and a muted subtitle."""
        out = [_box(x, y, bw, bh, fill=fill)]
        out.append(
            svg_text(
                x + bw / 2,
                y + 18,
                words[0],
                weight="bold",
                size=12.5,
                math=False,
            )
        )
        out.append(
            svg_text(
                x + bw / 2, y + 38, words[1], size=12, color=MUTED, math=False
            )
        )
        return "".join(out)

    # --- top band: Part I -------------------------------------------------
    body.append(_box(6, 14, 748, 56, fill=HIGHLIGHT_FILL))
    body.append(
        svg_text(380, 33, lab["part1"], weight="bold", size=13.5, math=False)
    )
    body.append(
        svg_text(380, 53, lab["part1_sub"], size=12, color=MUTED, math=False)
    )

    # --- frame of the topic parts and its six boxes -----------------------
    body.append(
        _box(6, 90, 748, 288, stroke=MUTED, width=1.2, dashed=True, rx=10)
    )
    body.append(
        svg_text(
            740,
            108,
            lab["theory"],
            anchor="end",
            size=13,
            color=MUTED,
            math=False,
        )
    )
    body.append(group(cols[0], row_b, lab["conv"], fill=HIGHLIGHT_FILL))
    body.append(group(cols[1], row_b, lab["lp"]))
    body.append(group(cols[2], row_b, lab["net"]))
    body.append(group(cols[0], row_c, lab["nlp"]))
    body.append(group(cols[1], row_c, lab["mip"]))
    body.append(group(cols[2], row_c, lab["dp"]))

    # --- dependencies (solid arrows, tip on the receiving box) -----------
    # Part I -> convex optimization and Part I -> linear programming
    body.append(svg_line(mids[0], 70, mids[0], row_b, arrow=True))
    body.append(svg_line(mids[1], 70, mids[1], row_b, arrow=True))
    # convex optimization -> LP (geometry of LP, LP duality)
    ym = row_b + bh / 2
    body.append(svg_line(cols[0] + bw, ym, cols[1], ym, arrow=True))
    # LP -> network optimization
    body.append(svg_line(cols[1] + bw, ym, cols[2], ym, arrow=True))
    # convex optimization -> NLP, LP -> linearization
    body.append(svg_line(mids[0], row_b + bh, mids[0], row_c, arrow=True))
    body.append(svg_line(mids[1], row_b + bh, mids[1], row_c, arrow=True))

    # --- same problem from two sides (dashed, no arrow head) -------------
    # LP duality <-> KKT conditions (NLP)
    body.append(
        svg_line(290, row_b + bh, 220, row_c, color=MUTED, dashed=True)
    )
    # two-line label between the dashed line and the LP -> linearization arrow
    body.append(
        svg_text(
            274,
            238,
            lab["kkt"][0],
            anchor="start",
            size=12,
            color=MUTED,
            math=False,
        )
    )
    body.append(
        svg_text(
            274,
            254,
            lab["kkt"][1],
            anchor="start",
            size=12,
            color=MUTED,
            math=False,
        )
    )
    # networks <-> DP: shortest paths
    body.append(
        svg_line(mids[2], row_b + bh, mids[2], row_c, color=MUTED, dashed=True)
    )
    body.append(
        svg_text(
            mids[2] + 8,
            246,
            lab["sp"],
            anchor="start",
            size=12,
            color=MUTED,
            math=False,
        )
    )

    # --- Parts VIII and IX as bands across the whole width ---------------
    body.append(
        _box(6, 394, 748, 40, fill=HIGHLIGHT_FILL, stroke=MUTED, width=1.2)
    )
    body.append(svg_text(380, 414, lab["part8"], size=13, math=False))
    body.append(_box(6, 446, 748, 40, stroke=MUTED, width=1.2))
    body.append(svg_text(380, 466, lab["part9"], size=13, math=False))
    return write_svg("01-03-roadmap", "\n".join(body), W, H, outdir)


FIGURES = {
    "01-03-modelling-process": modelling_process,
    "01-03-roadmap": roadmap,
}
