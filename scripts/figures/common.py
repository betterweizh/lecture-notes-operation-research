"""Shared style and drawing helpers for the book's generated figures.

Every figure module ``scripts/figures/pPPcCC.py`` imports what it needs from
here, e.g.::

    import matplotlib.pyplot as plt
    from .common import PALETTE, lp_region, draw_line, save_mpl
    from .common import network_svg, write_svg, svg_text

Contents
--------
Style
    ``INK``, ``MUTED``, ``GRID``, ``PALETTE`` (``BLUE``, ``ORANGE``,
    ``GREEN``, ...), ``HIGHLIGHT``, ``FONT_FAMILY``. The matplotlib style
    (``MPL_STYLE``) is applied when this module is imported.
Output
    ``save_mpl(fig, name, outdir)`` saves a matplotlib figure as SVG.
    ``write_svg(name, body, w, h, outdir)`` wraps hand-written SVG elements
    in a complete SVG file with arrow markers.
Hand-written SVG
    ``svg_math(tex)`` turns a small TeX-like label (``x_{12}``, ``\\pi_1``)
    into SVG text with italic letters and subscripts.
    ``svg_text(...)``, ``svg_line(...)`` and ``network_svg(...)`` build SVG
    elements as strings.
Matplotlib
    ``lp_region(ax, constraints, xlim, ylim)`` shades the feasible region
    of a two-variable LP and returns its vertices; ``draw_line(...)`` draws
    one line a1*x1 + a2*x2 = b (constraints, objective contours).

Conventions
-----------
* A figure function has the signature ``func(outdir, lang="en")``, writes
  exactly one file ``<name>.svg`` into ``outdir`` and returns its path.
* Names are ``PP-CC-descriptive-name`` (part ``PP``, chapter ``CC`` of the
  chapter that first shows the figure), registered in the module's
  ``FIGURES`` dict.
* Hand-written SVG uses pixel coordinates with the origin at the top left
  and y pointing *down*.
"""

from __future__ import annotations

import html
import math as _math
import os
import time
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")  # no display needed; must come before pyplot

import matplotlib.pyplot as plt  # noqa: E402
from cycler import cycler  # noqa: E402

# ---------------------------------------------------------------------------
# Colours and fonts
# ---------------------------------------------------------------------------

INK = "#2b2b2b"  # text, axes, lines, node outlines
MUTED = "#6b6b6b"  # secondary text, de-emphasised lines
GRID = "#d9d9d9"  # light grid lines
WHITE = "#ffffff"  # node fill, label halo (becomes dark in dark mode)

#: Series colours, in order of use. The first three are the house palette;
#: use the others only when a plot needs more than three colours.
PALETTE = [
    "#2a78d6",  # blue
    "#eb6834",  # orange
    "#1baf7a",  # green
    "#8e5bd0",  # purple
    "#d4a20b",  # amber
    "#d64f8f",  # pink
    "#3a9fb0",  # teal
    "#7a7a7a",  # grey
]
BLUE, ORANGE, GREEN, PURPLE, AMBER, PINK, TEAL, GREY = PALETTE

#: Colour of highlighted arcs and nodes (optimal path, chosen tree, ...).
HIGHLIGHT = BLUE
#: Light tint used to fill highlighted nodes.
HIGHLIGHT_FILL = "#e2edfa"

#: Font stack of every figure (matplotlib and hand-written SVG). Latin text
#: uses DejaVu Sans; Chinese labels of the zh edition fall back to a CJK font
#: available on the reader's system (Noto Sans CJK / Source Han Sans,
#: Microsoft YaHei on Windows, PingFang SC on macOS).
FONT_FAMILY = (
    "DejaVu Sans, Noto Sans CJK SC, Source Han Sans SC, Microsoft YaHei, "
    "PingFang SC, Helvetica, Arial, sans-serif"
)
#: Default font size of hand-written SVG text, in px.
FONT_SIZE = 14

#: Matplotlib settings used while the figures of the Chinese edition are
#: drawn (``make_figures.py`` applies them for ``lang="zh"``): a CJK font
#: comes first, so that Chinese labels, also those mixed with ``$…$`` math,
#: are measured and written with real glyphs. The first installed font of
#: the list is used; GitHub Actions installs Noto Sans CJK.
ZH_RC = {
    "font.sans-serif": [
        "Noto Sans CJK SC",
        "Source Han Sans SC",
        "Microsoft YaHei",
        "PingFang SC",
        "SimHei",
        "Noto Sans CJK TC",
        "DejaVu Sans",
    ],
}

# ---------------------------------------------------------------------------
# Matplotlib style (applied on import)
# ---------------------------------------------------------------------------

MPL_STYLE = {
    # Keep text as <text> elements (searchable, crisp, small files).
    "svg.fonttype": "none",
    # Deterministic element ids, so regenerated SVGs do not change.
    "svg.hashsalt": "lecture-notes",
    "figure.figsize": (6.0, 4.0),
    "figure.dpi": 100,
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.transparent": True,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "font.family": "sans-serif",
    "font.sans-serif": [
        "DejaVu Sans",
        # CJK fallbacks for Chinese labels (matplotlib >= 3.6 falls back per
        # glyph); the first installed one is used to measure text.
        "Noto Sans CJK SC",
        "Noto Sans CJK TC",
        "Source Han Sans SC",
        "Microsoft YaHei",
        "PingFang SC",
        "SimHei",
        "Helvetica",
        "Arial",
        "Liberation Sans",
    ],
    "font.size": 11,
    "mathtext.fontset": "dejavusans",
    "text.color": INK,
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "axes.titlecolor": INK,
    "axes.titlesize": 12,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "axes.prop_cycle": cycler(color=PALETTE),
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "xtick.color": INK,
    "ytick.color": INK,
    "legend.frameon": False,
    "lines.linewidth": 1.8,
}
plt.rcParams.update(MPL_STYLE)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def _target(
    name: str, outdir: str | os.PathLike, default_suffix: str = ".svg"
) -> Path:
    """Return ``outdir/name`` with ``default_suffix`` unless ``name`` has one."""
    path = Path(outdir) / name
    if path.suffix.lower() not in (".svg", ".png", ".jpg", ".jpeg", ".pdf"):
        path = Path(outdir) / f"{name}{default_suffix}"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _replace_atomically(tmp: Path, path: Path) -> None:
    """Move the finished ``tmp`` over ``path`` (readers never see half a file).

    On Windows the move fails while another program (an image viewer, a
    browser, OneDrive syncing) holds ``path`` open; retry for a moment, then
    fail with a message that names the file.
    """
    for attempt in range(10):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            if attempt == 9:
                tmp.unlink(missing_ok=True)
                raise PermissionError(
                    f"cannot overwrite {path}: the file is open in another "
                    "program (close it and draw the figure again)"
                ) from None
            time.sleep(0.2)


def save_mpl(fig, name: str, outdir: str | os.PathLike, **savefig_kw) -> Path:
    """Save a matplotlib figure as ``outdir/<name>.svg`` and close it.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    name : str
        Figure name, e.g. ``"03-02-product-mix-region"`` (the ``FIGURES``
        key).
        An explicit ``.png`` suffix saves a raster image at 200 dpi instead.
    outdir : path-like
        Output folder (the ``outdir`` argument of the figure function).
    **savefig_kw
        Extra arguments for ``fig.savefig`` (rarely needed).

    Returns
    -------
    pathlib.Path
        Path of the written file.

    Example
    -------
    ::

        fig, ax = plt.subplots(figsize=(5.5, 4))
        ax.plot(x, y)
        return save_mpl(fig, "06-04-golden-section-geometry", outdir)
    """
    path = _target(name, outdir)
    fmt = path.suffix.lstrip(".").lower()
    options = {"format": fmt, "transparent": True, "bbox_inches": "tight"}
    if fmt == "svg":
        options["metadata"] = {"Date": None}  # reproducible output
    else:
        options["dpi"] = 200
    options.update(savefig_kw)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        # A binary file: matplotlib then writes "\n" line ends on every
        # platform (a text-mode path would get "\r\n" on Windows).
        with open(tmp, "wb") as handle:
            fig.savefig(handle, **options)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    finally:
        plt.close(fig)
    _replace_atomically(tmp, path)
    return path


#: Arrow heads are drawn by markers: ``marker-end="url(#arrow)"`` (ink),
#: ``url(#arrow-hl)`` (highlight colour) and ``url(#arrow-muted)``. The tip
#: sticks out ``ARROW_TIP`` px beyond the end of the line, so a line that
#: should touch a point P must stop ``ARROW_TIP`` px before P (``svg_line``
#: and ``network_svg`` do this for you).
ARROW_TIP = 6.3
_ARROW_LEN = 9.0  # marker size in px (path is drawn in a 10 x 10 box)


def _arrow_marker(marker_id: str, color: str) -> str:
    # Notched arrow head; refX=3 is the notch, so the line ends inside it.
    return (
        f'<marker id="{marker_id}" viewBox="0 0 10 10" refX="3" refY="5" '
        f'markerWidth="{_ARROW_LEN}" markerHeight="{_ARROW_LEN}" '
        f'markerUnits="userSpaceOnUse" orient="auto-start-reverse">'
        f'<path d="M0,0 L10,5 L0,10 L3,5 z" fill="{color}"/></marker>'
    )


ARROW_DEFS = (
    _arrow_marker("arrow", INK)
    + _arrow_marker("arrow-hl", HIGHLIGHT)
    + _arrow_marker("arrow-muted", MUTED)
)


def write_svg(
    name: str,
    body: str,
    w: float,
    h: float,
    outdir: str | os.PathLike,
    *,
    extra_defs: str = "",
) -> Path:
    """Write hand-written SVG elements as a complete file ``outdir/<name>.svg``.

    The root element sets the house font, the default text size
    (``FONT_SIZE``) and ``fill=INK``, and ``<defs>`` holds the arrow
    markers ``arrow``, ``arrow-hl`` and ``arrow-muted`` (see ``ARROW_TIP``).
    The background stays transparent.

    Parameters
    ----------
    name : str
        Figure name, e.g. ``"04-02-dijkstra-network"``.
    body : str
        SVG elements (strings from ``network_svg``, ``svg_text``,
        ``svg_line`` or written by hand), concatenated.
    w, h : float
        Width and height in px; the viewBox is ``0 0 w h``. Choose them so
        that the drawing has a margin of about 10 px on every side.
    outdir : path-like
        Output folder (the ``outdir`` argument of the figure function).
    extra_defs : str, optional
        More ``<defs>`` content (gradients, patterns, other markers).

    Returns
    -------
    pathlib.Path
        Path of the written file.

    Example
    -------
    ::

        body = network_svg(nodes, arcs) + svg_text(20, 20, "(a)")
        return write_svg("04-01-sewage-network", body, 420, 260, outdir)
    """
    path = _target(name, outdir)
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{_num(w)}" height="{_num(h)}" '
        f'viewBox="0 0 {_num(w)} {_num(h)}" font-family="{FONT_FAMILY}" '
        f'font-size="{FONT_SIZE}" fill="{INK}">\n'
        f"<defs>{ARROW_DEFS}{extra_defs}</defs>\n"
        f"{body}\n</svg>\n"
    )
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    # Bytes, not text mode: the same "\n" line ends on every platform.
    tmp.write_bytes(svg.encode("utf-8"))
    _replace_atomically(tmp, path)
    return path


# ---------------------------------------------------------------------------
# Hand-written SVG elements
# ---------------------------------------------------------------------------


def _num(value: float) -> str:
    """Format a coordinate compactly (``12``, ``12.5``, ``12.35``)."""
    text = f"{value:.2f}".rstrip("0").rstrip(".")
    return "0" if text in ("-0", "") else text


_GREEK = {
    "alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "epsilon": "ε",
    "theta": "θ", "lambda": "λ", "mu": "μ", "pi": "π", "rho": "ρ",
    "sigma": "σ", "tau": "τ", "phi": "φ", "psi": "ψ", "omega": "ω",
    "Delta": "Δ", "Phi": "Φ", "Omega": "Ω", "infty": "∞", "le": "≤",
    "ge": "≥", "ne": "≠", "to": "→", "cdot": "·", "times": "×", "pm": "±",
    "ldots": "…", "prime": "′",
}  # fmt: skip


def svg_math(tex: str) -> str:
    """Convert a small TeX-like label to SVG text markup (already escaped).

    Supported: single Latin letters in italic (``x``, ``V``), words of two
    or more letters upright (``Stage``, ``max``), digits and punctuation
    upright, subscripts and superscripts ``x_1``, ``x_{ij}``, ``x^*``,
    ``c^{(k)}``, Greek letters and a few symbols written as ``\\pi``,
    ``\\lambda``, ``\\infty``, ``\\le``, ``\\ge``, ``\\to``, ``\\cdot``;
    bold upright ``\\mathbf{x}`` and the book's macros ``\\vx``, ``\\mA``.
    ``$`` signs are ignored, so ``"$x_1$"`` and ``"x_1"`` give the same.

    Example: ``svg_text(10, 20, svg_math("V(j) = 7"), raw=True)``.
    ``network_svg`` applies it to node, arc and note labels by default.
    """
    return "".join(_math_chunks(tex, scripts=True))


# Baseline offsets of subscripts and superscripts, in units of the normal
# font size (positive = down), and the size of scripts.
_SUB, _SUP, _SCRIPT_SCALE = 0.22, -0.42, 0.72


def _math_tokens(tex: str) -> list[tuple[str, str]]:
    """Split a label into (kind, text): kind is "word", "letter", "other",
    "sub" or "sup" (the latter two carry the raw script text)."""
    tokens: list[tuple[str, str]] = []
    i, n = 0, len(tex)

    def read_group(j: int) -> tuple[str, int]:
        # "{...}", a "\command", or a single character after _ or ^.
        if j < n and tex[j] == "{":
            depth, k = 1, j + 1
            while k < n and depth:
                depth += {"{": 1, "}": -1}.get(tex[k], 0)
                k += 1
            return tex[j + 1 : k - 1], k
        if j < n and tex[j] == "\\":
            k = j + 1
            while k < n and tex[k].isalpha():
                k += 1
            return tex[j:k], k
        return tex[j : j + 1], j + 1

    while i < n:
        ch = tex[i]
        if ch in "${}":
            i += 1
        elif ch in "_^":
            group, i = read_group(i + 1)
            tokens.append(("sub" if ch == "_" else "sup", group))
        elif ch == "\\":
            k = i + 1
            while k < n and tex[k].isalpha():
                k += 1
            word = tex[i + 1 : k]
            if word in ("mathbf", "boldsymbol"):
                # \mathbf{x}: the group is drawn bold and upright
                group, k = read_group(k)
                tokens.append(("bold", _GREEK.get(group.lstrip("\\"), group)))
            elif (
                word not in _GREEK
                and len(word) == 2
                and word[0] in "vm"
                and word[1].isalpha()
            ):
                # the book's macros \vx (bold vector), \mA (bold matrix)
                tokens.append(("bold", word[1]))
            elif word in _GREEK:
                symbol = _GREEK[word]
                tokens.append(
                    (
                        (
                            "letter"
                            if symbol.isalpha() and symbol.islower()
                            else "other"
                        ),
                        symbol,
                    )
                )
            elif word:
                tokens.append(("word", word))  # unknown command: its name
            else:  # escaped character such as "\{", "\%" or a space "\,"
                k = i + 2
                char = tex[i + 1 : k]
                tokens.append(
                    ("other", " " if char in (",", ";", " ", "!") else char)
                )
            i = k
        elif ch.isalpha():
            k = i
            while k < n and tex[k].isalpha():
                k += 1
            word = tex[i:k]
            tokens.append(("letter" if len(word) == 1 else "word", word))
            i = k
        else:
            k = i
            while k < n and not (tex[k].isalpha() or tex[k] in "_^\\${}"):
                k += 1
            tokens.append(("other", tex[i:k]))
            i = k
    return tokens


def _math_chunks(tex: str, *, scripts: bool) -> list[str]:
    """SVG markup pieces for ``svg_math``; scripts inside scripts are flat."""
    out: list[str] = []
    level = 0.0  # current baseline offset, in normal-size em
    for kind, text in _math_tokens(tex):
        if kind in ("sub", "sup") and scripts:
            target = _SUB if kind == "sub" else _SUP
            # dy is measured in the script's own (smaller) font size.
            dy = (target - level) / _SCRIPT_SCALE
            inner = "".join(_math_chunks(text, scripts=False))
            out.append(
                f'<tspan dy="{dy:.3f}em" font-size="{_SCRIPT_SCALE:.0%}">{inner}</tspan>'
            )
            level = target
            continue
        if kind in ("sub", "sup"):  # nested script: write it inline
            text = "".join(_math_chunks(text, scripts=False))
            out.append(text)
            continue
        piece = html.escape(text, quote=False)
        if kind == "letter":
            piece = f'<tspan font-style="italic">{piece}</tspan>'
        elif kind == "bold":
            piece = f'<tspan font-weight="bold">{piece}</tspan>'
        if level:
            # Back to the baseline after a sub/superscript.
            piece = f'<tspan dy="{-level:.3f}em">{piece}</tspan>'
            level = 0.0
        out.append(piece)
    return out


def _plain(tex: str) -> str:
    """The characters ``svg_math`` draws for ``tex`` (``\\infty`` -> ``∞``)."""
    return "".join(
        _plain(text) if kind in ("sub", "sup") else text
        for kind, text in _math_tokens(tex)
    )


def _text_width(text: str, size: float) -> float:
    """Rough width of a label in px (good enough for placement)."""
    # Measure the drawn glyphs, not the TeX source ("\\infty" is one).
    plain = _plain(text)
    # CJK characters (and full-width punctuation) are about one em wide.
    wide = sum(1 for ch in plain if ord(ch) >= 0x2E80)
    return (
        size * (0.6 * (len(plain) - wide) + 1.0 * wide)
        if plain
        else 0.6 * size
    )


def svg_text(
    x: float,
    y: float,
    text: str,
    *,
    size: float | None = None,
    anchor: str = "middle",
    color: str | None = None,
    weight: str | None = None,
    math: bool = True,
    halo: bool = False,
    raw: bool = False,
    rotate: float | None = None,
) -> str:
    """Return a ``<text>`` element centred (vertically) at ``(x, y)``.

    Parameters
    ----------
    x, y : float
        Anchor point in px.
    text : str
        The label. With ``math=True`` (default) it goes through
        ``svg_math`` (italic single letters, ``_`` subscripts); with
        ``math=False`` it is plain escaped text; with ``raw=True`` it is
        inserted as is (already SVG markup).
    size : float, optional
        Font size in px (default ``FONT_SIZE``).
    anchor : {"middle", "start", "end"}
        Horizontal alignment relative to ``x``.
    color : str, optional
        Fill colour (default ``INK``).
    weight : str, optional
        ``"bold"`` for bold text.
    halo : bool
        Draw a white outline behind the glyphs so the label stays legible
        where it crosses lines.
    rotate : float, optional
        Rotation in degrees about ``(x, y)`` (clockwise, SVG convention).
    """
    if raw:
        content = text
    elif math:
        content = svg_math(text)
    else:
        content = html.escape(text, quote=False)
    attrs = [
        f'x="{_num(x)}"',
        f'y="{_num(y)}"',
        f'text-anchor="{anchor}"',
        'dominant-baseline="central"',
    ]
    if size:
        attrs.append(f'font-size="{_num(size)}"')
    if color:
        attrs.append(f'fill="{color}"')
    if weight:
        attrs.append(f'font-weight="{weight}"')
    if halo:
        attrs.append(
            f'stroke="{WHITE}" stroke-width="4" stroke-linejoin="round" paint-order="stroke"'
        )
    if rotate:
        attrs.append(f'transform="rotate({_num(rotate)} {_num(x)} {_num(y)})"')
    return f"<text {' '.join(attrs)}>{content}</text>"


def svg_line(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    *,
    color: str = INK,
    width: float = 1.4,
    arrow: bool = False,
    highlight: bool = False,
    dashed: bool = False,
    shorten_start: float = 0.0,
    shorten_end: float = 0.0,
) -> str:
    """Return a straight ``<line>``, optionally with an arrow head at the end.

    With ``arrow=True`` the line is shortened by ``ARROW_TIP`` so that the
    arrow tip lands exactly on ``(x2, y2)`` (after ``shorten_end``).
    ``highlight=True`` uses ``HIGHLIGHT`` colour, a thicker stroke and the
    ``arrow-hl`` marker. ``shorten_start``/``shorten_end`` pull the ends in
    by that many px (e.g. a node radius).
    """
    length = _math.hypot(x2 - x1, y2 - y1) or 1.0
    ux, uy = (x2 - x1) / length, (y2 - y1) / length
    end_cut = shorten_end + (ARROW_TIP if arrow else 0.0)
    ax_, ay_ = x1 + ux * shorten_start, y1 + uy * shorten_start
    bx, by = x2 - ux * end_cut, y2 - uy * end_cut
    if highlight:
        color, width = HIGHLIGHT, max(width, 2.6)
    attrs = [
        f'x1="{_num(ax_)}" y1="{_num(ay_)}" x2="{_num(bx)}" y2="{_num(by)}"',
        f'stroke="{color}" stroke-width="{_num(width)}" stroke-linecap="butt"',
    ]
    if dashed:
        attrs.append('stroke-dasharray="5 4"')
    if arrow:
        marker = (
            "arrow-hl"
            if highlight
            else ("arrow-muted" if color == MUTED else "arrow")
        )
        attrs.append(f'marker-end="url(#{marker})"')
    return f"<line {' '.join(attrs)}/>"


# ---------------------------------------------------------------------------
# Networks
# ---------------------------------------------------------------------------

Arc = tuple  # (u, v) or (u, v, label)


def network_svg(
    nodes: Mapping[str, tuple[float, float]],
    arcs: Iterable[Arc],
    *,
    directed: bool = True,
    node_labels: Mapping[str, str] | None = None,
    node_radius: float = 16,
    node_size: float = 15,
    label_size: float = 13,
    highlight_arcs: Iterable[tuple[str, str]] = (),
    highlight_nodes: Iterable[str] = (),
    dashed_arcs: Iterable[tuple[str, str]] = (),
    muted_arcs: Iterable[tuple[str, str]] = (),
    bend: Mapping[tuple[str, str], float] | None = None,
    antiparallel_bend: float = 0.16,
    label_gap: float = 4,
    label_pos: Mapping[tuple[str, str], float] | None = None,
    label_side: Mapping[tuple[str, str], int] | None = None,
    label_shift: Mapping[tuple[str, str], tuple[float, float]] | None = None,
    node_notes: Mapping[str, str] | None = None,
    note_pos: Mapping[str, str | tuple[float, float]] | None = None,
    highlight_color: str = HIGHLIGHT,
    offset: tuple[float, float] = (0.0, 0.0),
    math: bool = True,
) -> str:
    """Draw a network (graph) as SVG elements; return them as a string.

    Pass the result to ``write_svg`` (alone, or together with other panels
    and labels). Arcs are drawn first, then nodes on top, then labels.

    Parameters
    ----------
    nodes : dict
        ``{name: (x, y)}`` in px, origin top left, y down. The order does
        not matter. Example: ``{"s": (30, 100), "1": (130, 40), ...}``.
    arcs : iterable
        ``(u, v)`` or ``(u, v, label)``; the label (str or number) is drawn
        next to the arc's midpoint, e.g. ``("s", "1", 4)`` or
        ``("1", "2", "[5, 2]")``.
    directed : bool
        Draw arrow heads at ``v``. If both ``(u, v)`` and ``(v, u)`` are
        arcs, the two are curved apart automatically (antiparallel pair).
        For undirected networks, ``highlight_arcs`` etc. match both orders.
    node_labels : dict, optional
        Text inside each circle (default: the node name), e.g.
        ``{"1": "v_1"}``. Use ``""`` for an empty circle.
    node_radius, node_size, label_size : float
        Circle radius, font size inside nodes, font size of arc labels (px).
    highlight_arcs, highlight_nodes : iterable
        Arcs ``(u, v)`` and nodes drawn thick in ``highlight_color`` (an
        optimal path, a spanning tree, an augmenting path).
    dashed_arcs, muted_arcs : iterable
        Arcs drawn dashed, or in grey (e.g. alternatives, removed arcs).
    bend : dict, optional
        ``{(u, v): b}`` curves arc ``(u, v)`` sideways by ``b`` times its
        length (``0.2`` is a gentle curve; negative bends the other way).
        Use it to route an arc around a node that lies on the straight line.
    antiparallel_bend : float
        Curvature used for antiparallel pairs.
    label_gap : float
        Space between an arc and its label, in px. The label is pushed
        clear of the line automatically, taking its width into account.
    label_pos : dict, optional
        ``{(u, v): t}`` puts the label at fraction ``t`` of the arc
        (default 0.5); move labels apart where arcs cross at the middle.
    label_side : dict, optional
        ``{(u, v): +1 or -1}`` flips the label to the other side of a
        straight arc (default: above it, or right of a vertical arc).
    label_shift : dict, optional
        ``{(u, v): (dx, dy)}`` final manual nudge of an arc label in px.
    node_notes : dict, optional
        ``{node: text}`` written next to a node: labels ``[V, P]``,
        supplies ``b_i``, DP values, e.g. ``{"3": "[7, 1]"}``.
    note_pos : dict, optional
        ``{node: "above" | "below" | "left" | "right" | (dx, dy)}``
        (default ``"above"``).
    highlight_color : str
        Colour of highlighted arcs and nodes (default ``HIGHLIGHT``).
        Arrow heads use the matching ``arrow-hl`` marker only for the
        default colour; other colours get ink arrow heads.
    offset : (float, float)
        Translate the whole drawing, to place several panels in one SVG.
    math : bool
        Pass labels through ``svg_math`` (italic letters, ``_`` subscripts).

    Returns
    -------
    str
        SVG elements wrapped in ``<g>``.

    Example
    -------
    ::

        nodes = {"s": (30, 90), "a": (150, 30), "b": (150, 150), "t": (270, 90)}
        arcs = [("s", "a", 4), ("s", "b", 2), ("a", "b", 1),
                ("b", "a", 3), ("a", "t", 5), ("b", "t", 6)]
        body = network_svg(nodes, arcs, highlight_arcs=[("s", "a"), ("a", "t")])
        return write_svg("04-03-small-network", body, 300, 180, outdir)
    """
    fmt = svg_math if math else (lambda s: html.escape(s, quote=False))

    # Node names may be given as str or int; compare them as strings.
    def arc_map(table) -> dict:
        return {(str(k[0]), str(k[1])): v for k, v in (table or {}).items()}

    def node_map(table) -> dict:
        return {str(k): v for k, v in (table or {}).items()}

    bend = arc_map(bend)
    label_pos = arc_map(label_pos)
    label_side = arc_map(label_side)
    label_shift = arc_map(label_shift)
    node_labels = node_map(node_labels)
    node_notes = node_map(node_notes)
    note_pos = node_map(note_pos)
    nodes = node_map(nodes)

    def key(u: str, v: str) -> tuple:
        # Undirected networks: (u, v) and (v, u) are the same edge.
        return (u, v) if directed else tuple(sorted((u, v)))

    def lookup(table: Mapping, u: str, v: str, default=None):
        if (u, v) in table:
            return table[(u, v)]
        if not directed and (v, u) in table:
            return table[(v, u)]
        return default

    highlight = {key(str(a[0]), str(a[1])) for a in highlight_arcs}
    dashed = {key(str(a[0]), str(a[1])) for a in dashed_arcs}
    muted = {key(str(a[0]), str(a[1])) for a in muted_arcs}
    hl_nodes = {str(n) for n in highlight_nodes}

    arc_list = [
        (str(a[0]), str(a[1]), a[2] if len(a) > 2 else None) for a in arcs
    ]
    arc_set = {(u, v) for u, v, _ in arc_list}
    for u, v, _ in arc_list:
        for end in (u, v):
            if end not in nodes:
                raise KeyError(
                    f"arc ({u}, {v}) uses node {end!r}, which is not in nodes"
                )

    arc_elems: list[str] = []
    label_elems: list[str] = []
    for u, v, label in arc_list:
        (x1, y1), (x2, y2) = nodes[u], nodes[v]
        dx, dy = x2 - x1, y2 - y1
        length = _math.hypot(dx, dy)
        if length == 0:
            continue  # self-loops are not drawn
        # Unit normal of the arc; "+1" side is above (or right of vertical).
        nx, ny = -dy / length, dx / length
        b = lookup(bend, u, v, 0.0)
        if not b and directed and (v, u) in arc_set:
            b = antiparallel_bend  # same relative side -> pair curves apart

        k = key(u, v)
        is_hl = k in highlight
        color = highlight_color if is_hl else (MUTED if k in muted else INK)
        width = 2.8 if is_hl else 1.4
        if is_hl and highlight_color == HIGHLIGHT:
            marker = "arrow-hl"
        elif k in muted:
            marker = "arrow-muted"
        else:
            marker = "arrow"
        tip = ARROW_TIP if directed else 0.0

        if b:
            # Quadratic Bezier with control point beside the midpoint.
            cx, cy = (x1 + x2) / 2 + nx * b * length, (
                y1 + y2
            ) / 2 + ny * b * length
            sx, sy = _toward(x1, y1, cx, cy, node_radius)
            ex, ey = _toward(x2, y2, cx, cy, node_radius + tip)
            d = f"M{_num(sx)},{_num(sy)} Q{_num(cx)},{_num(cy)} {_num(ex)},{_num(ey)}"
            t = lookup(label_pos, u, v, 0.5)
            # Point of the drawn curve at parameter t.
            px = (1 - t) ** 2 * sx + 2 * (1 - t) * t * cx + t**2 * ex
            py = (1 - t) ** 2 * sy + 2 * (1 - t) * t * cy + t**2 * ey
            side = 1 if b > 0 else -1  # label on the outer side of the bulge
            lnx, lny = nx * side, ny * side
        else:
            sx, sy = _toward(x1, y1, x2, y2, node_radius)
            ex, ey = _toward(x2, y2, x1, y1, node_radius + tip)
            d = f"M{_num(sx)},{_num(sy)} L{_num(ex)},{_num(ey)}"
            t = lookup(label_pos, u, v, 0.5)
            px, py = x1 + t * dx, y1 + t * dy
            # Default side: above the arc; right of a vertical arc.
            lnx, lny = nx, ny
            if lny > 1e-9 or (abs(lny) <= 1e-9 and lnx < 0):
                lnx, lny = -lnx, -lny
            lnx, lny = lnx * lookup(label_side, u, v, 1), lny * lookup(
                label_side, u, v, 1
            )

        attrs = f'fill="none" stroke="{color}" stroke-width="{width}"'
        if k in dashed:
            attrs += ' stroke-dasharray="5 4"'
        if directed:
            attrs += f' marker-end="url(#{marker})"'
        arc_elems.append(f'<path d="{d}" {attrs}/>')

        if label is not None and str(label) != "":
            text = str(label)
            # Push the label clear of the line: half its extent along the normal.
            half_w = _text_width(text, label_size) / 2
            half_h = label_size * 0.55
            dist = label_gap + abs(lnx) * half_w + abs(lny) * half_h
            lx, ly = px + lnx * dist, py + lny * dist
            sdx, sdy = lookup(label_shift, u, v, (0.0, 0.0))
            label_elems.append(
                svg_text(
                    lx + sdx,
                    ly + sdy,
                    fmt(text),
                    size=label_size,
                    color=highlight_color if is_hl else None,
                    halo=True,
                    raw=True,
                )
            )

    node_elems: list[str] = []
    for name, (x, y) in nodes.items():
        is_hl = name in hl_nodes
        stroke = highlight_color if is_hl else INK
        fill = HIGHLIGHT_FILL if is_hl else WHITE
        sw = 2.4 if is_hl else 1.4
        node_elems.append(
            f'<circle cx="{_num(x)}" cy="{_num(y)}" r="{_num(node_radius)}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'
        )
        text = node_labels.get(name, name)
        if text != "":
            node_elems.append(
                svg_text(x, y, fmt(str(text)), size=node_size, raw=True)
            )

        if name in node_notes:
            pos = note_pos.get(name, "above")
            gap = node_radius + 4
            note = str(node_notes[name])
            if not isinstance(pos, str):  # an offset (dx, dy), tuple or list
                dx_, dy_ = pos
                nx_, ny_, anchor = x + dx_, y + dy_, "middle"
            elif pos == "below":
                nx_, ny_, anchor = x, y + gap + label_size * 0.5, "middle"
            elif pos == "left":
                nx_, ny_, anchor = x - gap, y, "end"
            elif pos == "right":
                nx_, ny_, anchor = x + gap, y, "start"
            elif pos == "above":
                nx_, ny_, anchor = x, y - gap - label_size * 0.5, "middle"
            else:
                raise ValueError(
                    f"note_pos[{name!r}] must be 'above', 'below', 'left', "
                    f"'right' or (dx, dy), not {pos!r}"
                )
            node_elems.append(
                svg_text(
                    nx_,
                    ny_,
                    fmt(note),
                    size=label_size,
                    anchor=anchor,
                    halo=True,
                    raw=True,
                )
            )

    group = "\n".join(arc_elems + node_elems + label_elems)
    ox, oy = offset
    transform = (
        f' transform="translate({_num(ox)},{_num(oy)})"' if (ox or oy) else ""
    )
    return f"<g{transform}>\n{group}\n</g>"


def _toward(
    x: float, y: float, tx: float, ty: float, dist: float
) -> tuple[float, float]:
    """Move ``dist`` px from ``(x, y)`` toward ``(tx, ty)``."""
    length = _math.hypot(tx - x, ty - y) or 1.0
    return x + (tx - x) / length * dist, y + (ty - y) / length * dist


# ---------------------------------------------------------------------------
# Two-variable linear programs (matplotlib)
# ---------------------------------------------------------------------------

Number = int | float | Fraction


def _frac(value: Number) -> Fraction:
    """Exact rational for ints/Fractions; floats are rounded sensibly."""
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    return Fraction(value).limit_denominator(10**6)


def _line_in_box(a1: Fraction, a2: Fraction, b: Fraction, xlim, ylim):
    """Endpoints of the segment of a1*x1 + a2*x2 = b inside the box, or None."""
    (x0, x1), (y0, y1) = (tuple(map(_frac, xlim)), tuple(map(_frac, ylim)))
    points = []
    if a2 != 0:  # intersections with the vertical sides
        for x in (x0, x1):
            y = (b - a1 * x) / a2
            if y0 <= y <= y1:
                points.append((x, y))
    if a1 != 0:  # intersections with the horizontal sides
        for y in (y0, y1):
            x = (b - a2 * y) / a1
            if x0 <= x <= x1:
                points.append((x, y))
    unique = sorted(set(points))
    if len(unique) < 2:
        return None
    # The two extreme points along the line.
    return unique[0], unique[-1]


def lin_label(a1: Number, a2: Number, b: Number, op: str = "=") -> str:
    """Mathtext label of ``a1*x1 + a2*x2 op b``, e.g. ``$6x_1 + 5x_2 = 60$``."""

    def number(value: Number) -> str:
        # Integers plainly, fractions as \frac, floats as written (0.15).
        if isinstance(value, float):
            return f"{value:g}"
        value = _frac(value)
        if value.denominator == 1:
            return str(value.numerator)
        return rf"\frac{{{value.numerator}}}{{{value.denominator}}}"

    def term(c: Number, var: str, first: bool) -> str:
        if c == 0:
            return ""
        negative = c < 0
        mag = -c if negative else c
        text = "" if mag == 1 else number(mag)
        if first:
            return f"{'-' if negative else ''}{text}{var}"
        return f" {'-' if negative else '+'} {text}{var}"

    lhs = term(a1, "x_1", True) + term(a2, "x_2", a1 == 0)
    rel = {"<=": r"\leq", ">=": r"\geq", "=": "="}[op]
    if b < 0:
        rhs = "-" + number(-b)
    else:
        rhs = number(b)
    return f"${lhs.strip()} {rel} {rhs}$"


def draw_line(
    ax,
    a1: Number,
    a2: Number,
    b: Number,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    *,
    label: str | None = None,
    label_at: float = 0.8,
    label_color: str | None = None,
    label_size: float | None = None,
    **plot_kw,
):
    """Draw the line ``a1*x1 + a2*x2 = b`` across the window; label it.

    Use it for constraint lines and objective contours (iso-profit lines).

    Parameters
    ----------
    ax : matplotlib Axes
    a1, a2, b : number
        Line coefficients.
    xlim, ylim : (float, float)
        Window in which the line is drawn.
    label : str, optional
        Text written along the line, above it (mathtext allowed).
    label_at : float
        Position of the label along the visible segment (0 = left/bottom
        end, 1 = right/top end).
    **plot_kw
        Passed to ``ax.plot`` (``color``, ``linestyle``, ``linewidth``...).

    Returns
    -------
    list of matplotlib artists (line and label), or [] if the line misses
    the window.
    """
    seg = _line_in_box(_frac(a1), _frac(a2), _frac(b), xlim, ylim)
    if seg is None:
        return []
    (px, py), (qx, qy) = seg
    px, py, qx, qy = float(px), float(py), float(qx), float(qy)
    plot_kw.setdefault("linewidth", 1.5)
    (line,) = ax.plot([px, qx], [py, qy], **plot_kw)
    artists = [line]
    if label:
        tx, ty = px + label_at * (qx - px), py + label_at * (qy - py)
        # Angle in data coordinates; matplotlib converts it to screen angle.
        angle = _math.degrees(_math.atan2(qy - py, qx - px))
        if angle > 90:
            angle -= 180
        elif angle <= -90:
            angle += 180
        text = ax.text(
            tx,
            ty,
            label,
            rotation=angle,
            rotation_mode="anchor",
            transform_rotates_text=True,
            ha="center",
            va="bottom",
            color=label_color or line.get_color(),
            fontsize=label_size,
        )
        artists.append(text)
    return artists


def _clip(
    polygon: list, a1: Fraction, a2: Fraction, b: Fraction, op: str
) -> list:
    """Sutherland-Hodgman: keep the part of ``polygon`` with a.x op b."""

    def inside(p) -> bool:
        value = a1 * p[0] + a2 * p[1]
        return value <= b if op == "<=" else value >= b

    def cross(p, q):
        # Intersection of segment pq with the line a.x = b (exact).
        fp, fq = a1 * p[0] + a2 * p[1] - b, a1 * q[0] + a2 * q[1] - b
        t = fp / (fp - fq)
        return (p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]))

    result = []
    for i, current in enumerate(polygon):
        previous = polygon[i - 1]
        if inside(current):
            if not inside(previous):
                result.append(cross(previous, current))
            result.append(current)
        elif inside(previous):
            result.append(cross(previous, current))
    # Remove consecutive duplicates.
    cleaned: list = []
    for p in result:
        if not cleaned or p != cleaned[-1]:
            cleaned.append(p)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    return cleaned


def _normalise_constraint(spec) -> dict:
    """Accept ``(a1, a2, op, b)``, ``(a1, a2, op, b, label)`` or a dict."""
    if isinstance(spec, Mapping):
        c = dict(spec)
        a1, a2 = c.pop("a")
        c.update(a1=a1, a2=a2)
    else:
        a1, a2, op, b, *rest = spec
        c = {"a1": a1, "a2": a2, "op": op, "b": b}
        if rest:
            c["label"] = rest[0]
    if c["op"] not in ("<=", ">=", "="):
        raise ValueError(
            f"constraint operator must be '<=', '>=' or '=', not {c['op']!r}"
        )
    return c


def lp_region(
    ax,
    constraints: Sequence,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    *,
    nonneg: bool = True,
    fill: str = BLUE,
    alpha: float = 0.15,
    colors: Sequence[str] | None = None,
    labels: bool = True,
    label_at: float = 0.8,
    label_size: float | None = None,
    mark_vertices: bool = False,
    axis_labels: tuple[str, str] | None = ("$x_1$", "$x_2$"),
    exact: bool = False,
) -> list[tuple]:
    """Shade the feasible region of a two-variable LP; return its vertices.

    Draws every constraint line across the window (coloured from
    ``PALETTE``, labelled with its equation or your label), shades the
    feasible region, sets the axis limits and labels, and returns the
    extreme points of the region.

    Parameters
    ----------
    ax : matplotlib Axes
    constraints : sequence
        One entry per constraint ``a1*x1 + a2*x2 (op) b``, with op one of
        ``"<="``, ``">="``, ``"="``:

        * ``(a1, a2, op, b)`` - labelled with its equation, e.g.
          ``(6, 5, "<=", 60)`` -> "6x1 + 5x2 = 60";
        * ``(a1, a2, op, b, label)`` - your own label (``None`` for none);
        * a dict ``{"a": (a1, a2), "op": "<=", "b": 60, "label": ...,
          "label_at": 0.3, "color": ORANGE, "linestyle": "--"}`` for
          per-constraint options.

        Use integers or ``Fraction`` for exact vertices.
    xlim, ylim : (float, float)
        The window. An unbounded region is shaded up to the window edge.
    nonneg : bool
        Add ``x1 >= 0`` and ``x2 >= 0`` (not drawn: they are the axes).
    fill, alpha : str, float
        Colour and opacity of the shaded region.
    colors : sequence of str, optional
        Line colours in constraint order (default ``PALETTE``).
    labels : bool
        Write labels on the lines.
    label_at : float
        Default label position along each visible line (0..1); move single
        labels with the dict form (``"label_at"``) when they collide.
    mark_vertices : bool
        Draw a small dot at every vertex.
    axis_labels : (str, str) or None
        Axis titles (default ``$x_1$``, ``$x_2$``).
    exact : bool
        Return vertices as ``Fraction`` instead of float.

    Returns
    -------
    list of (x1, x2)
        The vertices (extreme points) of the feasible region, in
        counter-clockwise order. Corners created only by the window edge
        of an unbounded region are not vertices and are left out. Empty if
        the region is empty inside the window.

    Example
    -------
    ::

        fig, ax = plt.subplots(figsize=(5.5, 4.5))
        vertices = lp_region(
            ax,
            [(6, 5, "<=", 60), (10, 20, "<=", 150), (1, 0, "<=", 8)],
            xlim=(0, 11), ylim=(0, 9), mark_vertices=True,
        )
        # vertices == [(0, 0), (8, 0), (8, 2.4), (6.43, 4.29), (0, 7.5)]
        # Objective contour 500 x1 + 450 x2 = 36000/7 through the optimum.
        draw_line(ax, 500, 450, Fraction(36000, 7), (0, 11), (0, 9),
                  color=MUTED, linestyle="--", label="$z = 36000/7$")
        return save_mpl(fig, "03-02-product-mix-region", outdir)
    """
    specs = [_normalise_constraint(c) for c in constraints]
    colors = list(colors or PALETTE)
    (x0, x1), (y0, y1) = (tuple(map(_frac, xlim)), tuple(map(_frac, ylim)))

    # Half-planes, as (a1, a2, op, b) with exact coefficients.
    halfplanes = []
    for c in specs:
        a1, a2, b = _frac(c["a1"]), _frac(c["a2"]), _frac(c["b"])
        if c["op"] == "=":
            halfplanes += [(a1, a2, "<=", b), (a1, a2, ">=", b)]
        else:
            halfplanes.append((a1, a2, c["op"], b))
    boundary_lines = [
        (_frac(c["a1"]), _frac(c["a2"]), _frac(c["b"])) for c in specs
    ]
    if nonneg:
        halfplanes += [(Fraction(1), Fraction(0), ">=", Fraction(0)),
                       (Fraction(0), Fraction(1), ">=", Fraction(0))]  # fmt: skip
        boundary_lines += [(Fraction(1), Fraction(0), Fraction(0)),
                           (Fraction(0), Fraction(1), Fraction(0))]  # fmt: skip

    # Start from the window and cut it by every half-plane (exact arithmetic).
    polygon = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    for a1, a2, op, b in halfplanes:
        if a1 == 0 and a2 == 0:
            continue
        polygon = _clip(polygon, a1, a2, b, op)
        if not polygon:
            break

    # Shade the region (a segment or point if it has no area).
    if polygon:
        xs = [float(p[0]) for p in polygon]
        ys = [float(p[1]) for p in polygon]
        area = (
            sum(xs[i - 1] * ys[i] - xs[i] * ys[i - 1] for i in range(len(xs)))
            / 2
        )
        if abs(area) > 1e-12:
            ax.fill(xs, ys, color=fill, alpha=alpha, linewidth=0, zorder=1)
        else:
            ax.plot(
                xs + xs[:1],
                ys + ys[:1],
                color=fill,
                linewidth=4,
                alpha=0.6,
                zorder=1,
            )

    # Constraint lines and labels.
    for i, c in enumerate(specs):
        style = {
            "color": c.get("color", colors[i % len(colors)]),
            "linestyle": c.get("linestyle", "-"),
            "zorder": 2,
        }
        if "linewidth" in c:
            style["linewidth"] = c["linewidth"]
        label = (
            c.get("label", lin_label(c["a1"], c["a2"], c["b"]))
            if labels
            else None
        )
        draw_line(
            ax,
            c["a1"],
            c["a2"],
            c["b"],
            xlim,
            ylim,
            label=label,
            label_at=c.get("label_at", label_at),
            label_size=label_size,
            **style,
        )

    # Vertices: polygon corners on at least two distinct constraint lines.
    def normalised(line):
        # Scale so the first non-zero coefficient is 1 (duplicates collapse).
        a1, a2, b = line
        lead = a1 if a1 != 0 else a2
        return (a1 / lead, a2 / lead, b / lead)

    lines = {
        normalised(ln) for ln in boundary_lines if ln[0] != 0 or ln[1] != 0
    }

    def lines_through(p) -> int:
        return sum(a1 * p[0] + a2 * p[1] == b for a1, a2, b in lines)

    vertices = [p for p in polygon if lines_through(p) >= 2]
    # Counter-clockwise order (the clipping keeps the window's orientation),
    # starting from the lexicographically smallest vertex.
    if vertices:
        start = vertices.index(min(vertices))
        vertices = vertices[start:] + vertices[:start]
    if mark_vertices and vertices:
        ax.plot(
            [float(p[0]) for p in vertices],
            [float(p[1]) for p in vertices],
            linestyle="none",
            marker="o",
            markersize=4.5,
            color=INK,
            zorder=4,
            clip_on=False,  # dots on the axes stay round
        )

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    if axis_labels:
        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])
    if exact:
        return vertices
    return [(float(p[0]), float(p[1])) for p in vertices]
