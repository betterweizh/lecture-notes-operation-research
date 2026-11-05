# Figure generators

Every generated figure of the book is drawn by a Python function in this
folder. This page is the guide for writing one.

```text
scripts/figures/
  common.py      shared style and helpers (do not edit; ask the maintainer)
  pPPcCC.py      figures of part PP, chapter CC (one module per chapter)
  _example.py    a template module; NOT part of the book (see below)
  preview.py     renders figures to PNG for a quick look
scripts/make_figures.py   draws the figures into docs/<lang>/figures/
```

## The rules

- **One module per chapter**: `scripts/figures/pPPcCC.py`, where `PP` and
  `CC` are the part and chapter numbers of the chapter that first shows
  the figure (`p03c02.py` for chapter 3.2). Only files named like this are
  loaded.
- **One function per figure**, `func(outdir, lang="en")`. It writes exactly
  one file, `<name>.svg`, into `outdir` (with `save_mpl` or `write_svg`)
  and returns its path. `lang` is the edition; use it only when a label
  has words that need translating (prefer symbols and numbers).
- **Register every function** in the module's `FIGURES` dict, at the end of
  the module: `FIGURES = {"03-02-feasible-region": feasible_region, ...}`.
  The name is the file name without `.svg`, starts with the part and
  chapter numbers (`03-02-`), and is unique in the whole book
  (`make_figures.py` stops with an error on a duplicate name).
- **Comment the code**: a docstring saying what the figure shows and where
  it is used, comments on the data (constraint lines, coordinates, arc
  lists) and on each drawing step.
- Use the **book's notation** in every label, and the style and helpers of
  `common.py`. No random data without a fixed seed: regenerated figures
  must not change.
- Generated files in `docs/<lang>/figures/` are build products (git-ignored).
  Never edit them by hand.

In the chapter, cite the figure by its file name and label:

````markdown
```{figure} figures/03-02-feasible-region.svg
:name: fig-03-02-feasible-region
:width: 520px
:class: dark-light
:alt: Feasible region of the product-mix LP with the optimal vertex marked

Caption that explains what is shown and what to notice.
```
````

## Commands

Figures need numpy and matplotlib (`/usr/local/bin/python3` in the
sandbox; any interpreter with `requirements.txt` installed elsewhere).

```bash
python3 scripts/make_figures.py 03-02-feasible-region       # draw one (or several) figures
python3 scripts/make_figures.py                          # draw all figures
python3 scripts/make_figures.py --list                   # all names and their modules
python3 scripts/figures/preview.py 03-02-feasible-region    # -> /tmp/fig-preview/03-02-feasible-region.png
python3 scripts/figures/preview.py --make --dark --width 520 03-02-feasible-region
```

`preview.py` renders the SVG with Chromium (Playwright) onto a white page,
2x resolution. `--make` draws the figure first, `--width` shows it at the
`:width:` of the chapter, `--dark` adds `<name>.dark.png` as the website
shows it in dark mode. Open the PNG with an image viewer (or the Read
tool) and check: no overlapping labels, arrows point the right way, the
drawing matches the text and the data.

`scripts/build_docs.py` draws the figures that chapters reference when
they are missing or older than their module (or `common.py`); `build_docs.py
--figures-only` does only that. A reference to a figure that has no
generator (and is not in `docs/_shared/figures/`) stops the build with
the chapter, line and file name.

## Start from the example

`_example.py` shows a matplotlib figure (`lp_region`, `draw_line`) and a
network (`network_svg`). The leading underscore keeps it out of the book.
To look at its output:

```bash
python3 scripts/make_figures.py --example --outdir /tmp/fig-preview/example
python3 scripts/figures/preview.py /tmp/fig-preview/example/00-example-lp.svg \
    /tmp/fig-preview/example/00-example-network.svg
```

To start your module, copy it to `pPPcCC.py`, rename the functions and the
`FIGURES` names, and replace the data.

## Helpers in `common.py`

Import what you need with a relative import (style is applied on import):

```python
import matplotlib.pyplot as plt
from .common import BLUE, ORANGE, MUTED, lp_region, draw_line, save_mpl
from .common import network_svg, svg_text, svg_line, write_svg
```

### Style constants

| Name | Value | Use |
|---|---|---|
| `INK` | `#2b2b2b` | text, axes, lines, node outlines |
| `MUTED` | `#6b6b6b` | secondary text and lines |
| `GRID` | `#d9d9d9` | grid lines |
| `PALETTE` | `#2a78d6`, `#eb6834`, `#1baf7a`, then purple, amber, pink, teal, grey | series colours, in this order |
| `BLUE`, `ORANGE`, `GREEN`, … | the palette entries | |
| `HIGHLIGHT` | `#2a78d6` (blue) | highlighted arcs and nodes |
| `FONT_FAMILY`, `FONT_SIZE` | DejaVu Sans stack, 14 px | hand-written SVG text |

Matplotlib is set up with `svg.fonttype: none` (text stays text),
transparent background, no top/right spines, light grid, the palette as the
colour cycle, and font size 11. A figure of `figsize=(5.5, 4)` displays well
at `:width: 520px`.

### Saving

- `save_mpl(fig, name, outdir)` saves a matplotlib figure as
  `outdir/<name>.svg` (tight bounding box, no date stamp) and closes it.
- `write_svg(name, body, w, h, outdir, *, extra_defs="")` wraps SVG elements
  in a complete `w` x `h` px file with the house font and the arrow markers
  `url(#arrow)`, `url(#arrow-hl)` (highlight colour) and `url(#arrow-muted)`.

### Hand-written SVG

Coordinates are pixels, origin top left, **y down**.

- `network_svg(nodes, arcs, *, directed=True, node_labels=None,
  node_radius=16, highlight_arcs=(), highlight_nodes=(), dashed_arcs=(),
  muted_arcs=(), bend=None, label_pos=None, label_side=None,
  label_shift=None, node_notes=None, note_pos=None, offset=(0, 0), ...)`
  returns the SVG of a network. `nodes = {name: (x, y)}`; `arcs` are
  `(u, v)` or `(u, v, label)`. Node names inside the circles, arc labels
  placed beside the arc's midpoint clear of the line, arrow heads on the
  circle boundary. Antiparallel pairs `(u, v)`, `(v, u)` are curved apart
  automatically. Undirected networks: `directed=False`.
  - Fix a label collision with `label_pos={(u, v): 0.35}` (position along the
    arc), `label_side={(u, v): -1}` (other side) or
    `label_shift={(u, v): (dx, dy)}`.
  - Route an arc around a node in its way with `bend={(u, v): 0.25}`.
  - Labels next to nodes (`[V, P]`, supplies, DP values): `node_notes` and
    `note_pos` (`"above"`, `"below"`, `"left"`, `"right"` or `(dx, dy)`).
  - Several panels in one figure: call it once per panel with
    `offset=(dx, dy)` and add titles with `svg_text`.
- `svg_text(x, y, text, *, size=None, anchor="middle", color=None,
  weight=None, math=True, halo=False, raw=False, rotate=None)`: a label
  centred vertically at `(x, y)`.
- `svg_line(x1, y1, x2, y2, *, color=INK, width=1.4, arrow=False,
  highlight=False, dashed=False, shorten_start=0, shorten_end=0)`: a line;
  with `arrow=True` the tip lands exactly on `(x2, y2)`.
- `svg_math(tex)`: TeX-like label to SVG text: single letters italic,
  `x_1`, `x_{ij}`, `x^*`, Greek `\pi`, `\lambda`, symbols `\le`, `\ge`,
  `\infty`, `\to`. `network_svg` and `svg_text` apply it by default
  (`math=False` for plain words).

### Two-variable LPs (matplotlib)

- `lp_region(ax, constraints, xlim, ylim, *, nonneg=True, fill=BLUE,
  alpha=0.15, colors=None, labels=True, label_at=0.8, label_size=None,
  mark_vertices=False, axis_labels=("$x_1$", "$x_2$"), exact=False)`
  shades the feasible region, draws and labels every constraint line, sets
  limits and axis titles, and **returns the vertices** counter-clockwise
  (window corners of an unbounded region are not vertices). A constraint
  is `(a1, a2, op, b)` with `op` one of `"<="`, `">="`, `"="`, or
  `(a1, a2, op, b, label)`, or a dict
  `{"a": (a1, a2), "op": "<=", "b": 60, "label": ..., "label_at": 0.3,
  "color": ORANGE, "linestyle": "--"}`. With integer data and
  `exact=True` the vertices are exact `Fraction`s.
- `draw_line(ax, a1, a2, b, xlim, ylim, *, label=None, label_at=0.8,
  **plot_kw)` draws `a1*x1 + a2*x2 = b` across the window with a label
  along it: objective contours, extra lines.
- `lin_label(a1, a2, b, op="=")` returns the mathtext label
  `"$6x_1 + 5x_2 = 60$"`.

## Matplotlib labels

The MathJax macros of the chapters (`\vx`, `\R`, `\T`, ...) do not exist in
matplotlib. Write them out in mathtext: `$\mathbf{x}^\mathsf{T}$`,
`$\boldsymbol{\pi}$`, `$\mathbb{R}^n$`, `$\operatorname{conv}$`,
`$\mathrm{s.t.}$`.
