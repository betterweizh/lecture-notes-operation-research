# Lecture Notes on Operations Research · 运筹学学习笔记

**English** | [中文](README.zh-CN.md) | [Read online](https://betterweizh.github.io/lecture-notes-operation-research/en/) | [在线阅读](https://betterweizh.github.io/lecture-notes-operation-research/zh)

Lecture notes on operations research: linear programming, duality and
sensitivity analysis, network optimization, dynamic programming, nonlinear
programming and constraint linearization, with a preliminary part on linear
algebra, matrix calculus and convexity. Every method is derived from first
principles and illustrated by about a hundred worked examples; more than
160 further examples and problems come with complete, folded solutions; all
numbers were recomputed in code; and a final part implements the algorithms in
Python (NumPy, SciPy/HiGHS, PuLP, NetworkX, CVXPY).

The book is written in MyST Markdown, built with Sphinx and published with
GitHub Pages, in English and in Simplified Chinese.

## Language switch

| Edition | Source | Built site |
|---|---|---|
| English | `docs/en/` | `_build/html/en/` |
| 简体中文 | `docs/zh/` | `_build/html/zh/` |

The translate icon at the top right of every page opens the same chapter
in the other edition, and the root page `_build/html/index.html` remembers
the reader's choice. Both editions use identical file names, labels,
equations, code and figure names; figures are drawn for each edition with
its own labels. `scripts/glossary.py` is the single source of the
English–Chinese terminology and writes `glossary.md` in both editions.

## Contents

| Part | Chapters |
|---|---|
| **I — Preliminaries** | 1.1 Linear Algebra Essentials · 1.2 Vector and Matrix Calculus · 1.3 Introduction to Operations Research |
| **II — Convex Optimization** | 2.1 Convex Sets · 2.2 Convex Functions · 2.3 Convex Optimization and Lagrange Duality |
| **III — Linear Programming** | 3.1 Linear Programming: Model Building · 3.2 Geometry of Linear Programs · 3.3 The Simplex Method · 3.4 Artificial Variables: Two-Phase and Big-M Methods · 3.5 Duality Theory · 3.6 Revised and Dual Simplex Methods · 3.7 Sensitivity Analysis |
| **IV — Network Optimization** | 4.1 Networks and Minimum Spanning Trees · 4.2 Shortest-Path Problems · 4.3 Maximum Flow · 4.4 Minimum-Cost Flow |
| **V — Dynamic Programming** | 5.1 Dynamic Programming Fundamentals · 5.2 Equipment Replacement · 5.3 Resource Allocation and Capital Budgeting · 5.4 Shortest Paths by Dynamic Programming |
| **VI — Nonlinear Programming** | 6.1 Nonlinear Programming Models · 6.2 Optimality Conditions · 6.3 Karush–Kuhn–Tucker Conditions · 6.4 Line Search Methods · 6.5 Algorithms for Unconstrained Optimization · 6.6 Algorithms for Constrained Optimization |
| **VII — Linearization** | 7.1 Logical Constraints and the Big-M Method · 7.2 Linearization Techniques · 7.3 Piecewise-Linear Approximation |
| **VIII — Examples and Problems** | 8.1 LP Modeling and Geometry · 8.2 The Simplex Method · 8.3 Duality · 8.4 Revised and Dual Simplex, Sensitivity · 8.5 Network Optimization · 8.6 DP Fundamentals and Shortest Paths · 8.7 DP Applications · 8.8 Optimality Conditions and KKT · 8.9 Nonlinear Programming Algorithms |
| **IX — Implementation** | 9.1 Software Setup and Solvers · 9.2 Modeling LPs and MILPs in Python · 9.3 Linear Programming Algorithms in Python · 9.4 Network Algorithms in Python · 9.5 Dynamic Programming in Python · 9.6 Nonlinear Programming in Python |
| **Appendix** | Notation · References |

## Repository layout

```text
docs/
  _shared/            configuration, extensions, templates, CSS shared by all editions
    conf_common.py    Sphinx settings, MathJax macros, theme
    _ext/             part-prefixed chapter numbering, language switch
    figures/          existing (non-generated) figures, with SOURCES.md
    locales/          Chinese names of the proof blocks (定义, 定理, 例, …)
  en/                 English edition: conf.py, index.md, PP-CC-*.md chapters (part, chapter),
                      notation.md, glossary.md, references.md
  zh/                 Simplified Chinese edition: the same files, translated
scripts/
  build_docs.py       builds the website (figures first)
  make_figures.py     draws the generated figures
  figures/            one figure module per chapter (pPPcCC.py) + common.py
  run_code_cells.py   runs the Python blocks of the code chapters
  glossary.py         English–Chinese terminology; writes glossary.md of both editions
.github/workflows/    build and deploy to GitHub Pages
requirements.txt      Sphinx toolchain and the packages of the code chapters
```

`sources/` (the course materials the notes are based on) is kept locally
and is not part of the repository (it is listed in `.gitignore`).

## Figures

All figures are generated: each is drawn by a commented Python function in
`scripts/figures/pPPcCC.py` (part `PP`, chapter `CC` of the chapter that
first shows it), registered under a name that starts with the part and
chapter numbers (`03-02-product-mix-region`), and written to `docs/en/figures/` when the site
is built. Generated files are git-ignored; never edit them by hand.
Existing (non-generated) figures would go to `docs/_shared/figures/` and be
credited in its `SOURCES.md`. See
[`scripts/figures/README.md`](scripts/figures/README.md) for the helpers.

## Build locally

```bash
python -m venv .venv
source .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install -r requirements.txt

python scripts/build_docs.py --serve        # build and open http://127.0.0.1:8000/
python scripts/build_docs.py --strict       # warnings are errors (as on GitHub)
python scripts/build_docs.py --clean        # rebuild from scratch
python scripts/build_docs.py --figures-only # only draw the figures
python scripts/build_docs.py --lang zh       # build one edition only
python scripts/run_code_cells.py docs/en/09-*.md   # run all code blocks
```

The site goes to `_build/html/`. Figures are drawn automatically when they
are missing or older than their generator.

## Editing workflow

1. Edit the chapter in `docs/en/` and the same-named file in `docs/zh/`
   (and the figure module in `scripts/figures/` if a figure changes; a
   label with words has an English and a Chinese version). Keep labels
   (`eq-03-02-…`, `thm-03-05-…`, `fig-04-03-…`) stable and identical in
   both editions: other chapters refer to them.
2. A new chapter needs an entry in the matching toctree of both
   `docs/en/index.md` and `docs/zh/index.md`; a new term goes into
   `scripts/glossary.py` (then run it).
3. Run `python scripts/build_docs.py --strict` and, for code chapters,
   `scripts/run_code_cells.py`.
4. Push to `main`: GitHub Actions builds the site and publishes it to
   GitHub Pages (enable Pages with source "GitHub Actions" once).

## Writing conventions

- Vectors are bold lower case ($\mathbf{x}$), matrices bold upper case
  ($\mathbf{A}$), written with the MathJax macros of
  `docs/_shared/conf_common.py` (`\vx`, `\mA`, `\T`, …); the notation page
  lists every symbol.
- LP tableaux put the objective row last, labelled $-z$, with reduced costs
  $\bar c_j$; dual variables are $\boldsymbol{\pi}$.
- Chapters follow one skeleton: Overview, Learning objectives, sections,
  Summary and Common Mistakes, References. Examples use `{prf:example}`
  followed by a folded `{dropdown} Solution`; proofs are folded too.
- Chapters are numbered within their part (3.2 = Part III, chapter 2);
  figures, tables, examples and equations carry that prefix (3.2.1).
- Labels carry the chapter number; cross-chapter links use `{doc}` and
  labels, never heading anchors.
- References follow APA 7; each chapter lists what it cites.

## Author and license

Better W.Z. — released under the [MIT License](LICENSE). The course
materials the notes are based on and all cited works remain the property
of their owners.

## How to cite

Better W.Z. (2026). *Lecture Notes on Operations Research*. GitHub.
https://github.com/betterweizh/lecture-notes-operation-research

```bibtex
@misc{betterwz2026lecture,
  author       = {Better W.Z.},
  title        = {Lecture Notes on Operations Research},
  year         = {2026},
  howpublished = {\url{https://github.com/betterweizh/lecture-notes-operation-research}},
  note         = {Lecture notes}
}
```
