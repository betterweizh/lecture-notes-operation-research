# Lecture Notes on Operations Research

These notes develop the core models and methods of operations research: linear programming together with its duality theory and sensitivity analysis, network optimization, dynamic programming, nonlinear programming, and linearization techniques that convert logical and nonlinear relations into mixed-integer linear models. For each method, the assumptions are stated first, and the method is then derived from first principles. The linear algebra and matrix calculus needed are reviewed in a preliminaries section, while convex sets, convex functions, and Lagrangian duality are treated in a separate chapter.

Some methods are accompanied by worked examples and exercises with complete solutions, and all numerical results have been recomputed and verified with code. The final part implements the algorithms in Python and demonstrates how to solve the same models with SciPy (HiGHS), PuLP, NetworkX, and CVXPY.

This study note is based on lecture notes from the course *CE6001 Operations and Management of Infrastructure Systems* at the National University of Singapore (lecturer: Prof. MENG Qiang; the author took the course in the 2021/22 academic year). It is further supplemented by the author’s personal study notes and classic reference books, in particular Taha’s *Operations Research: An Introduction*, Boyd and Vandenberghe’s Convex Optimization, and Bazaraa, Sherali, and Shetty’s Nonlinear Programming. All material has been reorganized into a single study note with unified notation.

In November 2026, the content was reviewed, revised, and organized with the assistance of Claude Opus 5.5 (Anthropic), translated into Simplified Chinese, and a website was built for it.

The author’s knowledge is limited, and this note inevitably contains omissions and errors. Corrections from readers are sincerely welcome. If you find any problems, please submit feedback via a [GitHub issue](https://github.com/betterweizh/lecture-notes-operation-research/issues) to help improve these notes.

:::{admonition} 中文版
:class: tip

A complete Simplified Chinese edition is available. Click the translate icon
{fas}`language` at the top right of any page to open the same chapter in the
other language.
:::

## How the book is organized

- **Part I — Preliminaries** reviews the linear algebra (linear systems,
  pivoting, eigenvalues, definiteness) and the vector and matrix calculus
  (gradients, Hessians, Taylor expansions) used throughout, and introduces
  operations research: its history, the modeling process and the problem
  classes of this book.
- **Part II — Convex Optimization** covers convex sets, convex functions,
  convex optimization problems and Lagrange duality.
- **Part III — Linear Programming** develops model building, the geometry
  of linear programs, the simplex method with artificial variables,
  duality theory, the revised and dual simplex methods and sensitivity
  analysis.
- **Part IV — Network Optimization** treats minimum spanning trees,
  shortest paths, maximum flow and minimum-cost flow.
- **Part V — Dynamic Programming** presents the principle of optimality and
  its applications to equipment replacement, resource allocation and
  shortest paths.
- **Part VI — Nonlinear Programming** covers models, optimality conditions
  and the Karush–Kuhn–Tucker (KKT) conditions, line search, and algorithms for unconstrained and
  constrained problems.
- **Part VII — Linearization** shows how logical conditions, products and
  piecewise-linear functions become mixed-integer linear constraints.
- **Part VIII — Examples and Problems** works through the examples and
  problems of Parts III–VI topic by topic.
- **Part IX — Implementation** sets up the software, compares solvers and
  modeling tools, and implements the algorithms in Python, each checked
  against a library solver.

Chapters are numbered within their part: chapter 3.2 is the second chapter
of Part III, and its figures, tables, examples and equations are numbered
3.2.1, 3.2.2, … . Proofs and the solutions of examples and problems are
folded: click a **Proof** or **Solution** bar to open it, and try each
problem before you look.

The appendix lists the {doc}`notation` used across chapters, a
{doc}`glossary` of English and Chinese terms, and all {doc}`references`.

## License

The text, figures and code of these notes are released under the MIT
License: you may use, copy, modify, merge, publish, distribute, sublicense
and sell copies, provided that the copyright notice "Copyright (c) 2026
Better W.Z." and the license text are kept in all copies or substantial
portions. The notes are provided "as is", without warranty of any kind. The
course materials they are based on, and all cited works, remain the property
of their owners.

## How to cite

Better W.Z. (2026). *Lecture Notes on Operations Research*.
GitHub. https://github.com/betterweizh/lecture-notes-operation-research

```bibtex
@misc{betterwz2026lecture,
  author       = {Better W.Z.},
  title        = {Lecture Notes on Operations Research},
  year         = {2026},
  howpublished = {\url{https://github.com/betterweizh/lecture-notes-operation-research}},
  note         = {Lecture notes}
}
```

When you use a result, an example or a figure, please also cite the original
source given in the chapter's references.

## Table of contents

```{toctree}
:caption: Part I — Preliminaries
:numbered:
:maxdepth: 1

01-01-linear-algebra
01-02-matrix-calculus
01-03-introduction
```

```{toctree}
:caption: Part II — Convex Optimization
:numbered:
:maxdepth: 1

02-01-convex-sets
02-02-convex-functions
02-03-convex-optimization-duality
```

```{toctree}
:caption: Part III — Linear Programming
:numbered:
:maxdepth: 1

03-01-lp-model-building
03-02-lp-geometry
03-03-simplex-method
03-04-two-phase-big-m
03-05-lp-duality
03-06-revised-dual-simplex
03-07-sensitivity-analysis
```

```{toctree}
:caption: Part IV — Network Optimization
:numbered:
:maxdepth: 1

04-01-networks-mst
04-02-shortest-paths
04-03-maximum-flow
04-04-minimum-cost-flow
```

```{toctree}
:caption: Part V — Dynamic Programming
:numbered:
:maxdepth: 1

05-01-dp-fundamentals
05-02-dp-equipment-replacement
05-03-dp-resource-allocation
05-04-dp-shortest-paths
```

```{toctree}
:caption: Part VI — Nonlinear Programming
:numbered:
:maxdepth: 1

06-01-nlp-models
06-02-optimality-conditions
06-03-kkt-conditions
06-04-line-search
06-05-unconstrained-algorithms
06-06-constrained-algorithms
```

```{toctree}
:caption: Part VII — Linearization
:numbered:
:maxdepth: 1

07-01-logical-constraints
07-02-linearization-techniques
07-03-piecewise-linear-approximation
```

```{toctree}
:caption: Part VIII — Examples and Problems
:numbered:
:maxdepth: 1

08-01-examples-lp-modeling
08-02-examples-simplex
08-03-examples-duality
08-04-examples-revised-sensitivity
08-05-examples-networks
08-06-examples-dp-fundamentals
08-07-examples-dp-applications
08-08-examples-nlp-optimality
08-09-examples-nlp-algorithms
```

```{toctree}
:caption: Part IX — Implementation
:numbered:
:maxdepth: 1

09-01-code-setup-solvers
09-02-code-lp-modeling
09-03-code-lp-algorithms
09-04-code-networks
09-05-code-dynamic-programming
09-06-code-nlp
```

```{toctree}
:caption: Appendix
:maxdepth: 1

notation
glossary
references
```
