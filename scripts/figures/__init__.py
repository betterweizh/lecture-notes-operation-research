"""Generators of the book's figures.

One module per chapter, ``pPPcCC.py``, each with a ``FIGURES`` dict that maps
figure names (``"NN-descriptive-name"``) to functions ``func(outdir,
lang="en")``. ``scripts/make_figures.py`` imports every ``pPPcCC.py`` module
and draws the figures; shared style and helpers are in ``common.py``.
See ``README.md`` in this folder.
"""
