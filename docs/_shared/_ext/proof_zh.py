"""Chinese names for sphinx-proof blocks (定义, 定理, 例, …).

sphinx-proof's own zh_CN catalog writes "例子 3.2.1"; Chinese textbooks
write "例 3.2.1". This extension registers the book's own catalog in
``docs/_shared/locales/zh_CN/LC_MESSAGES/proof.mo`` (source: ``proof.po``
next to it) for sphinx-proof's message domain. Load it **before**
``sphinx_proof``: Sphinx keeps the first registered catalog and uses later
ones only as fallbacks. It has no effect on other languages.

It also gives sphinx-proof's domain the ``merge_domaindata`` method that
parallel builds (``sphinx-build -j N``) need.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from sphinx.application import Sphinx

LOCALES = Path(__file__).resolve().parent.parent / "locales"


def _merge_no_domaindata(self: Any, docnames: Any, otherdata: Any) -> None:
    """sphinx-proof keeps its data in ``env.proof_list`` (merged by its own
    ``env-merge-info`` handler); its domain has nothing to merge."""


def setup(app: Sphinx) -> dict[str, Any]:
    # Same message domain as sphinx-proof ("proof"); registered first, so
    # its translations take precedence over sphinx-proof's.
    app.add_message_catalog("proof", str(LOCALES))
    # sphinx-proof (0.4) declares itself parallel-read safe, but its domain
    # lacks merge_domaindata, so "sphinx-build -j N" crashes with
    # NotImplementedError. Supply the missing (empty) merge.
    try:
        from sphinx_proof.domain import ProofDomain
    except ImportError:  # sphinx-proof missing: Sphinx reports it itself
        pass
    else:
        if "merge_domaindata" not in vars(ProofDomain):
            ProofDomain.merge_domaindata = _merge_no_domaindata
    return {"version": "1.0", "parallel_read_safe": True, "parallel_write_safe": True}
