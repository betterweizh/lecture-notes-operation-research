"""Number chapters by part: "<part>.<chapter>", restarting in every part.

The front page has one ``:numbered:`` toctree per part. Sphinx already
restarts chapter numbers in every numbered toctree, so the 2nd chapter of
every part would be "2" and its figures "2.1", "2.2", ... in several parts
at once. This extension prefixes every number with the part number: the
2nd chapter of Part III is "3.2", its sections "3.2.1", its figures,
tables, equations and sphinx-proof blocks "3.2.1", "3.2.2", ...

After Sphinx has assigned its numbers (``env-get-updated``, priority 900,
i.e. after the built-in toctree collector at 500), it walks the numbered
toctrees of the root document in order (part p = 1, 2, ...) and, for the
c-th chapter of part p,

* replaces the first component c of ``env.toc_secnumbers[doc]`` (section
  numbers: headings, ``{numref}`` to sections) by ``p, c``;
* does the same for the ``secnumber`` of the document title
  (previous/next links) and of the reference nodes in ``env.tocs[doc]``
  (navigation and the toctrees on the front page);
* renumbers every entry of ``env.toc_fignumbers[doc]`` (figures, tables,
  code blocks, equations with ``math_numfig``, sphinx-proof blocks) as
  ``p, c, k`` with k counting from 1 in the chapter, per type and in
  document order. (Sphinx keys its counters on the local chapter number,
  so chapters with the same position in different parts would otherwise
  share one counter.)

Unnumbered toctrees (e.g. the appendix) are skipped.
"""

from __future__ import annotations

from typing import Any

from docutils import nodes
from sphinx import addnodes
from sphinx.application import Sphinx
from sphinx.environment import BuildEnvironment


def _numbered_parts(env: BuildEnvironment, root: str) -> list[list[str]]:
    """Return the documents of each numbered toctree of ``root``, in order."""
    try:
        doctree = env.get_doctree(root)
    except FileNotFoundError:
        return []
    parts: list[list[str]] = []
    for toctree in doctree.findall(addnodes.toctree):
        if not toctree.get("numbered"):
            continue
        chapters = []
        for _title, ref in toctree["entries"]:
            # Skip external links and "self"; keep only real documents.
            if ref == "self" or "://" in ref or ref not in env.tocs:
                continue
            chapters.append(ref)
        parts.append(chapters)
    return parts


def _prefixed(number: Any, part: int, chapter: int) -> Any:
    """Return ``number`` with its first component replaced by ``part, chapter``."""
    if not number:
        return number
    return type(number)([part, chapter, *number[1:]])


def renumber_chapters(app: Sphinx, env: BuildEnvironment) -> list[str]:
    """Give chapters part-prefixed numbers; return the documents that changed."""
    changed: list[str] = []
    for part, chapters in enumerate(
        _numbered_parts(env, app.config.root_doc), start=1
    ):
        for chapter, doc in enumerate(chapters, start=1):
            secnums = env.toc_secnumbers.get(doc)
            if not secnums:
                continue
            # Already prefixed (e.g. an incremental rebuild): leave it.
            local = secnums.get("") or next(iter(secnums.values()), ())
            if len(local) >= 2 and tuple(local[:2]) == (part, chapter):
                pass
            else:
                env.toc_secnumbers[doc] = {
                    anchor: _prefixed(number, part, chapter)
                    for anchor, number in secnums.items()
                }
                title = env.titles.get(doc)
                if title is not None and title.get("secnumber"):
                    title["secnumber"] = _prefixed(
                        title["secnumber"], part, chapter
                    )
                for reference in env.tocs[doc].findall(nodes.reference):
                    if reference.get("secnumber"):
                        reference["secnumber"] = _prefixed(
                            reference["secnumber"], part, chapter
                        )
                changed.append(doc)

            # Figures, tables, equations, proof blocks: (part, chapter, k).
            for numbers in env.toc_fignumbers.get(doc, {}).values():
                ordered = sorted(
                    (tuple(number), node_id)
                    for node_id, number in numbers.items()
                    if number
                )
                for rank, (_number, node_id) in enumerate(ordered, start=1):
                    renumbered = type(numbers[node_id])([part, chapter, rank])
                    if numbers[node_id] != renumbered:
                        numbers[node_id] = renumbered
                        if doc not in changed:
                            changed.append(doc)
    return changed


def setup(app: Sphinx) -> dict[str, Any]:
    app.connect("env-get-updated", renumber_chapters, priority=900)
    return {
        "version": "1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
