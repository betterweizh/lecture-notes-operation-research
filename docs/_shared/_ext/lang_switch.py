"""Language switch between the editions of the book.

For every HTML page this extension computes the URL of the same page in
each other edition and exposes three template variables, used by
``_templates/lang-switch.html``:

``lang_switch_current``
    code of this edition (``"en"``).
``lang_switch_current_label``
    its label (``"English"``).
``lang_switch_alternates``
    list of ``{"code", "label", "short", "url"}`` for the other editions.

Only editions that exist (``docs/<code>/conf.py``) are offered. If the page
has no counterpart in the other edition, the link goes to that edition's
front page. ``genindex`` maps to ``index``.

Configuration (set by ``apply_language()`` in ``conf_common.py``):
``lang_switch_current`` and ``lang_switch_languages`` (``LANGUAGES``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from sphinx.application import Sphinx


def _counterpart(docs_dir: Path, code: str, page: str) -> str:
    """Return the page name to link to in edition ``code``."""
    edition = docs_dir / code
    if any(
        (edition / f"{page}{suffix}").exists() for suffix in (".md", ".rst")
    ):
        return page
    return "index"


def add_alternates(
    app: Sphinx, pagename: str, templatename: str, context: dict, doctree: Any
) -> None:
    current = app.config.lang_switch_current
    languages = app.config.lang_switch_languages or {}
    docs_dir = Path(app.confdir).resolve().parent  # docs/
    page = "index" if pagename == "genindex" else pagename
    up = "../" * (pagename.count("/") + 1)

    alternates = []
    for code, info in languages.items():
        if code == current or not (docs_dir / code / "conf.py").exists():
            continue
        # Pages built from a source file must exist in the other edition;
        # generated pages (search) exist in every edition.
        target = (
            _counterpart(docs_dir, code, page) if doctree is not None else page
        )
        alternates.append(
            {
                "code": code,
                "label": info.get("label", code),
                "short": info.get("short", code),
                "url": f"{up}{code}/{target}.html",
            }
        )

    context["lang_switch_current"] = current
    context["lang_switch_current_label"] = languages.get(current, {}).get(
        "label", current
    )
    context["lang_switch_alternates"] = alternates


def setup(app: Sphinx) -> dict[str, Any]:
    app.add_config_value("lang_switch_current", "en", "html", types=[str])
    app.add_config_value("lang_switch_languages", {}, "html", types=[dict])
    app.connect("html-page-context", add_alternates, priority=100)
    return {
        "version": "1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
