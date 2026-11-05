"""Sphinx configuration of the Simplified Chinese edition.

All settings are shared in ``docs/_shared/conf_common.py``. A copy of this
file is kept in ``scripts/templates/conf_zh.py``; ``scripts/build_docs.py``
recreates ``docs/zh/conf.py`` from it when ``docs/zh/index.md`` exists but
``conf.py`` is missing.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_shared"))

from conf_common import *  # noqa: E402,F403

project = "运筹学学习笔记"
html_title = project

apply_language("zh", globals())  # noqa: F405
