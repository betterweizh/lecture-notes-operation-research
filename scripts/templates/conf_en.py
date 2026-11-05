"""Sphinx configuration of the English edition.

All settings are shared in ``docs/_shared/conf_common.py``. A copy of this
file is kept in ``scripts/templates/conf_en.py``; ``scripts/build_docs.py``
recreates ``docs/en/conf.py`` from it when it is missing.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_shared"))

from conf_common import *  # noqa: E402,F403

project = "Lecture Notes on Operations Research"
html_title = project

apply_language("en", globals())  # noqa: F405
