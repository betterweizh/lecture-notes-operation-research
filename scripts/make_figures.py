#!/usr/bin/env python3
"""Draw the book's generated figures.

Every module ``scripts/figures/pPPcCC.py`` has a ``FIGURES`` dict that maps a
figure name (``"PP-CC-descriptive-name"``) to a function ``func(outdir,
lang="en")``. This script imports all of them, merges the dicts (a name
defined twice is an error) and draws all figures, or only the names given,
into ``docs/<lang>/figures/`` for every edition that has a ``conf.py``.

Usage (needs numpy and matplotlib)::

    python3 scripts/make_figures.py                      # all figures
    python3 scripts/make_figures.py 03-02-half-space     # only these
    python3 scripts/make_figures.py --list               # names and modules
    python3 scripts/make_figures.py --lang en NAME ...   # one edition
    python3 scripts/make_figures.py --outdir DIR NAME    # write into DIR
    python3 scripts/make_figures.py --example --outdir /tmp/fig-preview/example

``scripts/build_docs.py`` runs this script for the figures that chapters
reference and that are missing or older than their generator.

A module that fails to import is reported and skipped, so one broken
module does not stop the figures of other chapters.
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys
import traceback
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parent
DOCS = ROOT / "docs"
FIGURE_PACKAGE = SCRIPTS / "figures"
EXAMPLE_MODULE = "_example"

if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


def editions() -> list[str]:
    """Edition codes that have a conf.py, English first."""
    codes = [
        d.name
        for d in DOCS.iterdir()
        if d.is_dir()
        and not d.name.startswith("_")
        and (d / "conf.py").exists()
    ]
    return sorted(codes, key=lambda c: (c != "en", c))


def module_names(example: bool = False) -> list[str]:
    """Figure modules to load: every pPPcCC.py, or only the example module."""
    if example:
        return [EXAMPLE_MODULE]
    return sorted(
        p.stem for p in FIGURE_PACKAGE.glob("p[0-9][0-9]c[0-9][0-9].py")
    )


def load_registry(example: bool = False) -> tuple[dict, dict]:
    """Import the figure modules and merge their FIGURES dicts.

    Returns
    -------
    registry : dict
        ``{name: (function, "scripts/figures/pPPcCC.py")}``.
    broken : dict
        ``{"scripts/figures/pPPcCC.py": error message}`` for modules that
        could not be imported or have no valid FIGURES dict.
    """
    registry: dict[str, tuple] = {}
    broken: dict[str, str] = {}
    duplicates: list[str] = []
    for name in module_names(example):
        rel = f"scripts/figures/{name}.py"
        try:
            module = importlib.import_module(f"figures.{name}")
        except Exception as exc:  # noqa: BLE001 - report any import error
            broken[rel] = "".join(
                traceback.format_exception_only(type(exc), exc)
            ).strip()
            continue
        figures = getattr(module, "FIGURES", None)
        if not isinstance(figures, dict):
            broken[rel] = "no FIGURES dict"
            continue
        for fig_name, func in figures.items():
            if fig_name in registry:
                duplicates.append(
                    f"figure name {fig_name!r} is defined in {registry[fig_name][1]} and in {rel}"
                )
                continue
            if not callable(func):
                broken[rel] = f"FIGURES[{fig_name!r}] is not a function"
                continue
            registry[fig_name] = (func, rel)
            # Names start with the module's part and chapter: p03c02 ->
            # "03-02-" (the example module has no number).
            prefix = f"{name[1:3]}-{name[4:6]}-"
            if not example and not fig_name.startswith(prefix):
                print(
                    f"warning: {rel}: figure {fig_name!r} does not start with {prefix}",
                    file=sys.stderr,
                )
    if duplicates:
        for message in duplicates:
            print(f"error: {message}", file=sys.stderr)
        sys.exit(2)
    return registry, broken


def draw(
    registry: dict, names: list[str], outdir: Path, lang: str
) -> list[str]:
    """Draw ``names`` into ``outdir``; return the names that failed."""
    outdir.mkdir(parents=True, exist_ok=True)
    failed = []
    for name in names:
        func, rel = registry[name]
        # Files of this figure before drawing: a generator that writes
        # another name must not pass because an old file is still there
        # (the returned path also counts: mtimes can be coarse).
        before = {
            p.name: p.stat().st_mtime_ns for p in outdir.glob(f"{name}.*")
        }
        try:
            if lang == "zh":
                # Chinese labels need a CJK font first (see common.ZH_RC).
                import matplotlib as mpl

                from figures.common import ZH_RC

                with mpl.rc_context(ZH_RC):
                    result = func(outdir, lang=lang)
            else:
                result = func(outdir, lang=lang)
        except Exception:  # noqa: BLE001 - report and continue
            print(f"error: {name} ({rel}) failed:", file=sys.stderr)
            traceback.print_exc()
            failed.append(name)
            continue
        written = sorted(
            p.name
            for p in outdir.glob(f"{name}.*")
            if not p.name.endswith(".tmp")
            and (
                before.get(p.name) != p.stat().st_mtime_ns
                or (result and Path(result).resolve() == p.resolve())
            )
        )
        if not written:
            print(
                f"error: {name} ({rel}) did not write {name}.svg",
                file=sys.stderr,
            )
            failed.append(name)
        else:
            print(
                f"  {lang}: {outdir.relative_to(ROOT) if outdir.is_relative_to(ROOT) else outdir}/{written[0]}"
            )
    return failed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Draw the generated figures into docs/<lang>/figures/."
    )
    parser.add_argument("names", nargs="*", help="figure names (default: all)")
    parser.add_argument(
        "--lang", help="draw only for this edition (default: every edition)"
    )
    parser.add_argument(
        "--outdir",
        help="write into this folder instead of docs/<lang>/figures/",
    )
    parser.add_argument(
        "--list", action="store_true", help="list figure names and modules"
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="with --list: machine-readable output",
    )
    parser.add_argument(
        "--example",
        action="store_true",
        help="use only scripts/figures/_example.py (requires --outdir unless --list)",
    )
    args = parser.parse_args(argv)
    # Paths may hold characters a Windows console pipe cannot encode.
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")
    if args.lang and args.lang not in editions():
        # Labels exist only for the editions; never create docs/<typo>/.
        parser.error(
            f"--lang {args.lang!r}: no edition docs/{args.lang}/conf.py "
            f"(editions: {', '.join(editions())})"
        )

    registry, broken = load_registry(example=args.example)
    for rel, message in broken.items():
        print(f"warning: {rel} skipped: {message}", file=sys.stderr)

    if args.list:
        if args.json:
            data = {
                "figures": {
                    n: rel for n, (_, rel) in sorted(registry.items())
                },
                "broken": broken,
            }
            print(json.dumps(data, indent=1))
        else:
            for name, (_, rel) in sorted(registry.items()):
                print(f"{name:45s} {rel}")
        return 0

    if args.example and not args.outdir:
        parser.error(
            "--example needs --outdir (example figures never go into docs/)"
        )

    # A name given twice is drawn once.
    names = list(dict.fromkeys(args.names)) or sorted(registry)
    unknown = [n for n in names if n not in registry]
    for name in unknown:
        print(f"error: no generator for figure {name!r}", file=sys.stderr)
    if unknown and broken:
        print(
            "note: these modules could not be loaded: " + ", ".join(broken),
            file=sys.stderr,
        )
    # Draw the known names anyway; the exit code still reports the unknown ones.
    names = [n for n in names if n in registry]

    if args.lang:
        langs = [args.lang]
    elif args.outdir:
        langs = ["en"]
    else:
        langs = editions()
    failed: list[str] = []
    for lang in langs:
        outdir = (
            Path(args.outdir).resolve()
            if args.outdir
            else DOCS / lang / "figures"
        )
        failed += draw(registry, names, outdir, lang)
    if failed:
        print(
            f"error: {len(failed)} figure(s) failed: {', '.join(sorted(set(failed)))}",
            file=sys.stderr,
        )
    return 1 if failed or unknown else 0


if __name__ == "__main__":
    sys.exit(main())
