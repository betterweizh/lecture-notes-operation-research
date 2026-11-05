#!/usr/bin/env python3
"""Render figures to PNG so they can be looked at with an image viewer.

Usage (needs Python Playwright with Chromium)::

    python3 scripts/figures/preview.py 03-02-product-mix-region 04-02-dijkstra-network
    python3 scripts/figures/preview.py --make 03-02-product-mix-region   # draw first
    python3 scripts/figures/preview.py --width 520 --dark 03-02-product-mix-region
    python3 scripts/figures/preview.py path/to/any-figure.svg

Each argument is a figure name, looked up in ``docs/<lang>/figures/`` and
then ``docs/_shared/figures/`` (any extension), or a path to an image file.
The PNG goes to ``/tmp/fig-preview/<name>.png``: the figure on a white page
with a 16 px margin, at 2x resolution. With ``--dark`` a second PNG,
``<name>.dark.png``, shows it as the website shows line art in dark mode
(``invert(0.88) hue-rotate(180deg)`` on a dark background).

``--width`` displays the figure at that CSS width, like ``:width:`` in the
chapter. ``--make`` first runs ``make_figures.py`` for the given names.
An SVG that is not well-formed XML is reported with line and column.
"""

from __future__ import annotations

import argparse
import html
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
SHARED_FIGURES = DOCS / "_shared" / "figures"
DEFAULT_OUT = Path("/tmp/fig-preview")
IMAGE_SUFFIXES = (".svg", ".png", ".jpg", ".jpeg", ".gif", ".webp")

LIGHT_BG = "#ffffff"
DARK_BG = "#14181e"  # close to the book theme's dark background


def find_figure(arg: str, lang: str) -> Path | None:
    """Resolve a figure name or a path to an existing image file."""
    path = Path(arg)
    if path.suffix.lower() in IMAGE_SUFFIXES and path.exists():
        return path.resolve()
    for folder in (DOCS / lang / "figures", SHARED_FIGURES):
        if (folder / arg).is_file():
            return (folder / arg).resolve()
        for suffix in IMAGE_SUFFIXES:
            candidate = folder / f"{arg}{suffix}"
            if candidate.exists():
                return candidate.resolve()
    return None


def check_svg(path: Path) -> str | None:
    """Return an error message if the SVG is not well-formed XML."""
    try:
        ET.parse(path)
    except ET.ParseError as exc:
        line, column = exc.position
        return f"{path}: not well-formed SVG (line {line}, column {column}): {exc}"
    return None


def wrapper_html(image: Path, width: int | None, dark: bool) -> str:
    """HTML page showing one image, as the book shows it."""
    background = DARK_BG if dark else LIGHT_BG
    size = f"width: {width}px; height: auto;" if width else ""
    image_filter = "filter: invert(0.88) hue-rotate(180deg);" if dark else ""
    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><style>
  html, body {{ margin: 0; background: {background}; }}
  #frame {{ display: inline-block; padding: 16px; background: {background}; }}
  img {{ display: block; max-width: none; {size} {image_filter} }}
</style></head>
<body><div id="frame"><img src="{html.escape(image.as_uri())}" alt=""></div></body></html>
"""


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Render figures to PNG with Playwright."
    )
    parser.add_argument(
        "figures", nargs="+", help="figure names or image paths"
    )
    parser.add_argument(
        "--lang", default="en", help="edition whose figures/ folder to use"
    )
    parser.add_argument(
        "--out",
        default=str(DEFAULT_OUT),
        help="output folder (default /tmp/fig-preview)",
    )
    parser.add_argument(
        "--width",
        type=int,
        help="display width in CSS px (default: natural size)",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=2.0,
        help="device scale factor (default 2)",
    )
    parser.add_argument(
        "--dark", action="store_true", help="also write the dark-mode view"
    )
    parser.add_argument(
        "--make",
        action="store_true",
        help="run make_figures.py for the names first",
    )
    args = parser.parse_args(argv)

    errors = 0
    if args.make:
        # make_figures.py draws every known name even if others fail.
        names = [
            a
            for a in args.figures
            if Path(a).suffix.lower() not in IMAGE_SUFFIXES
        ]
        if names:
            command = [
                sys.executable,
                str(ROOT / "scripts" / "make_figures.py"),
                "--lang",
                args.lang,
                *names,
            ]
            if subprocess.run(command, check=False).returncode != 0:
                errors += 1
                print(
                    "warning: drawing failed for some figures (see above); "
                    "a preview may show an older file",
                    file=sys.stderr,
                )

    # Absolute: the wrapper pages are opened by file:// URI.
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)

    jobs: list[tuple[str, Path]] = []
    for arg in args.figures:
        path = find_figure(arg, args.lang)
        if path is None:
            print(
                f"error: {arg}: not found in docs/{args.lang}/figures/ or docs/_shared/figures/ "
                "(draw it first: make_figures.py NAME, or use --make)",
                file=sys.stderr,
            )
            errors += 1
            continue
        if path.suffix.lower() == ".svg":
            message = check_svg(path)
            if message:
                print(f"error: {message}", file=sys.stderr)
                errors += 1
                continue
        jobs.append((path.name.rsplit(".", 1)[0], path))

    if jobs:
        try:
            from playwright.sync_api import sync_playwright
        except ImportError:
            print(
                "error: Python Playwright is not installed (pip install playwright)",
                file=sys.stderr,
            )
            return 1
        with sync_playwright() as pw:
            browser = pw.chromium.launch()
            page = browser.new_page(
                viewport={"width": 1800, "height": 1400},
                device_scale_factor=args.scale,
            )
            for stem, path in jobs:
                variants = [("", False)] + (
                    [(".dark", True)] if args.dark else []
                )
                for suffix, dark in variants:
                    page_file = out / f".{stem}{suffix}.html"
                    page_file.write_text(
                        wrapper_html(path, args.width, dark), encoding="utf-8"
                    )
                    try:
                        page.goto(page_file.as_uri())
                        page.wait_for_function(
                            "document.querySelector('img').complete"
                        )
                        natural = page.evaluate(
                            "document.querySelector('img').naturalWidth"
                        )
                        if not natural:
                            print(
                                f"error: {path}: the browser could not render it",
                                file=sys.stderr,
                            )
                            errors += 1
                            continue
                        png = out / f"{stem}{suffix}.png"
                        page.locator("#frame").screenshot(path=str(png))
                        box = page.locator("img").bounding_box() or {}
                        print(
                            f"{png}  ({box.get('width', 0):.0f} x {box.get('height', 0):.0f} CSS px)"
                        )
                    finally:
                        page_file.unlink(missing_ok=True)
            browser.close()
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
