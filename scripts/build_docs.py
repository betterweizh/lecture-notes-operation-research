#!/usr/bin/env python3
"""Build the book website.

Usage::

    python scripts/build_docs.py                  # every edition -> _build/html/
    python scripts/build_docs.py --strict         # warnings are errors (-W --keep-going)
    python scripts/build_docs.py --lang en        # one edition only
    python scripts/build_docs.py --outdir /tmp/build-me   # private output folder
    python scripts/build_docs.py --figures-only   # only assemble the figures
    python scripts/build_docs.py --clean          # rebuild the editions from scratch
    python scripts/build_docs.py --serve          # build, then serve on http://127.0.0.1:8000/
    python scripts/build_docs.py --prune-legacy   # delete renamed chapter files

Steps:

1. Recreate missing generated files from ``scripts/templates/``: the
   ``conf.py`` of an edition that has an ``index.md``, and the GitHub
   workflow.
2. Assemble the figures of every edition that has a ``conf.py``: copy every
   image in ``docs/_shared/figures/`` into ``docs/<lang>/figures/``, then
   find every ``figures/<name>`` reference in the edition's Markdown files
   and run ``make_figures.py`` for the referenced figures that are missing
   or older than their generator (``scripts/figures/pPPcCC.py`` or
   ``common.py``). A reference with neither a shared file nor a generator is
   an error naming the chapter, line and figure.
3. Build each edition with ``python -m sphinx -b html docs/<lang>
   <outdir>/<lang>``.
4. Write the root page ``<outdir>/index.html`` (redirect to the reader's
   edition, or straight to ``en/`` with a single edition) and
   ``.nojekyll``; report legacy chapter files.

Figures need numpy and matplotlib. If this interpreter lacks them, the
figures are drawn with ``$FIGURES_PYTHON``, or the first ``python3`` /
``python`` on PATH that has them.
"""

from __future__ import annotations

import argparse
import functools
import html
import http.server
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
SHARED = DOCS / "_shared"
SHARED_FIGURES = SHARED / "figures"
SCRIPTS = ROOT / "scripts"
TEMPLATES = SCRIPTS / "templates"
MAKE_FIGURES = SCRIPTS / "make_figures.py"
FIGURE_COMMON = SCRIPTS / "figures" / "common.py"
WORKFLOW = ROOT / ".github" / "workflows" / "docs.yml"
DEFAULT_OUTDIR = ROOT / "_build" / "html"

IMAGE_SUFFIXES = (".svg", ".png", ".jpg", ".jpeg", ".gif", ".webp", ".pdf")
FIGURE_REF = re.compile(
    r"(?<![\w./-])figures/([A-Za-z0-9_.\-]+\.(?:svg|png|jpe?g|gif|webp|pdf))\b",
    re.IGNORECASE,
)
FENCE = re.compile(r"^\s*(`{3,}|~{3,})\s*(.*)$")

sys.path.insert(0, str(SHARED))
import conf_common  # noqa: E402  (LANGUAGES, LANG_STORAGE_KEY, LEGACY_CHAPTER_FILES)

LANGUAGES: dict = conf_common.LANGUAGES
STORAGE_KEY: str = conf_common.LANG_STORAGE_KEY
LEGACY_CHAPTER_FILES: list = list(conf_common.LEGACY_CHAPTER_FILES)


def fail(message: str) -> None:
    print(f"error: {message}", file=sys.stderr)
    sys.exit(1)


def rel(path: Path) -> str:
    """Path relative to the repository root when possible."""
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


# ---------------------------------------------------------------------------
# Step 1: templates
# ---------------------------------------------------------------------------


def recreate_missing_files() -> None:
    """Recreate conf.py files and the workflow from scripts/templates/."""
    for code in LANGUAGES:
        edition = DOCS / code
        template = TEMPLATES / f"conf_{code}.py"
        if (edition / "index.md").exists() and not (
            edition / "conf.py"
        ).exists():
            if template.exists():
                shutil.copy2(template, edition / "conf.py")
                print(
                    f"recreated {rel(edition / 'conf.py')} from {rel(template)}"
                )
    template = TEMPLATES / "github-workflow-docs.yml"
    if not WORKFLOW.exists() and template.exists():
        WORKFLOW.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(template, WORKFLOW)
        print(f"recreated {rel(WORKFLOW)} from {rel(template)}")


def editions() -> list[str]:
    """Edition codes that have a conf.py, in LANGUAGES order."""
    codes = [
        d.name
        for d in DOCS.iterdir()
        if d.is_dir() and (d / "conf.py").exists()
    ]
    codes = [c for c in codes if not c.startswith("_")]
    order = list(LANGUAGES)
    return sorted(
        codes, key=lambda c: (order.index(c) if c in order else len(order), c)
    )


# ---------------------------------------------------------------------------
# Step 2: figures
# ---------------------------------------------------------------------------


@functools.lru_cache(maxsize=1)
def figure_python() -> str:
    """An interpreter that can import numpy and matplotlib."""
    candidates = [
        os.environ.get("FIGURES_PYTHON"),
        sys.executable,
        shutil.which("python3"),
        shutil.which("python"),
    ]
    tried = []
    for exe in candidates:
        if not exe or exe in tried:
            continue
        tried.append(exe)
        try:
            ok = (
                subprocess.run(
                    [exe, "-c", "import numpy, matplotlib"],
                    capture_output=True,
                    timeout=120,
                    check=False,
                ).returncode
                == 0
            )
        except (OSError, subprocess.SubprocessError):
            ok = False
        if ok:
            return exe
    fail(
        "drawing figures needs numpy and matplotlib, and none of these interpreters has "
        f"them: {', '.join(tried)}. Install requirements.txt or set FIGURES_PYTHON."
    )
    return ""  # unreachable


@functools.lru_cache(maxsize=1)
def figure_registry() -> dict:
    """``{"figures": {name: module path}, "broken": {module: error}}``."""
    if not any((SCRIPTS / "figures").glob("p[0-9][0-9]c[0-9][0-9].py")):
        return {"figures": {}, "broken": {}}
    result = subprocess.run(
        [figure_python(), str(MAKE_FIGURES), "--list", "--json"],
        capture_output=True,
        text=True,
        check=False,
    )
    sys.stderr.write(result.stderr)
    if result.returncode != 0:
        fail("make_figures.py --list failed (see above)")
    return json.loads(result.stdout)


def legacy_names() -> set[str]:
    return {Path(name).name for name in LEGACY_CHAPTER_FILES}


def source_files(lang: str) -> list[Path]:
    """Markdown files of an edition that Sphinx builds."""
    edition = DOCS / lang
    files = []
    for path in sorted(edition.rglob("*.md")):
        parts = path.relative_to(edition).parts
        if parts[0] in ("figures", "_build") or path.name == "README.md":
            continue
        if " " in path.name or path.name in legacy_names():
            continue
        files.append(path)
    return files


def figure_references(lang: str) -> dict[str, list[tuple[Path, int]]]:
    """``{file name: [(markdown file, line number), ...]}`` outside code blocks."""
    refs: dict[str, list[tuple[Path, int]]] = {}
    for path in source_files(lang):
        stack: list[tuple[str, int, bool]] = (
            []
        )  # (fence char, length, is code)
        for number, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(), 1
        ):
            match = FENCE.match(line)
            if match:
                fence, info = match.group(1), match.group(2).strip()
                top = stack[-1] if stack else None
                closes = (
                    top
                    and not info
                    and fence[0] == top[0]
                    and len(fence) >= top[1]
                )
                if closes:
                    stack.pop()
                    continue
                if not (top and top[2]):  # fences inside a code block are text
                    stack.append(
                        (fence[0], len(fence), not info.startswith("{"))
                    )
                    if info.startswith("{"):
                        pass  # directive line, e.g. ```{figure} figures/x.svg: scan it
                    else:
                        continue
            if stack and stack[-1][2]:
                continue  # inside a code block
            for ref in FIGURE_REF.findall(line):
                refs.setdefault(ref, []).append((path, number))
    return refs


def newer_than(target: Path, *sources: Path) -> bool:
    """True if ``target`` exists and is at least as new as every source."""
    if not target.exists():
        return False
    stamp = target.stat().st_mtime
    return all(not s.exists() or s.stat().st_mtime <= stamp for s in sources)


def assemble_figures(langs: list[str]) -> None:
    """Copy shared figures and draw the referenced generated ones."""
    shared = {}
    if SHARED_FIGURES.is_dir():
        shared = {
            p.name: p
            for p in sorted(SHARED_FIGURES.iterdir())
            if p.is_file()
            and p.suffix.lower() in IMAGE_SUFFIXES
            and not p.name.startswith(".")
        }

    errors: list[str] = []
    for lang in langs:
        figures_dir = DOCS / lang / "figures"
        figures_dir.mkdir(exist_ok=True)

        # Existing figures, copied when missing or changed.
        for name, src in shared.items():
            dst = figures_dir / name
            if (
                not newer_than(dst, src)
                or dst.stat().st_size != src.stat().st_size
            ):
                shutil.copy2(src, dst)

        # Generated figures that the chapters reference.
        to_draw: dict[str, str] = {}
        references = figure_references(lang)
        for file_name, places in sorted(references.items()):
            if file_name in shared:
                continue
            stem = file_name.rsplit(".", 1)[0]
            generators = figure_registry()["figures"]
            if stem not in generators:
                for path, line in places:
                    errors.append(
                        f"{rel(path)}:{line}: figures/{file_name} is neither in "
                        "docs/_shared/figures/ nor drawn by a generator in scripts/figures/pPPcCC.py"
                    )
                continue
            module = ROOT / generators[stem]
            if not newer_than(figures_dir / file_name, module, FIGURE_COMMON):
                to_draw[stem] = file_name

        if to_draw:
            print(
                f"drawing {len(to_draw)} figure(s) for {lang}: {', '.join(sorted(to_draw))}"
            )
            command = [
                figure_python(),
                str(MAKE_FIGURES),
                "--lang",
                lang,
                *sorted(to_draw),
            ]
            if subprocess.run(command, check=False).returncode != 0:
                errors.append(f"make_figures.py failed for {lang} (see above)")
            for stem, file_name in sorted(to_draw.items()):
                if not (figures_dir / file_name).exists():
                    path, line = references[file_name][0]
                    errors.append(
                        f"{rel(path)}:{line}: generator {stem!r} did not write "
                        f"figures/{file_name} (check the file extension)"
                    )

    if errors:
        broken = (
            figure_registry()["broken"]
            if any("neither" in e for e in errors)
            else {}
        )
        for message in errors:
            print(f"error: {message}", file=sys.stderr)
        for module, message in broken.items():
            print(
                f"note: {module} could not be loaded: {message}",
                file=sys.stderr,
            )
        sys.exit(1)

    # One name, one figure: warn when a shared file is also generated.
    if shared and any((SCRIPTS / "figures").glob("p[0-9][0-9]c[0-9][0-9].py")):
        generators = figure_registry()["figures"]
        for name in shared:
            stem = name.rsplit(".", 1)[0]
            if stem in generators:
                print(
                    f"warning: figure {stem!r} exists in docs/_shared/figures/ and is also "
                    f"generated by {generators[stem]}; the shared file is used",
                    file=sys.stderr,
                )


# ---------------------------------------------------------------------------
# Steps 3 and 4: Sphinx, root page
# ---------------------------------------------------------------------------


def build_edition(
    lang: str, outdir: Path, *, strict: bool, clean: bool, quiet: bool
) -> int:
    target = outdir / lang
    source = (DOCS / lang).resolve()
    if target == source or target in source.parents:
        fail(f"--outdir {outdir} would write into the sources docs/{lang}/")
    if clean and target.exists():
        # Delete only a Sphinx output folder, never a folder of the user's
        # that --outdir happens to name.
        if any(target.iterdir()) and not (
            (target / ".buildinfo").exists() or (target / ".doctrees").exists()
        ):
            fail(
                f"--clean: {target} is not a Sphinx output folder "
                "(no .buildinfo); not deleting it"
            )
        shutil.rmtree(target)
    command = [
        sys.executable,
        "-m",
        "sphinx",
        "-b",
        "html",
        str(DOCS / lang),
        str(target),
    ]
    if strict:
        command += ["-W", "--keep-going"]
    if quiet:
        command.append("-q")
    print(f"building {lang}: {' '.join(command[1:])}")
    return subprocess.run(command, check=False).returncode


def edition_title(lang: str) -> str:
    """The ``project`` string of an edition's conf.py (or template)."""
    for conf in (DOCS / lang / "conf.py", TEMPLATES / f"conf_{lang}.py"):
        if conf.exists():
            match = re.search(
                r'^project\s*=\s*["\'](.+?)["\']',
                conf.read_text("utf-8"),
                re.M,
            )
            if match:
                return match.group(1)
    return "Lecture Notes"


def hreflang(code: str) -> str:
    return LANGUAGES.get(code, {}).get("sphinx", code).replace("_", "-")


def root_page(built: list[str]) -> str:
    """HTML of the root index.html for the editions in ``built``."""
    default = "en" if "en" in built else built[0]
    titles = [edition_title(code) for code in built]
    title = html.escape(" · ".join(titles))
    if len(built) == 1:
        target = f"{default}/index.html"
        return f"""<!DOCTYPE html>
<html lang="{hreflang(default)}">
<head>
<meta charset="utf-8">
<title>{title}</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta http-equiv="refresh" content="0; url={target}">
<link rel="canonical" href="{target}">
<script>location.replace("{target}" + location.hash);</script>
</head>
<body>
<p><a href="{target}">{title}</a></p>
</body>
</html>
"""
    alternates = "\n".join(
        f'<link rel="alternate" hreflang="{hreflang(c)}" href="{c}/index.html">'
        for c in built
    )
    links = " | ".join(
        f'<a href="{c}/index.html" hreflang="{hreflang(c)}">'
        f'{html.escape(LANGUAGES.get(c, {}).get("label", c))}</a>'
        for c in built
    )
    codes = json.dumps(built)
    return f"""<!DOCTYPE html>
<html lang="{hreflang(default)}">
<head>
<meta charset="utf-8">
<title>{title}</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
{alternates}
<link rel="alternate" hreflang="x-default" href="{default}/index.html">
<noscript><meta http-equiv="refresh" content="0; url={default}/index.html"></noscript>
<script>
(function () {{
  var KEY = {json.dumps(STORAGE_KEY)}, editions = {codes}, lang = null;
  try {{ lang = window.localStorage.getItem(KEY); }} catch (e) {{ lang = null; }}
  if (editions.indexOf(lang) < 0) {{
    var nav = (navigator.languages && navigator.languages[0]) || navigator.language || "";
    lang = /^zh/i.test(nav) && editions.indexOf("zh") >= 0 ? "zh" : {json.dumps(default)};
  }}
  location.replace(lang + "/index.html");
}})();
</script>
<style>
  :root {{ color-scheme: light dark; }}
  body {{ background: #ffffff; color: #2b2b2b; font-family: system-ui, sans-serif;
         margin: 0; display: grid; min-height: 100vh; place-items: center; }}
  a {{ color: #2a6db0; }}
  @media (prefers-color-scheme: dark) {{
    body {{ background: #14181e; color: #e3e3e3; }}
    a {{ color: #79b0ef; }}
  }}
</style>
</head>
<body>
<main>
<h1>{title}</h1>
<p>{links}</p>
</main>
</body>
</html>
"""


def write_root(outdir: Path) -> None:
    built = [
        code for code in editions() if (outdir / code / "index.html").exists()
    ]
    if not built:
        return
    (outdir / "index.html").write_text(root_page(built), encoding="utf-8")
    (outdir / ".nojekyll").write_text("", encoding="utf-8")
    print(f"wrote {rel(outdir / 'index.html')} ({', '.join(built)})")


# Figure modules named before the nine-part restructure (ch01.py to
# ch39.py); they are no longer loaded (only pPPcCC.py is) and
# --prune-legacy deletes them. An explicit list, so that no other file of
# the user's is ever deleted.
LEGACY_FIGURE_MODULES = [f"ch{n:02d}.py" for n in range(1, 40)]


def report_legacy(prune: bool) -> None:
    legacy_modules = [SCRIPTS / "figures" / n for n in LEGACY_FIGURE_MODULES]
    for path in [p for p in legacy_modules if p.is_file()]:
        if prune:
            path.unlink()
            print(f"deleted legacy file {rel(path)}")
        else:
            print(
                f"note: legacy file {rel(path)} is not used; delete it with --prune-legacy"
            )
    for code in LANGUAGES:
        for name in LEGACY_CHAPTER_FILES:
            path = DOCS / code / name
            if not path.exists():
                continue
            if prune:
                path.unlink()
                print(f"deleted legacy file {rel(path)}")
            else:
                print(
                    f"note: legacy file {rel(path)} is excluded from the build; "
                    "delete it with --prune-legacy"
                )


class _Server(http.server.ThreadingHTTPServer):
    # On Windows SO_REUSEADDR lets a second server bind a port that is
    # already in use, and the browser may then reach the old server.
    allow_reuse_address = os.name != "nt"


def open_server(outdir: Path, port: int) -> _Server:
    """Bind the port (before building, so a busy port fails at once)."""
    handler = functools.partial(
        http.server.SimpleHTTPRequestHandler, directory=str(outdir)
    )
    try:
        return _Server(("127.0.0.1", port), handler)
    except OSError as exc:
        fail(
            f"cannot serve on port {port} ({exc.strerror or exc}); "
            "is another server running? Choose another one with --port"
        )
        raise  # unreachable


def serve(server: _Server, outdir: Path, port: int) -> None:
    with server:
        print(
            f"serving {rel(outdir)} on http://127.0.0.1:{port}/ (Ctrl+C to stop)"
        )
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build the book website.")
    parser.add_argument(
        "--lang", choices=list(LANGUAGES), help="build only this edition"
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="warnings are errors (-W --keep-going)",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="delete the edition's output first",
    )
    parser.add_argument(
        "--serve", action="store_true", help="serve the site after building"
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8000,
        help="port for --serve (default 8000)",
    )
    parser.add_argument(
        "--prune-legacy",
        action="store_true",
        help="delete LEGACY_CHAPTER_FILES and the chNN.py figure modules",
    )
    parser.add_argument(
        "--figures-only", action="store_true", help="only assemble the figures"
    )
    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTDIR),
        help="output folder (default _build/html); editions go to <outdir>/<lang>/",
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true", help="pass -q to Sphinx"
    )
    args = parser.parse_args(argv)

    recreate_missing_files()
    available = editions()
    if not available:
        fail("no edition found (docs/<lang>/conf.py)")
    if args.lang and args.lang not in available:
        fail(f"edition {args.lang!r} has no docs/{args.lang}/conf.py")
    langs = [args.lang] if args.lang else available

    if args.prune_legacy:
        report_legacy(prune=True)

    assemble_figures(langs)
    if args.figures_only:
        print("figures assembled for " + ", ".join(langs))
        return 0

    outdir = Path(args.outdir).resolve()
    server = open_server(outdir, args.port) if args.serve else None
    status = 0
    for lang in langs:
        code = build_edition(
            lang,
            outdir,
            strict=args.strict,
            clean=args.clean,
            quiet=args.quiet,
        )
        if code:
            print(
                f"error: Sphinx failed for {lang} (exit code {code})",
                file=sys.stderr,
            )
            status = code
    write_root(outdir)
    report_legacy(prune=False)
    if status == 0:
        for lang in langs:
            print(f"built {lang}: {rel(outdir / lang / 'index.html')}")
    if server is not None:
        if status == 0:
            serve(server, outdir, args.port)
        else:
            server.server_close()
    return status


if __name__ == "__main__":
    sys.exit(main())
