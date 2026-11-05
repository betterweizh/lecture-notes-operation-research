#!/usr/bin/env python3
"""Run the Python code blocks of a code chapter and report failures.

Usage::

    python3 scripts/run_code_cells.py docs/en/09-02-code-lp-modeling.md
    python3 scripts/run_code_cells.py --skip-missing docs/en/09-*.md
    python3 scripts/run_code_cells.py --quiet FILE     # output only for failures
    python3 scripts/run_code_cells.py --list FILE      # show the blocks, run nothing

Every ```` ```python ```` block of a Markdown file (also indented ones, and
blocks inside directives such as ``{dropdown}``) is run in order, in one
namespace per file, like the cells of a notebook. Each block's output is
printed under a header with its line numbers; a failing block prints its
traceback with the line numbers of the Markdown file. Later blocks still
run. The exit code is 1 if any block failed.

``--skip-missing`` skips what needs a package that is not installed
instead of failing: imports of missing packages are removed from a block,
and a block that uses a name such an import would define (or a name
defined by an earlier skipped block) is skipped with a note; its imports
of installed packages still run, so later blocks can use them.

Blocks run in a temporary working directory (``--cwd`` to choose one), with
the matplotlib backend set to ``Agg`` so ``plt.show()`` does nothing.
"""

from __future__ import annotations

import argparse
import ast
import contextlib
import importlib.util
import io
import os
import re
import sys
import tempfile
import traceback
from dataclasses import dataclass, field
from pathlib import Path

FENCE = re.compile(r"^(\s*)(`{3,}|~{3,})\s*(.*)$")
PYTHON_INFO = {"python", "python3", "py"}


@dataclass
class Block:
    """One code block: its source and the Markdown lines it came from."""

    number: int  # 1-based position among the file's Python blocks
    first_line: int  # Markdown line of the first code line
    last_line: int
    source: str
    skipped_reason: str = ""
    removed_imports: list[str] = field(default_factory=list)


def extract_blocks(path: Path) -> list[Block]:
    """Return the ```python blocks of a Markdown file, in order."""
    lines = path.read_text(encoding="utf-8").splitlines()
    blocks: list[Block] = []
    # Stack of open fences: (char, length, is_code, indent, is_python, start).
    stack: list[tuple[str, int, bool, int, bool, int]] = []
    body: list[str] = []
    for number, line in enumerate(lines, 1):
        match = FENCE.match(line)
        top = stack[-1] if stack else None
        if match:
            indent, fence, info = (
                len(match.group(1)),
                match.group(2),
                match.group(3).strip(),
            )
            closes = (
                top
                and not info
                and fence[0] == top[0]
                and len(fence) >= top[1]
            )
            if closes and top is not None:
                stack.pop()
                if top[4]:  # end of a Python block
                    blocks.append(
                        Block(
                            len(blocks) + 1,
                            top[5],
                            number - 1,
                            "\n".join(body) + "\n",
                        )
                    )
                continue
            if not (
                top and top[2]
            ):  # a fence inside a code block is just text
                language = info.split()[0].lower() if info else ""
                is_code = not info.startswith("{")
                stack.append(
                    (
                        fence[0],
                        len(fence),
                        is_code,
                        indent,
                        language in PYTHON_INFO,
                        number + 1,
                    )
                )
                body = []
                continue
        if top and top[4]:
            # Remove the fence's indentation (blocks inside list items).
            strip = top[3]
            body.append(
                line[strip:] if line[:strip].strip() == "" else line.lstrip()
            )
    if stack and stack[-1][4]:
        # An unclosed fence runs to the end of the file (as in Markdown).
        start = stack[-1][5]
        print(
            f"warning: {path}:{start - 1}: the ```python block is not "
            "closed; it runs to the end of the file",
            file=sys.stderr,
        )
        blocks.append(
            Block(len(blocks) + 1, start, len(lines), "\n".join(body) + "\n")
        )
    return blocks


# ---------------------------------------------------------------------------
# --skip-missing: find imports of packages that are not installed
# ---------------------------------------------------------------------------

_available: dict[str, bool] = {}


def is_installed(module: str) -> bool:
    top = module.split(".")[0]
    if top not in _available:
        try:
            _available[top] = importlib.util.find_spec(top) is not None
        except (ImportError, ValueError):
            _available[top] = False
    return _available[top]


class _DropMissingImports(ast.NodeTransformer):
    """Remove imports of missing packages; remember the names they bind."""

    def __init__(self) -> None:
        self.unbound: dict[str, str] = {}  # name -> missing package
        self.packages: set[str] = set()

    def visit_Import(self, node: ast.Import):
        keep = []
        for alias in node.names:
            if is_installed(alias.name):
                keep.append(alias)
            else:
                package = alias.name.split(".")[0]
                self.packages.add(package)
                self.unbound[alias.asname or package] = package
        return self._replace(node, keep)

    def visit_ImportFrom(self, node: ast.ImportFrom):
        if node.level or not node.module or is_installed(node.module):
            return node
        package = node.module.split(".")[0]
        self.packages.add(package)
        for alias in node.names:
            self.unbound[alias.asname or alias.name] = package
        return self._replace(node, [])

    @staticmethod
    def _replace(node, keep):
        if len(keep) == len(node.names):
            return node
        if keep:
            node.names = keep
            return node
        return ast.copy_location(ast.Pass(), node)


def _bound_names(tree: ast.Module) -> set[str]:
    """Names a block defines at module level (not inside functions)."""
    names: set[str] = set()

    def visit(node: ast.AST) -> None:
        if isinstance(
            node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
        ):
            names.add(node.name)
            return  # do not look inside bodies
        if isinstance(node, ast.Lambda):
            return
        if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            names.add(node.id)
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            for alias in node.names:
                names.add(alias.asname or alias.name.split(".")[0])
        for child in ast.iter_child_nodes(node):
            visit(child)

    visit(tree)
    return names


def _used_names(tree: ast.AST) -> set[str]:
    return {
        n.id
        for n in ast.walk(tree)
        if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)
    }


# ---------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------


def run_file(
    path: Path, *, skip_missing: bool, quiet: bool, stop: bool, cwd: Path
) -> tuple[int, int, int]:
    """Run every block of ``path``; return (passed, failed, skipped)."""
    blocks = extract_blocks(path)
    print(f"=== {path} ({len(blocks)} Python block(s))")
    namespace: dict = {"__name__": "__main__", "__file__": str(path)}
    unavailable: set[str] = set()  # names that skipped code would have defined
    missing_packages: dict[str, set[str]] = {}
    passed = failed = skipped = 0

    for block in blocks:
        header = f"--- block {block.number} (lines {block.first_line}-{block.last_line})"
        # Pad with newlines so tracebacks show the Markdown line numbers.
        padded = "\n" * (block.first_line - 1) + block.source
        try:
            tree = ast.parse(padded, filename=str(path))
        except SyntaxError as exc:
            print(
                f"*** block {block.number} (lines {block.first_line}-{block.last_line}): FAILED (syntax error)"
            )
            print(
                "".join(
                    traceback.format_exception_only(type(exc), exc)
                ).rstrip()
            )
            failed += 1
            if stop:
                break
            continue

        if skip_missing:
            dropper = _DropMissingImports()
            tree = ast.fix_missing_locations(dropper.visit(tree))
            for name, package in dropper.unbound.items():
                missing_packages.setdefault(name, set()).add(package)
            unavailable |= set(dropper.unbound)
            needed = sorted(_used_names(tree) & unavailable)
            if needed:
                packages = sorted(
                    set().union(
                        *(missing_packages.get(n, set()) for n in needed)
                    )
                )
                reason = (
                    f"needs {', '.join(packages)} (not installed)"
                    if packages
                    else ""
                )
                reason = (
                    reason or f"uses {', '.join(needed)} from a skipped block"
                )
                print(f"{header}: skipped, {reason}")
                # Still run the block's imports of installed packages (a
                # setup block often uses the missing package too): later
                # blocks need them.
                imports = [
                    node
                    for node in tree.body
                    if isinstance(node, (ast.Import, ast.ImportFrom))
                ]
                bound = _bound_names(
                    ast.Module(
                        body=[n for n in tree.body if n not in imports],
                        type_ignores=[],
                    )
                )
                if imports:
                    module = ast.Module(body=imports, type_ignores=[])
                    try:
                        exec(compile(module, str(path), "exec"), namespace)
                    except Exception:  # noqa: BLE001 - then skip them too
                        bound |= _bound_names(module)
                unavailable |= bound
                for name in bound:
                    missing_packages.setdefault(name, set()).update(packages)
                skipped += 1
                continue
            if dropper.packages:
                print(
                    f"{header}: note: removed imports of {', '.join(sorted(dropper.packages))}"
                )

        print(header)
        buffer = io.StringIO()
        sink = buffer if quiet else sys.stdout
        old_cwd = os.getcwd()
        error: BaseException | None = None
        try:
            os.chdir(cwd)
            with contextlib.redirect_stdout(sink), contextlib.redirect_stderr(
                sink
            ):
                exec(compile(tree, str(path), "exec"), namespace)
        except SystemExit as exc:
            if exc.code not in (None, 0):
                error = exc
        except BaseException as exc:  # noqa: BLE001 - report everything
            if isinstance(exc, KeyboardInterrupt):
                raise
            error = exc
        finally:
            os.chdir(old_cwd)

        if error is None:
            passed += 1
            continue
        failed += 1
        if quiet:
            sys.stdout.write(buffer.getvalue())
        print(
            f"*** block {block.number} (lines {block.first_line}-{block.last_line}): FAILED"
        )
        # Drop this runner's own frame from the traceback.
        tb = error.__traceback__.tb_next if error.__traceback__ else None
        traceback.print_exception(type(error), error, tb, file=sys.stdout)
        if stop:
            break

    print(f"=== {path}: {passed} passed, {failed} failed, {skipped} skipped")
    return passed, failed, skipped


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run the ```python blocks of Markdown files."
    )
    parser.add_argument("files", nargs="+", type=Path, help="Markdown files")
    parser.add_argument(
        "--skip-missing",
        action="store_true",
        help="skip code that needs packages that are not installed",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="show output only for failures",
    )
    parser.add_argument(
        "-x",
        "--stop",
        action="store_true",
        help="stop a file at its first failure",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="print the blocks instead of running them",
    )
    parser.add_argument(
        "--cwd", type=Path, help="working directory (default: a temporary one)"
    )
    args = parser.parse_args(argv)

    if args.cwd is not None and not args.cwd.is_dir():
        parser.error(f"--cwd {args.cwd}: no such folder")
    # Paths and block output may hold characters the console encoding
    # lacks (a cp1252 pipe on Windows): never fail on printing them.
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")
    # Agg even if MPLBACKEND names an interactive backend: plt.show()
    # must not open windows and wait.
    os.environ["MPLBACKEND"] = "Agg"
    totals = [0, 0, 0]
    # ignore_cleanup_errors: on Windows a file a block left open cannot be
    # deleted, which must not turn a passing run into a traceback.
    with tempfile.TemporaryDirectory(
        prefix="code-cells-", ignore_cleanup_errors=True
    ) as tmp:
        cwd = (args.cwd or Path(tmp)).resolve()
        for path in args.files:
            if not path.is_file():
                print(f"error: {path}: no such file", file=sys.stderr)
                totals[1] += 1
                continue
            if args.list:
                for block in extract_blocks(path):
                    print(
                        f"--- {path} block {block.number} (lines {block.first_line}-{block.last_line})"
                    )
                    print(block.source, end="")
                continue
            result = run_file(
                path.resolve(),
                skip_missing=args.skip_missing,
                quiet=args.quiet,
                stop=args.stop,
                cwd=cwd,
            )
            totals = [t + r for t, r in zip(totals, result)]
    if not args.list and len(args.files) > 1:
        print(
            f"=== total: {totals[0]} passed, {totals[1]} failed, {totals[2]} skipped"
        )
    return 1 if totals[1] else 0


if __name__ == "__main__":
    sys.exit(main())
