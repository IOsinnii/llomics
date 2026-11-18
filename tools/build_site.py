#!/usr/bin/env python3
"""Assemble the browser demo into one static folder, ready for GitHub Pages.

    python tools/build_site.py                 # writes _site/
    python tools/build_site.py --serve         # writes _site/ and serves it at http://localhost:8000
    python tools/build_site.py --out other/    # another folder
    python tools/build_site.py --no-pyodide    # skip Pyodide (the page will not start; for tests)

The page (web/index.html, app.js, worker.js, style.css) goes to the top of the folder. The Python
it runs goes under repo/, in the repository's own layout, so web/bridge.py, benchmark/score.py and
the judge script find each other exactly as they do in a checkout. repo/manifest.json lists those
files for the worker to load. Nothing is copied that the page does not use.

Pyodide (CPython compiled to WebAssembly, about 13 MB) is served from the site itself under
pyodide/, not from a CDN at run time: the page then depends on no other host, and a visitor's
browser talks to GitHub Pages only. The core files of one pinned release are downloaded once from
the jsDelivr CDN into .cache/ and checked against the hashes below before they are copied.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import http.server
import json
import shutil
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAGE = ["index.html", "app.js", "worker.js", "style.css"]
REPO_GLOBS = [
    "web/bridge.py",
    "python/llomics/*.py",
    "python/llomics/*.json",
    "skills/llomics-judge/scripts/judge_bedrock.py",
    "skills/llomics-limma/references/questions.contrast.json",
    "benchmark/score.py",
    "benchmark/scenarios.json",
    "benchmark/payloads/*.llomics.json",
    "benchmark/answers/*/*.answers.json",
    "benchmark/data/gene_symbols.txt",
    "examples/pbmc10k/payload.json",
]

# To move to another Pyodide release: change the version, delete the hashes, run the build once to
# see the new ones, and check the sizes against https://www.jsdelivr.com/package/npm/pyodide.
PYODIDE_VERSION = "314.0.8"
PYODIDE_URL = f"https://cdn.jsdelivr.net/pyodide/v{PYODIDE_VERSION}/full/"
PYODIDE_FILES = {
    "pyodide.mjs": "1bb8cc14f9405aaa7cd4e00c12edd0a8539de9cc026618d6ee5a9663f367dc70",
    "pyodide.asm.mjs": "22514b8be6980e28a02c8b1d6387ce8dfe666fc4edf8d547cff5dce29ff4d972",
    "pyodide.asm.wasm": "e14d49ceca9b4a0e274ea918f3a898de2d3880b6e5573a6cfe1a2267ecc70329",
    "python_stdlib.zip": "2bd3394a145bfbbe5bbc9d9d6ab9136abbb986c6fef6c142b3a2962227aa5c28",
    "pyodide-lock.json": "5dc2fc119108bc148c7457dc86e7675b5c87e1cafd420b9c34c1eaef7b36c010",
}
CACHE = ROOT / ".cache" / "pyodide" / PYODIDE_VERSION


def repo_files() -> list[str]:
    files = sorted({p.relative_to(ROOT).as_posix() for pattern in REPO_GLOBS for p in ROOT.glob(pattern) if p.is_file()})
    return [f for f in files if "__pycache__" not in f]


def pyodide_files() -> Path:
    """The pinned Pyodide core files, downloaded into the cache if needed and checked against their hashes."""
    CACHE.mkdir(parents=True, exist_ok=True)
    for name, expected in PYODIDE_FILES.items():
        path = CACHE / name
        if not path.exists():
            print(f"downloading {PYODIDE_URL}{name}")
            with urllib.request.urlopen(PYODIDE_URL + name, timeout=120) as response:
                data = response.read()
            tmp = path.with_suffix(path.suffix + ".part")
            tmp.write_bytes(data)
            tmp.replace(path)
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            path.unlink()
            raise SystemExit(f"{name}: sha256 {actual} does not match the pinned {expected}; the cached copy was removed")
    return CACHE


def build(out: Path, with_pyodide: bool = True) -> list[str]:
    if out.resolve() == ROOT.resolve():
        raise SystemExit("refusing to build into the repository root")
    source = pyodide_files() if with_pyodide else None
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    for name in PAGE:
        shutil.copy2(ROOT / "web" / name, out / name)
    files = repo_files()
    for rel in files:
        target = out / "repo" / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / rel, target)
    (out / "repo" / "manifest.json").write_text(json.dumps({"files": files}, indent=1) + "\n", encoding="utf-8")
    if source:
        (out / "pyodide").mkdir()
        for name in PYODIDE_FILES:
            shutil.copy2(source / name, out / "pyodide" / name)
    (out / ".nojekyll").write_text("", encoding="utf-8")  # serve files as they are, no Jekyll pass
    return files


class Handler(http.server.SimpleHTTPRequestHandler):
    # Module scripts need a JavaScript type, and WebAssembly streams only as application/wasm.
    extensions_map = {**http.server.SimpleHTTPRequestHandler.extensions_map,
                      ".js": "text/javascript", ".mjs": "text/javascript", ".wasm": "application/wasm"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=ROOT / "_site")
    parser.add_argument("--serve", action="store_true", help="serve the folder after building it")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--no-pyodide", action="store_true", help="leave Pyodide out (no download)")
    args = parser.parse_args()

    files = build(args.out, with_pyodide=not args.no_pyodide)
    size = sum(p.stat().st_size for p in args.out.rglob("*") if p.is_file())
    print(f"Built {args.out}: {len(files)} files for Python, {size / 1e6:.1f} MB in all")
    if args.serve:
        handler = functools.partial(Handler, directory=str(args.out))
        with http.server.ThreadingHTTPServer(("127.0.0.1", args.port), handler) as server:
            print(f"Serving at http://localhost:{args.port}  (Ctrl+C to stop)")
            try:
                server.serve_forever()
            except KeyboardInterrupt:
                pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
