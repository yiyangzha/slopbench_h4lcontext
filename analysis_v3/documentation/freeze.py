"""Freeze record of the code, the configuration and the results before unblinding (AGENTS.md "Freeze and unblind").

    pixi run py -- analysis_v3/documentation/freeze.py --results production_v3/results/v5/r1b --out freeze_v1

The sha256 of every source file of analysis_v3/ (code and configuration), of the frozen programs used by the final
stages and of every file of the results directory, written to <results>/<out>.json (never overwritten); the
sha256 of that record goes into experiment_log.md by hand.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
SOURCE_SUFFIXES = {".py", ".cpp", ".h", ".json", ".sh", ".md", ".txt"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    results = args.results if args.results.is_absolute() else REPO / args.results
    target = results / f"{args.out}.json"
    if target.exists():
        raise SystemExit(f"{target} exists")
    code = {str(p.relative_to(REPO)): sha256(p) for p in sorted((REPO / "analysis_v3").rglob("*"))
            if p.is_file() and p.suffix in SOURCE_SUFFIXES and "__pycache__" not in p.parts}
    binaries = {str(p.relative_to(REPO)): sha256(p) for p in sorted((REPO / "analysis_v3").rglob("bin/*")) if p.is_file()}
    outputs = {str(p.relative_to(REPO)): sha256(p) for p in sorted(results.rglob("*")) if p.is_file() and p != target}
    record = {"schema": "h4l_v3_freeze/1", "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
              "results": str(results), "code": code, "binaries": binaries, "results_files": outputs,
              "documents": {name: sha256(REPO / name) for name in ("AGENTS.md", "PLAN.md", "README.md")}}
    target.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[freeze] {len(code)} source files, {len(binaries)} binaries, {len(outputs)} result files\n[freeze] {target} "
          f"sha256 {sha256(target)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
