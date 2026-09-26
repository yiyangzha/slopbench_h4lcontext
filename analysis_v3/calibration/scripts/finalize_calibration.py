"""Final payload of a calibration run that met the iteration limit (stage 3a).

    pixi run py -- analysis_v3/calibration/scripts/finalize_calibration.py --version v3 --run nominal_r3 [--last 4]

run_calibration.py stops when every parameter's residual is below
converged_sigma of its uncertainty.  Parameters constrained by few, strongly
sculpted categories (soft legs) can instead oscillate between iterations at
the level of their statistical uncertainty.  For a run that reached
max_iterations this script defines the final payload as the mean of the
payloads applied in the last N iterations (a, b, c, d; the pT nodes of the
last iteration), and records per parameter the root-mean-square deviation of
those payloads from their mean ("iteration_rms"), to be added in quadrature to
the statistical uncertainty of the last iteration's joint least squares.  It
writes <run>/payload.json and <run>/summary_final.json; nothing is overwritten.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/calibration/scripts"))
from run_calibration import FLAVOURS, Model, write_new  # noqa: E402


def flatten(model: Model, payload: dict, eta_key: str, pt_key: str) -> np.ndarray:
    values = np.zeros(model.n_par)
    values[: model.n_eta] = payload[eta_key]
    for (r, p), k in model.b_index.items():
        values[k] = payload[pt_key][r][p]
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--last", type=int, default=4)
    args = parser.parse_args()
    base = PRODUCTION / "calibration" / args.version / args.run
    summary = json.loads((base / "summary.json").read_text(encoding="utf-8"))
    if summary.get("converged"):
        print("ERROR: the run converged; its payload.json is final", file=sys.stderr)
        return 1
    config = json.loads((base / "config.json").read_text(encoding="utf-8"))
    last = summary["iterations"][-1]["iteration"]
    # The payload applied in iteration n is iter_{n-1}/payload.json.
    applied = [last - k for k in range(args.last)]
    if min(applied) < 2:
        print("ERROR: too few iterations after the frozen category selection", file=sys.stderr)
        return 1
    payloads = [json.loads((base / f"iter_{n - 1:02d}" / "payload.json").read_text(encoding="utf-8")) for n in applied]
    final, rms = {}, {}
    for tag, name in FLAVOURS.items():
        model = Model(config["model"][name])
        stacks = {q: np.array([flatten(model, p[name], e, t) for p in payloads])
                  for q, (e, t) in (("scale", ("a", "b")), ("smear", ("c", "d")))}
        mean = {q: stacks[q].mean(axis=0) for q in stacks}
        rms[tag] = {q: np.sqrt(((stacks[q] - mean[q]) ** 2).mean(axis=0)).tolist() for q in stacks}
        entry = {"model": payloads[0][name]["model"], "node_pt": payloads[0][name]["node_pt"]}
        for q, (eta_key, pt_key) in (("scale", ("a", "b")), ("smear", ("c", "d"))):
            entry[eta_key] = mean[q][: model.n_eta].tolist()
            rows = [[0.0] * model.n_pt for _ in range(model.n_regions)]
            for (r, p), k in model.b_index.items():
                rows[r][p] = float(mean[q][k])
            entry[pt_key] = rows
        final[name] = entry
    solution_path = base / f"iter_{last:02d}" / "solution.json"
    write_new(base / "payload.json", {"schema": "h4l_v3_lepton_calibration_payload/2", "run": args.run,
                                      "finalized_by_average": True, "averaged_applied_iterations": applied,
                                      "solution": str(solution_path), "iteration_rms": rms, "flavours": final})
    write_new(base / "summary_final.json", {"iterations": summary["iterations"], "converged": False,
                                            "finalized_by_average": True, "averaged_applied_iterations": applied,
                                            "converged_iteration": last, "payload": str(base / "payload.json"),
                                            "finished_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")})
    worst = {tag: {q: float(max(v)) for q, v in r.items()} for tag, r in rms.items()}
    print(f"[finalize] {base / 'payload.json'}; largest iteration rms {worst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
