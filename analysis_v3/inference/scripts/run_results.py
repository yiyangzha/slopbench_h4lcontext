"""Every fit, scan and test behind the results of stage 8, run in parallel (idempotent).

    pixi run py -- analysis_v3/inference/scripts/run_results.py --select v5 --tag r1 --model-cat model_cat_v1 \
        --model-incl model_incl_v1 --model-pt4l model_pt4l_v1 --model-njets model_njets_v1 --model-ptj1 model_ptj1_v1 \
        --fiducial production_v3/signal_model/v5/gen_v1/fiducial_xsec.json --yr4 v1 [--workers 8] [--only fits|scans|slow]

The jobs (labels <tag>_<name> under production_v3/inference/<select>/<model>/, the Z -> 4l fit under
production_v3/inference/<select>/<tag>_z4l):
  categories model (m4l with or without the Z1 refit, D_bkg^kin, D_mass; set b unless noted):
    b3d          3D refit, data and Asimov (the benchmark mu and m_H; RESULT.json), set a alongside (a3d); the
                 post-fit Asimov at the data best fit (b3d_asimov_postfit, a3d_asimov_postfit)
    mh_fs        m_H per final state in one fit (mu common) for the mutual compatibility
    table 6      1D L(m4l), 2D L(m4l, D_mass), 3D x refit / no refit (b1d_*, b2dmass_*, b3d, b3d_norefit); the
                 L(m4l, D_bkg^kin) fits b2dkin_* as an extra row
    fs_mh_<fs>   m_H (and mu) of one final state's channels (3D, m_H floating)
    paper-style  paper_mu, fs_mu, category, mode, fv, stxs0: the paper's Eq. 10.1, L(m4l, D_bkg^kin) (m4l without
                 the Z1 refit) at m_H = 125.09 GeV, result set (a) (the benchmark mu is b3d, m_H floating)
    scans        mu and m_H (data and Asimov, total and statistical), m_H in 1D/2D/3D (refit), the
                 (mu_F, mu_V) plane
    width, gof   on-shell width (fit_width.py), saturated goodness of fit with toys (gof.py on b3d data)
  inclusive model (per final state; 1D m4l without the refit, as the paper's fiducial fits):
    fid_int_fix / fid_int_prof   integrated sigma_fid (final-state fractions floating), m_H 125.09 / profiled
    fid_fs_fix   sigma_fid per final state
  differential models: fid_<obs>_fix (pt4l, njets, ptj1; m_H 125.09)
  z4l          the Z -> 4l mass (z4l_mass.py)
Every job runs only if its output is missing.  Logs: production_v3/results/<select>/<tag>/logs/.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import subprocess
import sys
import time
from pathlib import Path

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
SCRIPTS = REPO / "analysis_v3/inference/scripts"
MH_PAPER = 125.09
FINAL_STATES = ("4mu", "4e", "2e2mu")


def jobs(args) -> list[dict]:
    """(name, kind, script, arguments, output) of every job."""
    out = []
    select = args.select

    def fit(model: str, name: str, extra: list, kind: str = "fits", after: str | None = None):
        label = f"{args.tag}_{name}"
        out.append({"name": name, "kind": kind, "script": "fit_model.py", "label": label, "model": model,
                    "args": ["--model", f"{select}/{model}", "--yr4", args.yr4, "--label", label] + extra,
                    "output": PRODUCTION / "inference" / select / model / label / "fit.json", "after": after})

    def both(model: str, name: str, extra: list):
        fit(model, f"{name}_data", extra)
        fit(model, f"{name}_asimov", extra + ["--dataset", "asimov"])

    def scan(model: str, name: str, extra: list, kind: str = "scans"):
        label = f"{args.tag}_{name}"
        out.append({"name": name, "kind": kind, "script": "scan_likelihood.py", "label": label, "model": model,
                    "args": ["--model", f"{select}/{model}", "--yr4", args.yr4, "--label", label, "--workers", str(args.scan_workers)]
                    + extra, "output": PRODUCTION / "inference" / select / model / label / "scan.json"})

    cat = args.model_cat
    both(cat, "b3d", [])
    fit(cat, "b3d_asimov_postfit", ["--dataset", "asimov_postfit", "--postfit-from", f"{args.tag}_b3d_data"],
        after=f"{args.tag}_b3d_data")
    both(cat, "a3d", ["--set", "a"])
    fit(cat, "a3d_asimov_postfit", ["--set", "a", "--dataset", "asimov_postfit", "--postfit-from", f"{args.tag}_a3d_data"],
        after=f"{args.tag}_a3d_data")
    both(cat, "mh_fs", ["--poi-scheme", "mh_fs"])
    # Table 6: 1D L(m4l), 2D L(m4l, D_mass) (the paper's) and 3D, each with and without the Z1 refit; the
    # L(m4l, D_bkg^kin) fit as an extra row.
    for dim, label in (("1D", "b1d"), ("2Dmass", "b2dmass"), ("2D", "b2dkin")):
        for refit in (True, False):
            both(cat, f"{label}_{'refit' if refit else 'norefit'}", ["--dimension", dim] + ([] if refit else ["--no-refit"]))
    both(cat, "b3d_norefit", ["--no-refit"])
    for fs in FINAL_STATES:
        both(cat, f"fs_mh_{fs}", ["--final-state", fs])
    # Signal strengths as in the paper (Eq. 10.1): L(m4l, D_bkg^kin) with m4l without the Z1 refit, at m_H = 125.09 GeV,
    # the paper-style result set (a).
    paper = ["--dimension", "2D", "--no-refit", "--set", "a", "--fix-mh", str(MH_PAPER)]
    both(cat, "paper_mu", paper)
    both(cat, "fs_mu", ["--poi-scheme", "final_state"] + paper)
    for scheme in ("category", "mode", "fv", "stxs0"):
        both(cat, scheme, ["--poi-scheme", scheme] + paper)
    for dataset in ("data", "asimov"):
        for stat in (False, True):
            suffix = f"{dataset}{'_stat' if stat else ''}"
            common = ["--dataset", dataset] + (["--stat-only"] if stat else [])
            scan(cat, f"scan_mu_{suffix}", ["--x", "mu", "0.2", "1.8", "33"] + common)
            scan(cat, f"scan_mh_{suffix}", ["--x", "mH", "123.5", "126.5", "31"] + common)
        for dim, label in (("1D", "1d"), ("2Dmass", "2dmass")):
            scan(cat, f"scan_mh_{label}_{dataset}", ["--x", "mH", "123.5", "126.5", "31", "--dimension", dim, "--dataset", dataset])
        scan(cat, f"scan_fv_{dataset}", ["--x", "mu_F", "0.0", "2.8", "15", "--y", "mu_V", "0.0", "6.0", "16", "--poi-scheme", "fv",
                                         "--dimension", "2D", "--no-refit", "--set", "a", "--fix-mh", str(MH_PAPER), "--dataset", dataset])
    out.append({"name": "width", "kind": "slow", "script": "fit_width.py", "label": f"{args.tag}_width", "model": cat,
                "args": ["--model", f"{select}/{cat}", "--yr4", args.yr4, "--label", f"{args.tag}_width"],
                "output": PRODUCTION / "inference" / select / cat / f"{args.tag}_width" / "width.json"})
    out.append({"name": "gof", "kind": "slow", "script": "gof.py", "label": f"{args.tag}_gof", "model": cat,
                "args": ["--model", f"{select}/{cat}", "--fit", f"{args.tag}_b3d_data", "--yr4", args.yr4, "--label",
                         f"{args.tag}_gof", "--toys", str(args.gof_toys), "--workers", str(args.scan_workers)],
                "output": PRODUCTION / "inference" / select / cat / f"{args.tag}_gof" / "gof.json", "after": f"{args.tag}_b3d_data"})
    fiducial = ["--fiducial", str(args.fiducial), "--dimension", "1D", "--no-refit"]
    both(args.model_incl, "fid_int_fix", ["--poi-scheme", "fid_int", "--fix-mh", str(MH_PAPER)] + fiducial)
    both(args.model_incl, "fid_int_prof", ["--poi-scheme", "fid_int"] + fiducial)
    both(args.model_incl, "fid_fs_fix", ["--poi-scheme", "fid_fs", "--fix-mh", str(MH_PAPER)] + fiducial)
    for obs, model in (("pt4l", args.model_pt4l), ("njets", args.model_njets), ("ptj1", args.model_ptj1)):
        both(model, f"fid_{obs}_fix", ["--poi-scheme", f"fid_{obs}", "--fix-mh", str(MH_PAPER)] + fiducial)
    out.append({"name": "z4l", "kind": "fits", "script": "z4l_mass.py", "label": f"{args.tag}_z4l", "model": None,
                "args": ["--select", select, "--label", f"{args.tag}_z4l"],
                "output": PRODUCTION / "inference" / select / f"{args.tag}_z4l" / "z4l_mass.json"})
    return out


def run(job: dict, log_dir: Path) -> dict:
    if job["output"].exists():
        return {"name": job["name"], "status": "exists"}
    log = log_dir / f"{job['label']}.log"
    started = time.time()
    with log.open("a", encoding="utf-8") as stream:
        result = subprocess.run([sys.executable, str(SCRIPTS / job["script"])] + job["args"], stdout=stream, stderr=subprocess.STDOUT,
                                cwd=REPO)
    status = "ok" if result.returncode == 0 and job["output"].exists() else f"failed ({result.returncode})"
    return {"name": job["name"], "status": status, "seconds": round(time.time() - started, 1), "log": str(log)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--select", required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--model-cat", required=True)
    parser.add_argument("--model-incl", required=True)
    parser.add_argument("--model-pt4l", required=True)
    parser.add_argument("--model-njets", required=True)
    parser.add_argument("--model-ptj1", required=True)
    parser.add_argument("--fiducial", type=Path, required=True)
    parser.add_argument("--yr4", required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--scan-workers", type=int, default=4)
    parser.add_argument("--rest-parallel", type=int, default=3, help="scans and slow jobs run at the same time")
    parser.add_argument("--gof-toys", type=int, default=500)
    parser.add_argument("--only", choices=["fits", "scans", "slow"], default=None)
    args = parser.parse_args()
    args.fiducial = args.fiducial if args.fiducial.is_absolute() else REPO / args.fiducial
    log_dir = PRODUCTION / "results" / args.select / args.tag / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    todo = [j for j in jobs(args) if args.only is None or j["kind"] == args.only]
    # Single fits in parallel (those that start from another fit's result after it); then the scans and the slow
    # jobs (each with its own workers), --rest-parallel at a time.
    fits = [j for j in todo if j["kind"] == "fits"]
    rest = [j for j in todo if j["kind"] != "fits"]
    results = []
    for phase in ([j for j in fits if not j.get("after")], [j for j in fits if j.get("after")]):
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            for r in pool.map(lambda j: run(j, log_dir), phase):
                print(f"[results] {r}", flush=True)
                results.append(r)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.rest_parallel) as pool:
        for r in pool.map(lambda j: run(j, log_dir), rest):
            print(f"[results] {r}", flush=True)
            results.append(r)
    failed = [r for r in results if r["status"].startswith("failed")]
    print(f"[results] {len(results)} jobs, {len(failed)} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
