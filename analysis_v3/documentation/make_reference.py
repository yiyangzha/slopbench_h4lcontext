"""REFERENCE.md of the GitHub branch `reference`: the layout of the branch and the key numbers, from the result files.

    pixi run py -- analysis_v3/documentation/make_reference.py --out <file> --an <AN directory name> \
        --main production_v3/results/v5/final_v2 --truth production_v3/results/v5/truth_v2/truth.json \
        --eval "half A=<run dir>" "half B=<run dir>" "full=<run dir>"
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")


def pois(result: dict) -> dict:
    return {p["name"]: p for p in result["pois"]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--an", required=True)
    parser.add_argument("--main", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--eval", nargs="+", required=True, help="label=<run directory> (branch results/<name> in order A, B, full)")
    args = parser.parse_args()
    main_result = json.loads((REPO / args.main / "RESULT.json").read_text(encoding="utf-8"))
    truth = json.loads((REPO / args.truth).read_text(encoding="utf-8"))["compare"]
    names = ["eval_test_10fb_A", "eval_test_10fb_B", "eval_full_20fb"]
    runs = []
    for (item, name) in zip(args.eval, names):
        label, path = item.split("=", 1)
        r = json.loads((REPO / path / "RESULT.json").read_text(encoding="utf-8"))
        full = json.loads((REPO / path / "results_full.json").read_text(encoding="utf-8"))
        runs.append((label, name, r, full))
    m = pois(main_result)
    lines = [
        "# Reference analysis: H -> ZZ* -> 4l on UL16 PFNano pseudo-data (branch `reference`)",
        "",
        "## Evaluation submission (EVAL_CONTRACT v0.4, task h4l_ntuple)",
        "",
        "    ./run.sh <dataset_dir> <output_dir>                  # from the repository root (forwards to my_analysis/run.sh)",
        "    my_analysis/run.sh <dataset_dir> <output_dir>        # the submission directory itself",
        "    tests/run_test.sh <dataset_dir> <output_dir>         # run + check_outputs.py (contract checks, pyhf build of MODEL.json)",
        "",
        "`my_analysis/` is the self-contained submission (run.sh, pixi.toml + pixi.lock with root_base and the C++ compiler,",
        "EVAL_CONTRACT.md with the per-flavour calibration block, the code in h4l_eval/ and src/).  It carries out the main",
        "analysis of this repository with the same methods and parameters in every step; the only simplification is the",
        "tag-and-probe nominal fit model (see my_analysis/README.md and the AN section \"Evaluation submission\").",
        "`<dataset_dir>` has the DATASET.md layout (data/part_*.root, mc/<process>/*.root + meta.json,",
        "mc/cross_sections.json, lumi.json, meta.json) with ROOT files in the UL16 PFNano format.  Outputs: RESULT.json,",
        "MODEL.json (pyhf workspace), MODEL_meta.json, summary.txt, results_full.json, logs/run.log (and the intermediate",
        "reader files in work/).",
        "",
        "`RESULT.json` and `MODEL.json` at the repository root are the results of the main analysis on the 20 fb^-1",
        "pseudo-data (identical to results/main_final/).",
        "",
        "## Contents",
        "",
        "* `my_analysis/` — the evaluation submission.",
        "* `analysis_v3/` — the main analysis code (C++ programs with Makefiles, Python stages, configurations); the binaries",
        "  are built with `make`.  `analysis_v3/eval/` builds the submission's method constants and the test datasets.",
        "* `results/main_final/` — the main analysis on the 20 fb^-1 pseudo-data (RESULT.json with the per-flavour",
        "  calibration, MODEL.json, results_full.json, summary.md, tables, plots, the freeze record).",
    ]
    for label, name, r, full in runs:
        t = full.get("timings_s", {}).get("total", float("nan"))
        lines.append(f"* `results/{name}/` — the submission run on {label} ({full['lumi_fb']:.2f} fb^-1 in the DATASET.md layout; "
                     f"{t / 60:.0f} min on 16 cores with the input on EOS).")
    lines += [
        "* `results/truth/` — the truth of the pseudo-data (injection record and generator profiles, make_truth.py).",
        f"* `deliverables/{args.an}/` — the analysis note (main.pdf, main.tex; regenerated with",
        "  analysis_v3/documentation/make_an.py), including the truth comparison and the evaluation submission.",
        "* `AGENTS.md`, `PLAN.md`, `README.md`, `experiment_log.md` — rules, plan, workflow and the full decision log.",
        "",
        "## Key numbers",
        "",
        f"Truth: mu = {truth['mu']:g}, m_H = {truth['mH']:g} GeV; scale_shift muon {truth['scale_shift muon']:.5f}, electron "
        f"{truth['scale_shift electron']:.5f}; smear {truth['smear muon']:g} (pT part; the injected angular smears are absorbed by",
        f"the measured smear); sel_eff muon {truth['sel_eff muon']:.5f}, electron {truth['sel_eff electron']:.5f} (data-weighted).",
        "",
        "| analysis | L [fb^-1] | mu | m_H [GeV] | Z obs (exp) | GoF p | scale_shift mu / e | smear mu / e | sel_eff mu / e |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    def row(label, lumi, r):
        p = pois(r)
        c = r["calibration"]
        q = r.get("quality", {})
        return (f"| {label} | {lumi:.2f} | {p['mu']['value']:.3f} +- {p['mu']['total']:.3f} | {p['mH']['value']:.2f} +- {p['mH']['total']:.2f} | "
                f"{r['significance_obs']:.2f} ({r.get('significance_exp', float('nan')):.2f}) | {q.get('gof_pvalue', {}).get('value', float('nan')):.2f} | "
                f"{c['scale_shift']['muon']['value']:.5f} / {c['scale_shift']['electron']['value']:.5f} | "
                f"{c['smear']['muon']['value']:.4f} / {c['smear']['electron']['value']:.4f} | "
                f"{c['sel_eff']['muon']['value']:.4f} / {c['sel_eff']['electron']['value']:.4f} |")
    lines.append(row("main analysis (3D, categories)", 20.0, main_result))
    for label, name, r, full in runs:
        lines.append(row(f"submission, {label}", full["lumi_fb"], r))
    lines += ["", f"Main analysis: mu = {m['mu']['value']:.3f} +- {m['mu']['total']:.3f}, m_H = {m['mH']['value']:.2f} +- {m['mH']['total']:.2f} GeV "
              "(set b, 3D fit with m_H floating); the submission fits m4l x D_mass binned without MELA categories.", ""]
    args.out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[reference] {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
