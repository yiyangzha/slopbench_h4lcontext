"""Z tag-and-probe efficiencies and scale factors (stage 3c).

    pixi run py -- analysis_v3/tnp/scripts/run_tnp.py --extract v2 --calibration v4/nominal --run v6 \
        [--config analysis_v3/tnp/config/tnp_ul16_v7.json] [--fit-label fits_i2]
    pixi run py -- analysis_v3/tnp/scripts/run_tnp.py --extract v2 --calibration v4/nominal --run closure_v1_k2 \
        --closure analysis_v3/tnp/config/tnp_closure_ul16_v1.json --closure-point k2

Inputs: the extraction scan production_v3/tnp/<extract>/extract/scan.json, the
calibration payload production_v3/calibration/<calibration>/payload.json, the
configuration (default tnp_ul16_v7.json) and the skim trigger bit order
(production_v3/skims/v2/plan.json).  The programs tnp_histograms and fit_tnp
are frozen by content under production_v3/program_inputs/tnp_v3/.
 1. histograms (production_v3/tnp/<extract>/<run>/hist/): one job per
    extraction output (data corrected and MC smeared with the calibration
    payload; MC weighted with sigma L / sum genEventSumw; the DY-MC jobs also
    fill the signal templates of the prompt-prompt pairs).  The jobs embed
    only the histogram part of the configuration (tag, probe, steps, bins,
    histogram, Z mass), so that the fit settings and the per-bin overrides
    can be iterated without redoing them; <run>/run.json binds the inputs,
    that part of the configuration and tnp_histograms.
 2. fits (<run>/<fit-label>/): one job per flavour and step, the inclusive
    bin first, then every probe bin, for the data and the MC, with the
    nominal and the two alternative models; <fit-label>/fit_run.json binds
    the histograms, the fit configuration (incl. overrides) and fit_tnp.  A
    new inspection round of the galleries is a new --fit-label.
 3. scale factors per step and bin (<run>/<fit-label>/sf.json): SF =
    eps_data / eps_MC, both from the tag-and-probe fits (nominal model;
    counting only where a fit is unusable; the MC generator truth never
    enters), statistical error from the fit errors, fit-model systematic =
    the largest |SF_alt - SF_nominal| (alternative models applied to data and
    MC alike); the full single-lepton SF = SF_id SF_sip SF_iso (the chain),
    with the directly measured full step as the cross-check; the full SF at
    pT = 45 GeV, |eta| = 1.2 by bilinear interpolation between the centres of
    the neighbouring bins.
Nothing is overwritten; finished identical jobs are reused.  An existing run
keeps the frozen tnp_histograms that its run.json records.

Efficiency closure (--closure, --closure-point): the data role is the MC with
even original file keys (every MC sample of the extraction, normalized to the
genEventSumw of its even files), the MC role the odd ones (normalized to
theirs; the DY-MC templates from this role), both calibrated as MC; in the data
role a passed identification is dropped with probability 1 - keep_id per
flavour (tnp_histograms closure mode), a known SF_id = SF_full = keep_id.  The
closure configuration sets the (coarser) bins; the per-bin fit overrides and
exclusions of the configuration refer to its own bins and are dropped.
<fit-label>/closure.json compares the recovered SFs with keep_id.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import sys
from pathlib import Path

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
sys.path.insert(0, str(REPO / "analysis_v3/calibration/scripts"))
from run_calibration import PRODUCTION, canonical, fail, parallel, run_job, sha256_file, utc, write_new  # noqa: E402

CONFIG = REPO / "analysis_v3/tnp/config/tnp_ul16_v7.json"
PROGRAMS = {name: REPO / "analysis_v3/tnp/bin" / name for name in ("tnp_histograms", "fit_tnp")}
FLAVOURS = {"mm": "muon", "ee": "electron"}
# The signal templates of the fits come from the prompt-prompt pairs of the DY MC (user decision 2026-09-24).
TEMPLATE_SAMPLE = "DYJetsToLL"
STEPS = ("id", "sip", "iso", "full")
CHAIN = ("id", "sip", "iso")
MODELS = ("nominal", "alt_signal", "alt_background")
REPORT_POINT = (45.0, 1.2)


def freeze(names: tuple) -> dict:
    digests = {name: sha256_file(PROGRAMS[name]) for name in names}
    key = hashlib.sha256("".join(digests[n] for n in sorted(digests)).encode()).hexdigest()[:12]
    directory = PRODUCTION / "program_inputs" / "tnp_v3" / key
    directory.mkdir(parents=True, exist_ok=True)
    frozen = {}
    for name in names:
        source = PROGRAMS[name]
        target = directory / name
        if not target.exists():
            temporary = target.with_name(name + f".partial.{os.getpid()}")
            shutil.copy2(source, temporary)
            if sha256_file(temporary) != digests[name]:
                fail(f"frozen copy of {name} differs")
            temporary.rename(target)
        elif sha256_file(target) != digests[name]:
            fail(f"frozen {target} differs from its content address")
        os.chmod(target, 0o755)
        frozen[name] = {"path": str(target), "sha256": digests[name], "source": str(source)}
    return {"key": key, "programs": frozen}


def bind(path: Path, record: dict) -> None:
    """Create the binding record, or check that an existing one binds the same things."""
    if path.exists():
        existing = json.loads(path.read_text(encoding="utf-8"))
        if canonical({k: v for k, v in existing.items() if k != "started_utc"}) != canonical(record):
            fail(f"{path} records other inputs, settings or programs; use a new name")
    else:
        write_new(path, dict(record, started_utc=utc()))


def interpolate(values: dict, errors: dict, pt_centres: list, eta_centres: list, pt: float, eta: float) -> dict:
    """Bilinear interpolation on the grid of bin centres (clamped at the outermost centres)."""
    def bracket(centres, x):
        if x <= centres[0]:
            return 0, 0, 0.0
        if x >= centres[-1]:
            return len(centres) - 1, len(centres) - 1, 0.0
        for i in range(len(centres) - 1):
            if centres[i] <= x <= centres[i + 1]:
                return i, i + 1, (x - centres[i]) / (centres[i + 1] - centres[i])
        raise ValueError("unsorted centres")
    i0, i1, tp = bracket(pt_centres, pt)
    j0, j1, te = bracket(eta_centres, eta)
    weights = {}
    for i, wi in ((i0, 1 - tp), (i1, tp)):
        for j, wj in ((j0, 1 - te), (j1, te)):
            weights[(i, j)] = weights.get((i, j), 0.0) + wi * wj
    value, variance = 0.0, 0.0
    for (i, j), w in weights.items():
        if w == 0:
            continue
        if values.get((i, j)) is None:
            return {"value": None}
        value += w * values[(i, j)]
        variance += (w * errors[(i, j)]) ** 2
    return {"value": value, "error": math.sqrt(variance), "weights": {f"{i},{j}": w for (i, j), w in weights.items() if w}}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--extract", required=True)
    parser.add_argument("--calibration", required=True, help="<version>/<run> of the calibration payload")
    parser.add_argument("--run", required=True)
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--fit-label", default="fits")
    parser.add_argument("--manifests", default="v2")
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--histograms-only", action="store_true", help="stop after the histogram jobs")
    parser.add_argument("--closure", type=Path, default=None, help="efficiency-closure configuration")
    parser.add_argument("--closure-point", default=None, help="the closure point (keep_id per flavour) to run")
    parser.add_argument("--condor", choices=["plan", "merge"], default=None,
                        help="plan: write the fits as a Condor plan of bin chunks (then stage it); merge: collect the "
                             "chunk outputs of that plan and assemble the scale factors")
    parser.add_argument("--chunk-bins", type=int, default=4, help="probe bins per Condor chunk (each chunk refits 'all' first)")
    args = parser.parse_args()
    config_path = args.config if args.config.is_absolute() else REPO / args.config
    base = PRODUCTION / "tnp" / args.extract / args.run
    config = json.loads(config_path.read_text(encoding="utf-8"))
    closure, point = None, None
    if args.closure:
        closure_path = args.closure if args.closure.is_absolute() else REPO / args.closure
        closure = json.loads(closure_path.read_text(encoding="utf-8"))
        if args.closure_point not in closure["points"]:
            fail(f"--closure-point must be one of {sorted(closure['points'])}")
        point = closure["points"][args.closure_point]
        config = dict(config, bins=closure["bins"])
        config["fit"] = {k: v for k, v in config["fit"].items() if k not in ("overrides", "exclude_alternatives")}
    elif args.closure_point:
        fail("--closure-point needs --closure")
    scan_path = PRODUCTION / "tnp" / args.extract / "extract" / "scan.json"
    scan = json.loads(scan_path.read_text(encoding="utf-8"))
    if canonical(config["extract"]) != canonical(scan["extract_config"]):
        fail("the extraction section of the configuration differs from the scanned extraction")
    payload_path = PRODUCTION / "calibration" / args.calibration / "payload.json"
    payload = json.loads(payload_path.read_text(encoding="utf-8"))["flavours"]
    skim_plan = json.loads((PRODUCTION / "skims" / "v2" / "plan.json").read_text(encoding="utf-8"))
    paths = skim_plan["tasks"][0]["analysis_config"]["triggers"]["paths"]
    trigger_bits = {name: index for index, name in enumerate(paths)}
    manifests = PRODUCTION / "manifests" / args.manifests
    lumi = scan["samples"]["pseudo_data"]["lumi_fb"]
    scales, half_scales = {}, {}
    for name, sample in scan["samples"].items():
        if sample["kind"] == "mc":
            manifest = json.loads((manifests / f"mc_{name}.json").read_text(encoding="utf-8"))
            sumw = sum(item["genEventSumw"] for item in manifest["files"])
            scales[name] = sample["sigma_eff_pb"] * 1000.0 * lumi / sumw
            # Closure halves: normalized to the genEventSumw of the files of their own original-file-key parity.
            half_sumw = [sum(item["genEventSumw"] for item in manifest["files"] if int(item["file_key"], 16) % 2 == parity)
                         for parity in (0, 1)]
            half_scales[name] = [sample["sigma_eff_pb"] * 1000.0 * lumi / s if s > 0 else None for s in half_sumw]
    hist_config = {key: config[key] for key in ("schema", "version", "tag", "probe", "steps", "bins", "histogram")}
    hist_config["fit"] = {"z_mass": config["fit"]["z_mass"]}
    run_record = base / "run.json"
    if run_record.exists():
        # An existing run keeps the frozen histogram program it was made with.
        frozen_hist = json.loads(run_record.read_text(encoding="utf-8"))["frozen"]
        for name, item in frozen_hist["programs"].items():
            if sha256_file(Path(item["path"])) != item["sha256"]:
                fail(f"the frozen {name} of {run_record} changed")
    else:
        frozen_hist = freeze(("tnp_histograms",))
    record = {"extract": args.extract, "extract_scan": str(scan_path), "extract_scan_sha256": sha256_file(scan_path),
              "calibration_payload": str(payload_path), "calibration_payload_sha256": sha256_file(payload_path),
              "histogram_config": hist_config, "frozen": frozen_hist, "mc_scales": scales}
    if closure:
        record["closure"] = {"config": str(closure_path), "config_sha256": sha256_file(closure_path), "point": args.closure_point,
                             "keep_id": point, "half_scales": half_scales}
    bind(run_record, record)

    print(f"[{utc()}] histograms", flush=True)
    tasks, files = [], {"data": [], "mc": []}
    for output in scan["outputs"]:
        inputs = [{"path": output["path"], "tnp_pairs": output["tnp_pairs"], "root_sha256": output["root_sha256"]}]
        if closure:
            if output["kind"] != "mc":
                continue
            specs = []
            for role, parity, keep in (("data", 0, point), ("mc", 1, {"muon": 1.0, "electron": 1.0})):
                scale = half_scales[output["sample"]][parity]
                if scale is None:
                    continue
                specs.append((role, {"tnp_config": hist_config, "role": "mc", "weight_scale": scale, "trigger_bits": trigger_bits,
                                     "payload": payload, "fill_templates": role == "mc" and output["sample"] == TEMPLATE_SAMPLE,
                                     "closure": {"parity": parity, "keep_id": keep}, "inputs": inputs}))
        else:
            role = "data" if output["kind"] == "data" else "mc"
            specs = [(role, {"tnp_config": hist_config, "role": role,
                             "weight_scale": 1.0 if role == "data" else scales[output["sample"]], "trigger_bits": trigger_bits,
                             "payload": payload, "fill_templates": output["sample"] == TEMPLATE_SAMPLE, "inputs": inputs})]
        for role, job in specs:
            name = f"{role}_{output['task_id']}"
            out = base / "hist" / f"{name}.root"
            tasks.append((run_job, frozen_hist["programs"]["tnp_histograms"], job, base / "jobs" / f"{name}.json", out,
                          base / "logs" / f"{name}.log"))
            files[role].append(str(out))
    records = parallel(tasks, args.workers)
    histograms = {Path(r["output"]).name: r["output_sha256"] for r in records}
    if args.histograms_only:
        print(f"[{utc()}] histograms done ({len(records)} outputs)", flush=True)
        return 0

    fit_dir = base / args.fit_label
    frozen_fit = freeze(("fit_tnp",))
    binding = {"histograms": histograms, "config": str(config_path), "config_sha256": sha256_file(config_path),
               "fit_config": config["fit"], "frozen": frozen_fit, "orchestrator_sha256": sha256_file(Path(__file__))}
    if args.condor:
        binding["backend"] = {"condor": True, "chunk_bins": args.chunk_bins}
    bind(fit_dir / "fit_run.json", binding)
    print(f"[{utc()}] fits ({args.fit_label})", flush=True)
    tasks, jobs = [], {}
    for tag, name in FLAVOURS.items():
        n_bins = len(config["bins"][name]["pt_edges"]) * len(config["bins"][name]["eta_edges"])
        bins = ["all"] + [str(b) for b in range(n_bins)]
        for step in STEPS:
            job = {"tnp_config": {"fit": config["fit"]}, "flavour": tag, "step": step, "bins": bins, "data": files["data"],
                   "mc": files["mc"], "curves": True, "n_eta": len(config["bins"][name]["eta_edges"]),
                   "pt_edges": config["bins"][name]["pt_edges"]}
            label = f"fit_{tag}_{step}"
            jobs[label] = job
            tasks.append((run_job, frozen_fit["programs"]["fit_tnp"], job, fit_dir / "jobs" / f"{label}.json",
                          fit_dir / "fits" / f"{label}.json", fit_dir / "logs" / f"{label}.log"))
    plan_path = fit_dir / "condor_plan.json"
    stage = f"tnp_fits_{args.run}_{args.fit_label}".replace("-", "_")
    if args.condor == "plan":
        # Every bin is independent; the per-bin fits start from the inclusive fit, so each chunk fits "all" first and
        # reproduces the serial fits of its bins exactly.
        plan_tasks = []
        for label, job in jobs.items():
            probe_bins = job["bins"][1:]
            for k in range(0, len(probe_bins), args.chunk_bins):
                task_id = f"{label}_c{k // args.chunk_bins:03d}"
                plan_tasks.append(dict(job, task_id=task_id, bins=["all"] + probe_bins[k:k + args.chunk_bins],
                                       outputs={"json": str(fit_dir / "chunks" / f"{task_id}.json")}))
        write_new(plan_path, {"schema": "h4l_v3_task_plan/1", "stage": stage, "program": "fit_tnp", "created_utc": utc(),
                              "fit_run": str(fit_dir / "fit_run.json"), "tasks": plan_tasks})
        print(f"[{utc()}] Condor plan {plan_path}: {len(plan_tasks)} chunks; stage {stage} with analysis_v3/tnp/bin/fit_tnp")
        return 0
    if args.condor == "merge":
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
        chunks = {}
        for task in plan["tasks"]:
            out = Path(task["outputs"]["json"])
            if not out.exists():
                fail(f"missing chunk output {out}")
            report = json.loads(out.read_text(encoding="utf-8"))
            if report.get("task_id") != task["task_id"] or report.get("bins") != task["bins"] or \
                    report.get("data") != task["data"] or report.get("mc") != task["mc"]:
                fail(f"{out} does not belong to its task")
            if (report.get("frozen_program") or {}).get("sha256") != frozen_fit["programs"]["fit_tnp"]["sha256"]:
                fail(f"{out} comes from another fit_tnp than the bound one")
            chunks.setdefault(task["task_id"].rsplit("_c", 1)[0], []).append((task, report))
        records = []
        for label, job in jobs.items():
            parts = sorted(chunks.get(label, []), key=lambda item: item[0]["task_id"])
            covered = [b for task, _ in parts for b in task["bins"][1:]]
            if covered != job["bins"][1:]:
                fail(f"the chunks of {label} do not cover its bins exactly")
            results = []
            for role in ("data", "mc"):
                # The inclusive fit of the first chunk; every chunk's own probe bins.
                for index, (task, report) in enumerate(parts):
                    for entry in report["results"]:
                        if entry["role"] == role and (entry["bin"] != "all" or index == 0):
                            results.append(entry)
            merged = {"schema": "h4l_v3_tnp_fits/3", "flavour": job["flavour"], "step": job["step"], "data": job["data"],
                      "mc": job["mc"], "program": parts[0][1].get("frozen_program"), "bins": job["bins"],
                      "merged_from": [str(Path(task["outputs"]["json"])) for task, _ in parts], "results": results,
                      "finished_utc": utc()}
            out = fit_dir / "fits" / f"{label}.json"
            write_new(out, merged)
            records.append({"output": str(out)})
    else:
        records = parallel(tasks, args.workers)

    results = {}
    for rec in records:
        report = json.loads(Path(rec["output"]).read_text(encoding="utf-8"))
        for entry in report["results"]:
            results[(report["flavour"], report["step"], entry["bin"], entry["role"])] = entry

    def efficiency(entry: dict | None, model: str):
        if entry is None:
            return None, None, "missing"
        if entry.get("status") == "empty":
            return None, None, "empty"
        if entry.get("status") == "counting_only" or entry["fits"][model]["usable"] is False:
            c = entry["counting"]
            return c["efficiency"], c["error"], "counting"
        fit = entry["fits"][model]
        return fit["efficiency"], fit["efficiency_error"], entry["status"]

    # Alternatives that the gallery inspection found not to describe a bin (recorded with a note) do not enter
    # the fit-model systematic of that bin.
    excluded = {(item["flavour"], item["step"], item["bin"], item["model"]): item.get("note", "")
                for item in config["fit"].get("exclude_alternatives", [])}
    sf = {"schema": "h4l_v3_tnp_sf/3", "run": args.run, "fit_label": args.fit_label, "config": str(config_path),
          "calibration_payload": str(payload_path), "chain": list(CHAIN), "report_point": list(REPORT_POINT), "flavours": {}}
    for tag, name in FLAVOURS.items():
        bins_cfg = config["bins"][name]
        n_pt, n_eta = len(bins_cfg["pt_edges"]), len(bins_cfg["eta_edges"])
        flavour = {"bins": bins_cfg, "steps": {}}
        for step in STEPS:
            rows = []
            for b in ["all"] + [str(k) for k in range(n_pt * n_eta)]:
                row = {"bin": b}
                for role in ("data", "mc"):
                    entry = results.get((tag, step, b, role))
                    row[role] = {m: dict(zip(("efficiency", "error", "status"), efficiency(entry, m))) for m in MODELS}
                    row[f"{role}_counting"] = entry["counting"] if entry else None
                    if entry:
                        row[f"{role}_mean_pt"] = entry.get("mean_pt")
                        row[f"{role}_mean_abs_eta"] = entry.get("mean_abs_eta")
                nominal_d, nominal_m = row["data"]["nominal"], row["mc"]["nominal"]
                if nominal_d["efficiency"] and nominal_m["efficiency"] and nominal_d["error"] is not None \
                        and nominal_m["error"] is not None:
                    value = nominal_d["efficiency"] / nominal_m["efficiency"]
                    stat = value * math.hypot(nominal_d["error"] / nominal_d["efficiency"], nominal_m["error"] / nominal_m["efficiency"])
                    systematic = 0.0
                    used, skipped = [], {}
                    for model in ("alt_signal", "alt_background"):
                        d, m = row["data"][model], row["mc"][model]
                        if (tag, step, b, model) in excluded:
                            skipped[model] = excluded[(tag, step, b, model)]
                            continue
                        # An alternative whose efficiency error is inflated beyond max_error_ratio times the nominal one
                        # (a degenerate or failed fit) does not measure a model difference.
                        ratio_limit = config["fit"].get("alternatives", {}).get("max_error_ratio")
                        if ratio_limit and any(r[model]["error"] is None or r["nominal"]["error"] is None
                                               or r[model]["error"] > ratio_limit * max(r["nominal"]["error"], 1e-6)
                                               for r in (row["data"], row["mc"])):
                            skipped[model] = f"alternative error above {ratio_limit} x the nominal error"
                            continue
                        if d["efficiency"] and m["efficiency"] and d["status"] != "counting" and m["status"] != "counting":
                            systematic = max(systematic, abs(d["efficiency"] / m["efficiency"] - value))
                            used.append(model)
                    row["sf"] = {"value": value, "stat": stat, "fit_model": systematic, "alternatives_used": used,
                                 "alternatives_excluded": skipped}
                rows.append(row)
            flavour["steps"][step] = rows
        # The full single-lepton SF as the product of the chain (per bin).
        product = []
        for index in range(n_pt * n_eta + 1):
            factors = [flavour["steps"][s][index].get("sf") for s in CHAIN]
            if all(factors):
                value = math.prod(f["value"] for f in factors)
                stat = value * math.sqrt(sum((f["stat"] / f["value"]) ** 2 for f in factors))
                model = value * math.sqrt(sum((f["fit_model"] / f["value"]) ** 2 for f in factors))
                product.append({"bin": flavour["steps"]["id"][index]["bin"], "value": value, "stat": stat, "fit_model": model,
                                "total": math.hypot(stat, model)})
            else:
                product.append({"bin": flavour["steps"]["id"][index]["bin"], "value": None})
        flavour["product"] = product
        # The full SF at the report point: bilinear interpolation between the bin centres (the unbounded
        # last pT bin at the mean pT of its data probes).
        values, errors = {}, {}
        for index in range(n_pt * n_eta):
            ip, ie = divmod(index, n_eta)
            row = product[index + 1]
            values[(ip, ie)] = row.get("value")
            errors[(ip, ie)] = row.get("total", 0.0)
        last = [flavour["steps"]["full"][1 + (n_pt - 1) * n_eta + ie].get("data_mean_pt") for ie in range(n_eta)]
        last = [m for m in last if m]
        edges = bins_cfg["pt_edges"]
        pt_centres = [0.5 * (a + b) for a, b in zip(edges[:-1], edges[1:])] + [sum(last) / len(last) if last else 1.25 * edges[-1]]
        eta_top = 2.4 if tag == "mm" else 2.5
        eta_edges = bins_cfg["eta_edges"]
        eta_centres = [0.5 * (a + b) for a, b in zip(eta_edges, eta_edges[1:] + [eta_top])]
        flavour["report_point"] = {"pt": REPORT_POINT[0], "abs_eta": REPORT_POINT[1], "method": "bilinear between bin centres",
                                   "sf_full": interpolate(values, errors, pt_centres, eta_centres, *REPORT_POINT),
                                   "pt_centres": pt_centres, "eta_centres": eta_centres}
        sf["flavours"][name] = flavour
        rp = flavour["report_point"]["sf_full"]
        print(f"[{utc()}] {name}: SF(full) at {REPORT_POINT} = {rp.get('value')} +- {rp.get('error')}")
    write_new(fit_dir / "sf.json", sf)
    print(f"[{utc()}] scale factors: {fit_dir / 'sf.json'}")
    if closure:
        summary = {"schema": "h4l_v3_tnp_closure/1", "point": args.closure_point, "keep_id": point, "flavours": {}}
        for tag, name in FLAVOURS.items():
            flavour = sf["flavours"][name]
            k = point[name]
            per_bin = []
            for index, row in enumerate(flavour["product"]):
                direct = flavour["steps"]["full"][index].get("sf")
                entry = {"bin": row["bin"], "injected": k}
                if row.get("value") is not None:
                    entry["chain"] = {"value": row["value"], "error": row["total"], "pull": (row["value"] - k) / row["total"]
                                      if row["total"] > 0 else None}
                if direct:
                    error = math.hypot(direct["stat"], direct["fit_model"])
                    entry["direct"] = {"value": direct["value"], "error": error, "pull": (direct["value"] - k) / error if error > 0 else None}
                per_bin.append(entry)
            rp = flavour["report_point"]["sf_full"]
            summary["flavours"][name] = {"injected": k, "report_point_chain": rp,
                                         "report_point_pull": (rp["value"] - k) / rp["error"] if rp.get("value") and rp.get("error") else None,
                                         "bins": per_bin}
            print(f"[{utc()}] closure {name}: injected {k}, recovered at {REPORT_POINT} {rp.get('value')} +- {rp.get('error')}")
        write_new(fit_dir / "closure.json", summary)
        print(f"[{utc()}] closure summary: {fit_dir / 'closure.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
