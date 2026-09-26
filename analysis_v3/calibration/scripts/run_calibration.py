"""Iterative lepton momentum-scale and resolution calibration (stage 3a),
factorized model (h4l/calibration.h, user decision 2026-09-24).

    pixi run py -- analysis_v3/calibration/scripts/run_calibration.py --version v3 --extract v2 --run nominal
    pixi run py -- analysis_v3/calibration/scripts/run_calibration.py --version v3 --extract v2 \\
        --run closure_null_halves --closure null_halves

Inputs: the extraction scan production_v3/calibration/<extract>/extract/scan.json
(outputs, row counts, sha256, per-sample normalization), the configuration
analysis_v3/calibration/config/calibration_ul16_v4.json (--config) and, for closure
runs, analysis_v3/calibration/config/closure_ul16_v1.json.  The programs
zpeak_histograms and template_fit are frozen by content under
production_v3/program_inputs/calibration_<version>/<sha12>/ and run from there.

Model per flavour: u = ln(1 + s) = a[e] + b_R(pT), v = r^2 = c[e] + d_R(pT),
e the fine |eta| bin, R the coarse region, b and d zero at the reference pT
bin.  Per iteration n, with the payload P_n (P_0: no correction):
 1. histograms, one job per extraction output: the data role (the pseudo-data
    corrected with P_n; delta-smeared md histograms) and the MC role (the MC
    smeared with P_n; TemplatePairs).  In closure runs the data role is MC
    with a known injection (halves: the even original files, the MC role the
    odd ones; same: all MC in both roles).
 2. at iteration 1: the MC role events play the data with no injection, a
    uniform scale and a uniform smear; their template fits give each
    category's response; categories with a response outside
    [min_response, max_response] are left out, the others enter with it.
 3. template fits of both category families (A: fine |eta| pairs with both
    legs above family_a_min_pt; B: (pT, region) pairs), in parallel chunks;
 4. one joint weighted least squares per flavour and quantity, e.g. for the
    scale of an A category (i, j):
       ln k = R [ (da_i + beta_i) + (da_j + beta_j) ] / 2,
       beta_s = sum_p f_sp db_{R(s), p}    (f: the pT composition of the legs)
    and of a B category (k, l):
       ln k = R [ (db_k + alpha_k) + (db_l + alpha_l) ] / 2,
       alpha_s = sum_e g_se da_e           (g: the |eta| composition);
    the smear likewise with E = [...] / 4.  Only categories whose data mode
    lies in the Z peak (z_mode_range) enter; outliers beyond outlier_pull are
    removed one by one (at most max_outlier_fraction) and recorded; a
    random-walk smoothness prior between neighbouring pT bins regularizes the
    pT terms the Z barely constrains; from iteration freeze_selection_after on
    the category set is frozen.  The covariance is inflated by chi2/ndf of the
    category residuals when that exceeds 1.  a, b, c, d are updated additively; the pT
    nodes are the mean pT of the corrected data legs in the Z peak;
 5. convergence: every parameter's residual below converged_sigma of its
    uncertainty.
Outputs in production_v3/calibration/<version>/<run>/: run.json, config.json,
iter_<n>/ (jobs, logs, hist, response_*, fits, response.json, solution.json,
payload.json), and after convergence payload.json and summary.json.  Nothing
is overwritten: an existing output is reused only if its recorded job equals
the new job.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import uproot

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
CONFIG = REPO / "analysis_v3/calibration/config/calibration_ul16_v4.json"
CLOSURE = REPO / "analysis_v3/calibration/config/closure_ul16_v1.json"
PROGRAMS = {name: REPO / "analysis_v3/calibration/bin" / name for name in ("zpeak_histograms", "template_fit")}
FLAVOURS = {"mm": "muon", "ee": "electron"}
FAMILIES = ("A", "B")


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def utc() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def write_new(path: Path, value) -> None:
    """Write JSON atomically; never overwrite (an identical existing file is accepted)."""
    text = json.dumps(value, indent=1, sort_keys=True) + "\n"
    if path.exists():
        if path.read_text(encoding="utf-8") != text:
            fail(f"{path} exists with a different content")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".partial.{os.getpid()}")
    temporary.write_text(text, encoding="utf-8")
    temporary.rename(path)


def freeze(version: str) -> dict:
    digests = {name: sha256_file(path) for name, path in PROGRAMS.items()}
    key = hashlib.sha256("".join(digests[name] for name in sorted(digests)).encode()).hexdigest()[:12]
    directory = PRODUCTION / "program_inputs" / f"calibration_{version}" / key
    directory.mkdir(parents=True, exist_ok=True)
    frozen = {}
    for name, source in PROGRAMS.items():
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


# ---------------------------------------------------------------------------
# Inputs and normalization.
def load_inputs(extract_version: str, manifests_version: str):
    scan_path = PRODUCTION / "calibration" / extract_version / "extract" / "scan.json"
    scan = json.loads(scan_path.read_text(encoding="utf-8"))
    manifests = PRODUCTION / "manifests" / manifests_version
    per_file = {}
    for name, sample in scan["samples"].items():
        if sample["kind"] == "mc":
            manifest = json.loads((manifests / f"mc_{name}.json").read_text(encoding="utf-8"))
            per_file[name] = {item["file_key"]: item["genEventSumw"] for item in manifest["files"]}
    return scan, per_file, scan_path


def mc_scale(scan: dict, per_file: dict, sample: str, half: int) -> float:
    """sigma L / sum genEventSumw over the files of the sample (of one key parity)."""
    info = scan["samples"][sample]
    lumi = scan["samples"]["pseudo_data"]["lumi_fb"]
    sumw = sum(w for key, w in per_file[sample].items() if half < 0 or int(key, 16) % 2 == half)
    if not sumw > 0:
        fail(f"{sample}: no genEventSumw for half {half}")
    return info["sigma_eff_pb"] * 1000.0 * lumi / sumw


# ---------------------------------------------------------------------------
# Jobs.
def run_job(program: dict, job: dict, job_path: Path, output: Path, log: Path) -> dict:
    """Run one frozen program on one job; reuse a finished identical output."""
    job = dict(job, program={"path": program["path"], "sha256": program["sha256"]})
    done = output.with_name(output.name + ".done.json")
    job_hash = hashlib.sha256(canonical(job).encode()).hexdigest()
    if done.exists():
        record = json.loads(done.read_text(encoding="utf-8"))
        if record["job_sha256"] != job_hash:
            fail(f"{output} was produced by another job")
        return record
    write_new(job_path, job)
    log.parent.mkdir(parents=True, exist_ok=True)
    if not output.exists():
        with log.open("a", encoding="utf-8") as stream:
            stream.write(f"# {utc()} {program['path']} --job {job_path} --out {output}\n")
            stream.flush()
            result = subprocess.run([program["path"], "--job", str(job_path), "--out", str(output)], stdout=stream,
                                    stderr=subprocess.STDOUT, env=dict(os.environ, PATH="/usr/bin:/bin"))
        if result.returncode != 0:
            fail(f"{Path(program['path']).name} failed ({result.returncode}) for {job_path}; see {log}")
    record = {"job": str(job_path), "job_sha256": job_hash, "output": str(output), "output_sha256": sha256_file(output),
              "finished_utc": utc()}
    write_new(done, record)
    return record


def parallel(tasks: list, workers: int) -> list:
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(*task) for task in tasks]
        return [future.result() for future in futures]


def histogram_jobs(ctx: dict, iteration_dir: Path, payload, label: str, specs: dict) -> dict:
    """Run the histogram jobs of one pass: specs maps the fit role ("data":
    delta-smeared md histograms; "mc": TemplatePairs) to its inputs, program
    role (data correction or MC smear), scales, half and injection."""
    cfg, scan = ctx["config"], ctx["scan"]
    tasks, files = [], {role: [] for role in specs}
    for output in scan["outputs"]:
        for role, spec in specs.items():
            if output["sample"] not in spec["samples"]:
                continue
            name = f"{role}_{output['task_id']}"
            job = {"calibration_config": cfg, "role": spec["role"], "weighted": output["kind"] == "mc",
                   "weight_scale": spec["samples"][output["sample"]], "half": spec.get("half", -1),
                   "delta_smear": role == "data", "write_template_pairs": role == "mc",
                   "write_full_pairs": spec.get("write_full_pairs", False),
                   "inputs": [{"path": output["path"], "calib_pairs": output["calib_pairs"], "root_sha256": output["root_sha256"]}]}
            if spec.get("inject"):
                job["inject"] = spec["inject"]
                job["inject_stream"] = 1
            if payload:
                job["payload"] = payload
            out = iteration_dir / label / "hist" / f"{name}.root"
            tasks.append((run_job, ctx["frozen"]["programs"]["zpeak_histograms"], job,
                          iteration_dir / label / "jobs" / f"{name}.json", out, iteration_dir / label / "logs" / f"{name}.log"))
            files[role].append(str(out))
    parallel(tasks, ctx["workers"])
    return files


def categories(files: list, family: str, tag: str) -> list:
    prefix = f"md{family}_{tag}_"
    names = set()
    for path in files:
        with uproot.open(path) as f:
            for key in f.keys(recursive=False):
                name = key.split(";")[0]
                if name.startswith(prefix):
                    names.add(name[len(prefix):])
    return sorted(names, key=lambda c: tuple(int(x) for x in c.split("_")))


def template_fits(ctx: dict, iteration_dir: Path, label: str, files: dict, chunks: int, curves: bool = True) -> dict:
    """Template fits of every category of both families; {family: {tag: [results]}}."""
    tasks = []
    for family in FAMILIES:
        for tag in FLAVOURS:
            cats = categories(files["data"], family, tag)
            parts = chunks if family == "B" else max(1, chunks // 2)
            for part in range(parts):
                subset = cats[part::parts]
                if not subset:
                    continue
                job = {"calibration_config": ctx["config"], "family": family, "flavour": tag, "data": files["data"],
                       "mc": files["mc"], "categories": subset, "curves": curves}
                name = f"tfit_{family}_{tag}_{part:02d}"
                tasks.append((run_job, ctx["frozen"]["programs"]["template_fit"], job,
                              iteration_dir / label / "jobs" / f"{name}.json", iteration_dir / label / "fits" / f"{name}.json",
                              iteration_dir / label / "logs" / f"{name}.log"))
    results = {family: {tag: [] for tag in FLAVOURS} for family in FAMILIES}
    for record in parallel(tasks, ctx["workers"]):
        report = json.loads(Path(record["output"]).read_text(encoding="utf-8"))
        results[report["family"]][report["flavour"]].extend(report["results"])
    return results


# ---------------------------------------------------------------------------
# Model bookkeeping.
class Model:
    """Binning of one flavour and the parameter indexing of the joint least squares."""

    def __init__(self, node: dict):
        self.node = node
        self.eta = node["eta_edges"]
        self.regions = node["region_edges"]
        self.pt = node["pt_edges"]
        self.ref = node["reference_pt_bin"]
        self.n_eta, self.n_regions, self.n_pt = len(self.eta), len(self.regions), len(self.pt)
        for edge in self.regions:
            if edge not in self.eta:
                fail("every region edge must be a fine |eta| edge")
        self.region_of_eta = [max(r for r, edge in enumerate(self.regions) if edge <= self.eta[e]) for e in range(self.n_eta)]
        self.b_index = {}
        for r in range(self.n_regions):
            for p in range(self.n_pt):
                if p != self.ref:
                    self.b_index[(r, p)] = self.n_eta + len(self.b_index)
        self.n_par = self.n_eta + len(self.b_index)

    def q_bin(self, q: int) -> tuple:
        """(pT bin, region) of a family-B bin index q = p * n_regions + r."""
        return divmod(q, self.n_regions)


def compositions(files: list, family: str, tag: str) -> dict:
    """Summed leg compositions per category of one family."""
    out = {}
    prefix = f"w{family}_{tag}_"
    for path in files:
        with uproot.open(path) as f:
            for key in f.keys(recursive=False):
                name = key.split(";")[0]
                if name.startswith(prefix):
                    values = f[key].values()
                    c = name[len(prefix):]
                    out[c] = out.get(c, 0) + values
    return out


def design_row(model: Model, family: str, category: str, composition, weight: float) -> np.ndarray:
    """The coefficients of one category's pair-level residual in the parameters."""
    row = np.zeros(model.n_par)
    i, j = (int(x) for x in category.split("_"))
    size = model.n_pt if family == "A" else model.n_eta
    sets = [composition[:size], composition[size:]]
    legs = [(i, sets[0])] if i == j else [(i, sets[0]), (j, sets[1])]
    per_leg = 2 * weight if i == j else weight  # both legs of a diagonal category share one bin
    for bin_index, comp in legs:
        total = comp.sum()
        if not total > 0:
            return None
        fractions = comp / total
        if family == "A":
            row[bin_index] += per_leg
            region = model.region_of_eta[bin_index]
            for p in range(model.n_pt):
                if p != model.ref and fractions[p] > 0:
                    row[model.b_index[(region, p)]] += per_leg * fractions[p]
        else:
            p, r = model.q_bin(bin_index)
            if p != model.ref:
                row[model.b_index[(r, p)]] += per_leg
            for e in range(model.n_eta):
                if fractions[e] > 0:
                    row[e] += per_leg * fractions[e]
    return row


def prior_rows(model: Model, current: list, tau: float) -> tuple:
    """Random-walk smoothness rows between neighbouring pT bins of each region,
    on the total terms (current value + residual): (t_{p+1} - t_p) = 0 +- tau."""
    rows, y = [], []
    for r in range(model.n_regions):
        for p in range(model.n_pt - 1):
            row = np.zeros(model.n_par)
            for q, sign in ((p + 1, 1.0), (p, -1.0)):
                if q != model.ref:
                    row[model.b_index[(r, q)]] += sign
            rows.append(row / tau)
            y.append(-(current[r][p + 1] - current[r][p]) / tau)
    return rows, y


def solve(model: Model, fits: dict, comps: dict, response: dict | None, limits: tuple, selection: dict,
          current: dict, allowed: dict | None) -> dict:
    """Joint least squares of both families for the scale (ln k) and the smear (E),
    with the Z-mode category selection, the iterative outlier rejection (unless
    the category set is frozen: allowed) and the pT smoothness prior."""
    out = {}
    z_low, z_high = selection["z_mode_range"]
    for quantity, key, weight in (("scale", "lnk", 0.5), ("smear", "E", 0.25)):
        rows, y, sigma, used, excluded = [], [], [], [], []
        for family in FAMILIES:
            for entry in fits[family]:
                if not entry.get("ok"):
                    continue
                label = f"{family}_{entry['category']}"
                if allowed is not None and label not in allowed[quantity]:
                    excluded.append(label)
                    continue
                if not (z_low <= entry["data_mode"] <= z_high):
                    excluded.append(label)
                    continue
                # A fit whose free normalization misses the data sum in its window did not describe the data.
                if abs(entry["fit"]["norm"]["value"] / entry["data_entries"] - 1.0) > selection["max_norm_deviation"]:
                    excluded.append(label)
                    continue
                r = 1.0
                if response is not None:
                    r = response.get(label, {}).get(quantity)
                    if r is None or not (limits[0] <= r <= limits[1]):
                        excluded.append(label)
                        continue
                comp = comps[family].get(entry["category"])
                if comp is None:
                    excluded.append(label)
                    continue
                row = design_row(model, family, entry["category"], comp, weight)
                if row is None:
                    excluded.append(label)
                    continue
                fit = entry["fit"]
                if not fit[key]["error"] > 0:
                    fail(f"{label}: a category fit has no uncertainty")
                rows.append(r * row)
                y.append(fit[key]["value"])
                sigma.append(fit[key]["error"])
                used.append(label)
        a = np.array(rows)
        y = np.array(y)
        sigma = np.array(sigma)
        terms = current["b" if quantity == "scale" else "d"]
        p_rows, p_y = prior_rows(model, terms, selection["pt_smoothness"][quantity])
        p_rows, p_y = np.array(p_rows), np.array(p_y)
        keep = np.ones(len(rows), dtype=bool)
        outliers = []
        while True:
            data_rows = a[keep] / sigma[keep][:, None]
            design = np.vstack([data_rows, p_rows])
            target = np.concatenate([y[keep] / sigma[keep], p_y])
            determined = np.flatnonzero(np.abs(data_rows).sum(axis=0) > 0)
            ad = design[:, determined]
            normal = ad.T @ ad
            if np.linalg.matrix_rank(normal) < len(determined):
                fail(f"{quantity}: the joint system is singular")
            cov = np.linalg.inv(normal)
            solution = cov @ ad.T @ target
            residual = target[: int(keep.sum())] - data_rows[:, determined] @ solution
            worst = int(np.argmax(np.abs(residual)))
            if allowed is not None or abs(residual[worst]) <= selection["outlier_pull"] or \
                    len(outliers) + 1 > selection["max_outlier_fraction"] * len(rows):
                break
            index = np.flatnonzero(keep)[worst]
            outliers.append({"category": used[index], "pull": float(residual[worst])})
            keep[index] = False
        used = [c for c, k in zip(used, keep) if k]
        prior_residual = p_y - p_rows[:, determined] @ solution
        result = {"value": [0.0] * model.n_par, "error": [math.nan] * model.n_par, "determined": [False] * model.n_par,
                  "categories": int(keep.sum()), "excluded": excluded, "outliers": outliers,
                  "prior_chi2": float(prior_residual @ prior_residual)}
        chi2 = float(residual @ residual)
        ndf = int(keep.sum()) - len(determined)
        inflation = max(1.0, chi2 / ndf) if ndf > 0 else 1.0
        cov = cov * inflation
        for position, k in enumerate(determined):
            result["value"][k] = float(solution[position])
            result["error"][k] = float(math.sqrt(cov[position, position]))
            result["determined"][k] = True
        result.update({"chi2": chi2, "ndf": ndf, "inflation": inflation,
                       "covariance": {"parameters": [int(k) for k in determined], "matrix": cov.tolist()},
                       "pulls": {c: float(p) for c, p in zip(used, residual)}})
        out[quantity] = result
    return out


def pt_nodes(files: list, tag: str, model: Model, previous: list | None) -> list:
    """Mean pT of the data legs in the Z peak per (region, pT bin)."""
    w = wpt = None
    for path in files:
        with uproot.open(path) as f:
            a, b = f[f"leg_{tag}_w"].values(), f[f"leg_{tag}_wpt"].values()
        w, wpt = (a, b) if w is None else (w + a, wpt + b)
    nodes = []
    for r in range(model.n_regions):
        row = []
        for p in range(model.n_pt):
            q = p * model.n_regions + r
            row.append(float(wpt[q] / w[q]) if w[q] > 0 else math.nan)
        if any(math.isnan(x) for x in row):
            if previous is None:
                fail(f"{tag}: a (region, pT) bin has no leg in the Z peak")
            row = [p if math.isnan(x) else x for x, p in zip(row, previous[r])]
        if any(b <= a for a, b in zip(row, row[1:])):
            fail(f"{tag}: pT nodes are not increasing: {row}")
        nodes.append(row)
    return nodes


def update_payload(model: Model, previous: dict | None, solution: dict, nodes: list) -> dict:
    """a, b, c, d updated additively by the solved residuals."""
    zero_rows = [[0.0] * model.n_pt for _ in range(model.n_regions)]
    old = previous or {"a": [0.0] * model.n_eta, "c": [0.0] * model.n_eta, "b": zero_rows, "d": zero_rows}
    new = {"model": model.node, "node_pt": nodes}
    for eta_key, pt_key, quantity in (("a", "b", "scale"), ("c", "d", "smear")):
        sol = solution[quantity]
        new[eta_key] = [old[eta_key][e] + (sol["value"][e] if sol["determined"][e] else 0.0) for e in range(model.n_eta)]
        rows = []
        for r in range(model.n_regions):
            row = []
            for p in range(model.n_pt):
                if p == model.ref:
                    row.append(0.0)
                    continue
                k = model.b_index[(r, p)]
                row.append(old[pt_key][r][p] + (sol["value"][k] if sol["determined"][k] else 0.0))
            rows.append(row)
        new[pt_key] = rows
    return new


def response_values(base: dict, variants: dict, injected: dict) -> dict:
    """Per category, (fit(variant) - fit(baseline)) / injected per-pair value."""
    out = {}
    for tag in FLAVOURS:
        values = {}
        for family in FAMILIES:
            base_fits = {e["category"]: e for e in base[family][tag] if e.get("ok")}
            for quantity, key in (("scale", "lnk"), ("smear", "E")):
                for entry in variants[quantity][family][tag]:
                    c = entry["category"]
                    if not entry.get("ok") or c not in base_fits:
                        continue
                    delta = entry["fit"][key]["value"] - base_fits[c]["fit"][key]["value"]
                    values.setdefault(f"{family}_{c}", {})[quantity] = delta / injected[quantity]
        out[tag] = values
    return out


# ---------------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="calibration product version (the configuration version)")
    parser.add_argument("--extract", required=True, help="version of the extraction whose scan is the input")
    parser.add_argument("--run", required=True)
    parser.add_argument("--manifests", default="v2")
    parser.add_argument("--closure", help="closure scenario of the closure configuration")
    parser.add_argument("--closure-config", type=Path, default=CLOSURE, help="closure configuration (default %(default)s)")
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--chunks", type=int, default=6, help="template-fit processes per flavour of family B")
    parser.add_argument("--config", type=Path, default=CONFIG, help="calibration configuration (default %(default)s)")
    args = parser.parse_args()
    base = PRODUCTION / "calibration" / args.version / args.run
    config_path = args.config if args.config.is_absolute() else REPO / args.config
    config = json.loads(config_path.read_text(encoding="utf-8"))
    closure_path = args.closure_config if args.closure_config.is_absolute() else REPO / args.closure_config
    closure_cfg = json.loads(closure_path.read_text(encoding="utf-8"))
    scan, per_file, scan_path = load_inputs(args.extract, args.manifests)
    if canonical(config["extract"]) != canonical(scan["extract_config"]):
        fail("the extraction section of the configuration differs from the scanned extraction")
    frozen = freeze(args.version)
    run_record = {"version": args.version, "run": args.run, "config": str(config_path), "config_sha256": sha256_file(config_path),
                  "closure": args.closure, "closure_config_sha256": sha256_file(closure_path) if args.closure else None,
                  "extract_scan": str(scan_path), "extract_scan_sha256": sha256_file(scan_path),
                  "extract_config_fnv1a64": scan["extract_config_fnv1a64"], "frozen": frozen, "manifests": args.manifests,
                  "orchestrator_sha256": sha256_file(Path(__file__))}
    run_path = base / "run.json"
    if run_path.exists():
        existing = json.loads(run_path.read_text(encoding="utf-8"))
        if canonical({k: v for k, v in existing.items() if k != "started_utc"}) != canonical(run_record):
            fail(f"{run_path} records another configuration or program; use a new --run name")
    else:
        write_new(run_path, dict(run_record, started_utc=utc()))
    write_new(base / "config.json", config)
    ctx = {"config": config, "scan": scan, "frozen": frozen, "workers": args.workers}

    mc_samples = [name for name, s in scan["samples"].items() if s["kind"] == "mc"]
    if args.closure:
        scenario = closure_cfg["scenarios"][args.closure]
        halves = scenario["mode"] == "halves"
        closure_samples = mc_samples
        if halves:
            # A sample whose files do not populate both key parities cannot be split into independent halves;
            # it is left out of both roles (only single-file samples, negligible in the Z peak).
            closure_samples = [s for s in mc_samples
                               if all(any(int(k, 16) % 2 == h for k in per_file[s]) for h in (0, 1))]
            left_out = sorted(set(mc_samples) - set(closure_samples))
            if left_out:
                print(f"[{utc()}] halves mode: samples without files in both halves left out: {left_out}", flush=True)
        data_spec = {"role": "data", "samples": {s: mc_scale(scan, per_file, s, 0 if halves else -1) for s in closure_samples},
                     "half": 0 if halves else -1, "inject": scenario["inject"]}
        mc_spec = {"role": "mc", "samples": {s: mc_scale(scan, per_file, s, 1 if halves else -1) for s in closure_samples},
                   "half": 1 if halves else -1}
    else:
        data_spec = {"role": "data", "samples": {"pseudo_data": 1.0}}
        mc_spec = {"role": "mc", "samples": {s: mc_scale(scan, per_file, s, -1) for s in mc_samples}}
    response_cfg = config["response"]
    # Every response variant carries the same decorrelating smear r0: a baseline whose data are the template
    # events themselves pulls D to its lower limit (the template reproduces the events' own fluctuations best
    # when it is sharpest), and a baseline at the limit voids both responses of the category.
    r0 = response_cfg["decorrelation_smear"]
    response_injections = {
        "baseline": {name: {"scale": {"a": 0.0}, "smear": {"a": r0}} for name in FLAVOURS.values()},
        "scale": {name: {"scale": {"a": response_cfg["scale"]}, "smear": {"a": r0}} for name in FLAVOURS.values()},
        "smear": {name: {"scale": {"a": 0.0}, "smear": {"a": math.hypot(r0, response_cfg["smear"])}}
                  for name in FLAVOURS.values()}}
    injected = {"scale": math.log1p(response_cfg["scale"]), "smear": 2 * response_cfg["smear"] ** 2 / 4}
    limits = (response_cfg["min_response"], response_cfg["max_response"])

    models = {tag: Model(config["model"][name]) for tag, name in FLAVOURS.items()}
    max_iterations = config["iteration"]["max_iterations"]
    converged_sigma = config["iteration"]["converged_sigma"]
    payload = None
    response = None
    frozen_sets = {}
    history = []
    for iteration in range(max_iterations + 1):
        idir = base / f"iter_{iteration:02d}"
        print(f"[{utc()}] iteration {iteration}: histograms", flush=True)
        files = histogram_jobs(ctx, idir, payload, "main", {"data": data_spec, "mc": mc_spec})
        if iteration == 1:
            print(f"[{utc()}] iteration {iteration}: response passes", flush=True)
            fits = {}
            for variant, injection in response_injections.items():
                spec = {"role": "mc", "samples": mc_spec["samples"], "half": mc_spec.get("half", -1)}
                if injection:
                    spec["inject"] = injection
                variant_files = histogram_jobs(ctx, idir, payload, f"response_{variant}", {"data": spec})
                variant_files["mc"] = files["mc"]
                fits[variant] = template_fits(ctx, idir, f"response_{variant}", variant_files, args.chunks, curves=False)
            response = response_values(fits["baseline"], {"scale": fits["scale"], "smear": fits["smear"]}, injected)
            write_new(idir / "response.json", {"injected_per_pair": injected, "limits": limits, "response": response})
        print(f"[{utc()}] iteration {iteration}: template fits", flush=True)
        fits = template_fits(ctx, idir, "main", files, args.chunks)
        solution, new_payload, worst = {}, {}, 0.0
        for tag, name in FLAVOURS.items():
            model = models[tag]
            comps = {family: compositions(files["data"], family, tag) for family in FAMILIES}
            flavour_fits = {family: fits[family][tag] for family in FAMILIES}
            zero = [[0.0] * model.n_pt for _ in range(model.n_regions)]
            current = {"b": payload[name]["b"] if payload else zero, "d": payload[name]["d"] if payload else zero}
            solution[tag] = solve(model, flavour_fits, comps, response[tag] if response else None, limits,
                                  {**config["categories"], **config["least_squares"]}, current, frozen_sets.get(tag))
            nodes = pt_nodes(files["data"], tag, model, payload[name]["node_pt"] if payload else None)
            new_payload[name] = update_payload(model, payload[name] if payload else None, solution[tag], nodes)
            for quantity in ("scale", "smear"):
                sol = solution[tag][quantity]
                for k in range(model.n_par):
                    if sol["determined"][k] and sol["error"][k] > 0:
                        worst = max(worst, abs(sol["value"][k]) / sol["error"][k])
        if iteration == config["least_squares"]["freeze_selection_after"]:
            frozen_sets = {tag: {q: set(solution[tag][q]["pulls"]) for q in ("scale", "smear")} for tag in FLAVOURS}
            write_new(idir / "frozen_selection.json", {tag: {q: sorted(v) for q, v in sets.items()}
                                                       for tag, sets in frozen_sets.items()})
        write_new(idir / "solution.json", solution)
        write_new(idir / "payload.json", new_payload)
        history.append({"iteration": iteration, "worst_residual_sigma": worst,
                        "chi2_ndf": {tag: {q: [solution[tag][q]["chi2"], solution[tag][q]["ndf"]] for q in ("scale", "smear")}
                                     for tag in FLAVOURS}})
        print(f"[{utc()}] iteration {iteration}: worst residual {worst:.2f} sigma; chi2/ndf "
              + ", ".join(f"{tag} {q} {solution[tag][q]['chi2']:.0f}/{solution[tag][q]['ndf']}"
                          for tag in FLAVOURS for q in ("scale", "smear")), flush=True)
        if iteration >= 1 and worst < converged_sigma:
            write_new(base / "payload.json", {"schema": "h4l_v3_lepton_calibration_payload/2", "run": args.run,
                                              "converged_iteration": iteration,
                                              "applied_payload": str(base / f"iter_{iteration - 1:02d}" / "payload.json"),
                                              "solution": str(idir / "solution.json"), "flavours": payload})
            break
        payload = new_payload
    else:
        write_new(base / "summary.json", {"iterations": history, "converged": False, "finished_utc": utc()})
        fail(f"no convergence after {max_iterations} iterations; see {base}")
    write_new(base / "summary.json", {"iterations": history, "converged": True, "converged_iteration": history[-1]["iteration"],
                                      "payload": str(base / "payload.json"), "finished_utc": utc()})
    print(f"[{utc()}] converged; {base / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
