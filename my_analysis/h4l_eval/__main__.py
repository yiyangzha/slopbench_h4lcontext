"""Evaluation run: python -m h4l_eval <dataset_dir> <output_dir>  (called by run.sh; see README.md)."""

from __future__ import annotations

import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

from . import analysis as A
from . import calibration as CAL
from . import config as C
from . import dataset as D
from . import model as M
from . import reader as R
from . import selection as S
from . import tnp as T
from . import zx as Z
from .outputs import write_outputs

MC_EVERY = {"DY": 2, "qqZZ": 2, "TT": 2}  # every n-th file (sorted names) of the large samples: a deterministic half
N_TOYS = 200
SEED = 20260926


class Log:
    def __init__(self, path: Path):
        self.stream = path.open("a", encoding="utf-8")
        self.t0 = time.time()

    def __call__(self, msg: str):
        line = f"[{time.time() - self.t0:7.1f}s] {msg}"
        print(line, flush=True)
        self.stream.write(line + "\n")
        self.stream.flush()


def plan(ds: dict) -> list:
    tasks = []
    for path in ds["data_files"]:
        tasks.append({"sample": "data", "role": "data", "path": path, "is_mc": False, "key": D.file_key(Path(path).name),
                      "want": "calib,tnp,rec_data"})
    for name, proc in ds["processes"].items():
        role = proc["role"]
        if role == "other":
            continue
        every = MC_EVERY.get(role, 1)
        want = {"DY": "calib,tnp", "TT": "calib", "qqZZ": "rec_mc,zl", "ggZZ": "rec_mc,zl", "signal": "rec_mc"}[role]
        for i, path in enumerate(proc["files"]):
            if i % every:
                continue
            tasks.append({"sample": name, "role": role, "path": path, "is_mc": True, "key": D.file_key(f"{name}/{Path(path).name}"),
                          "want": want})
    return tasks


def main(argv) -> int:
    if len(argv) != 2:
        print("usage: python -m h4l_eval <dataset_dir> <output_dir>", file=sys.stderr)
        return 2
    dataset_dir, out_dir = Path(argv[0]), Path(argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "logs").mkdir(exist_ok=True)
    log = Log(out_dir / "logs" / "run.log")
    timings = {}
    t_start = time.time()
    ds = D.discover(dataset_dir)
    log(f"dataset {dataset_dir}: lumi {ds['lumi_fb']} fb^-1, {len(ds['data_files'])} data files, processes "
        + ", ".join(f"{n} ({p['role']}, {len(p['files'])} files, xsec {p['xsec_pb']})" for n, p in ds["processes"].items()))
    tasks = plan(ds)
    workers = max(1, min(os.cpu_count() or 1, int(os.environ.get("H4L_WORKERS", "16"))))
    debug_load = os.environ.get("H4L_DEBUG_LOAD")
    if debug_load:
        # Development only: the merged compact outputs of an earlier run of the same dataset.
        import pickle
        with open(debug_load, "rb") as stream:
            samples = pickle.load(stream)
        log(f"loaded the compact outputs from {debug_load}")
    else:
        exe = R.compile_reader(out_dir / "build", log)
        groups = R.chunks(tasks)
        log(f"reading {len(tasks)} files in {len(groups)} chunks with {workers} parallel ROOT processes")
        outputs = R.run_chunks(exe, out_dir / "work", groups, workers, log)
        samples = R.collect(groups, outputs)
        debug_dump = os.environ.get("H4L_DEBUG_DUMP")
        if debug_dump:
            import pickle
            with open(debug_dump, "wb") as stream:
                pickle.dump(samples, stream, protocol=pickle.HIGHEST_PROTOCOL)
    timings["read"] = time.time() - t_start
    log(f"read done in {timings['read']:.0f} s")

    # ------------------------------------------------------------------ normalization (count based, as prescribed)
    norm = {}
    for name, proc in ds["processes"].items():
        if name not in samples:
            continue
        s = samples[name]
        if proc["xsec_pb"] is None:
            raise SystemExit(f"no cross section for {name} in mc/cross_sections.json")
        n_pre = s["n_preselection"]
        norm[name] = {"xsec_pb": proc["xsec_pb"], "n_preselection_processed": n_pre, "files_processed": len(s["files"]),
                      "files_total": len(proc["files"]), "entries_processed": s["entries"],
                      "weight": proc["xsec_pb"] * ds["lumi_fb"] * 1000.0 / n_pre, "meta_count": D.meta_count(proc["meta"])}
        log(f"  {name}: weight {norm[name]['weight']:.4g} per entry ({len(s['files'])}/{len(proc['files'])} files, "
            f"n_preselection {n_pre}; meta.json count {norm[name]['meta_count']})")

    # ------------------------------------------------------------------ lepton calibration
    t = time.time()
    data_calib = samples["data"]["outputs"]["calib"]
    dy_names = [n for n, p in ds["processes"].items() if p["role"] == "DY" and n in samples]
    tmpl_names = [n for n, p in ds["processes"].items() if p["role"] in ("DY", "TT") and n in samples and "calib" in samples[n]["outputs"]]
    mc_calib = {k: np.concatenate([samples[n]["outputs"]["calib"][k] for n in tmpl_names]) for k in data_calib}
    mc_calib["w"] = np.concatenate([np.full(len(samples[n]["outputs"]["calib"]["mass"]), norm[n]["weight"]) for n in tmpl_names])
    mc_calib["is_dy"] = np.concatenate([np.full(len(samples[n]["outputs"]["calib"]["mass"]), ds["processes"][n]["role"] == "DY")
                                        for n in tmpl_names])
    cal = CAL.calibrate(data_calib, mc_calib, log)
    models = cal["models"]
    cal_summary = CAL.data_weighted_summary(models, data_calib)
    timings["calibration"] = time.time() - t
    for fl, v in cal_summary.items():
        log(f"[calib] {fl}: scale_shift {v['scale_shift']['value']:+.5f} +- {v['scale_shift']['unc']:.5f}, smear "
            f"{v['smear']['value']:.5f} +- {v['smear']['unc']:.5f}")

    # ------------------------------------------------------------------ tag and probe
    t = time.time()
    data_tnp = samples["data"]["outputs"]["tnp"]
    mc_tnp = {k: np.concatenate([samples[n]["outputs"]["tnp"][k] for n in dy_names]) for k in data_tnp}
    tnp = T.measure(data_tnp, mc_tnp, models, log)
    sfs = T.ScaleFactors(tnp)
    sel_eff = T.data_weighted_sel_eff(tnp, data_calib, models)
    timings["tnp"] = time.time() - t
    for fl, v in sel_eff.items():
        log(f"[tnp] {fl}: sel_eff {v['value']:.5f} +- {v['unc']:.5f}")

    # ------------------------------------------------------------------ event selection
    t = time.time()
    lo, hi = C.WINDOW
    sel_d = A.select_events(samples["data"]["outputs"]["records"], models, False, control=True, zl=True)
    data_sr = A.with_kinematics(sel_d, "sr")
    cr3 = A.with_kinematics(sel_d, "cr3")
    cr2 = A.with_kinematics(sel_d, "cr2")
    mc_sr, mc_cr3 = {}, {}
    for name, proc in ds["processes"].items():
        if name not in samples or proc["role"] not in ("signal", "qqZZ", "ggZZ"):
            continue
        rec = samples[name]["outputs"].get("records")
        if rec is None:
            continue
        sel = A.select_events(rec, models, True, control=proc["role"] in ("qqZZ", "ggZZ"), zl=False)
        c = A.with_kinematics(sel, "sr")
        info = A.event_sf(sel, c, sfs)
        c["w"] = norm[name]["weight"] * info["sf"]
        c["sfinfo"] = info
        c["sel"] = sel
        mc_sr[name] = c
        if proc["role"] in ("qqZZ", "ggZZ"):
            c3 = A.with_kinematics(sel, "cr3")
            c3["w"] = norm[name]["weight"] * np.ones(len(c3["event"]))
            mc_cr3[name] = c3
    timings["selection"] = time.time() - t
    in_win = (data_sr["m4l_refit"] > lo) & (data_sr["m4l_refit"] < hi)
    fsd = data_sr["final_state"]
    log(f"[select] data: {len(data_sr['event'])} SR candidates, {int(in_win.sum())} in {lo}-{hi} GeV ("
        + ", ".join(f"{n} {int(np.sum(in_win & (fsd == k)))}" for k, n in enumerate(A.FS_NAMES)) + ")")

    # ------------------------------------------------------------------ Z + X
    t = time.time()
    zl = sel_d["zl"]
    lep_d, pt_d = sel_d["lep"], sel_d["pt"]
    pr = zl["probe"].astype(np.int64)
    data_zl = {"flavour": lep_d["flavour"][pr], "pt": pt_d[pr], "abs_eta": A.abs_eta_of(lep_d)[pr], "pass": zl["pass"].astype(bool),
               "sip_ok": zl["sip_ok"].astype(bool)}
    mc_zl = []
    for name, proc in ds["processes"].items():
        if name in samples and proc["role"] in ("qqZZ", "ggZZ") and "zl" in samples[name]["outputs"]:
            rows = samples[name]["outputs"]["zl"]
            mc_zl.append((rows, np.full(len(rows["pt"]), norm[name]["weight"])))
    rates = Z.fake_rates(data_zl, mc_zl)

    def cr_rows(c, n):
        win = (c["m4l_refit"] > lo) & (c["m4l_refit"] < hi)
        sub = {k: v[win] for k, v in c.items() if isinstance(v, np.ndarray) and v.shape[:1] == win.shape}
        fl, fp, fe = A.failing_legs(sel_d, sub, n)
        return {"final_state": sub["final_state"], "fail_flavour": fl, "fail_pt": fp, "fail_eta": fe}

    rows3 = cr_rows(cr3, 1)
    rows2 = cr_rows(cr2, 2)
    zz3 = []
    for name, c3 in mc_cr3.items():
        win = (c3["m4l_refit"] > lo) & (c3["m4l_refit"] < hi)
        zz3.append(({"final_state": c3["final_state"][win]}, c3["w"][win]))
    zx_est = Z.os_estimate(rows3, rows2, zz3, rates)
    timings["zx"] = time.time() - t
    for fs, v in zx_est.items():
        log(f"[zx] {fs}: {v['yield_window']:.3f} +- {v['stat']:.3f} (3P1F {v['n_3p1f']:.0f}, 2P2F {v['n_2p2f']:.0f}, ZZ in 3P1F "
            f"{v['n_zz_3p1f']:.2f})")

    # ------------------------------------------------------------------ templates and likelihood
    t = time.time()
    from .templates import build_channels
    channels, observed, meta = build_channels(ds, samples, mc_sr, data_sr, zx_est, cr3, cr2, models, cal_summary, sfs, norm, log)
    lik = M.Likelihood(channels, observed)
    fitter = M.Fitter(lik)
    timings["templates"] = time.time() - t

    t = time.time()
    from .inference import run_inference
    inference = run_inference(lik, fitter, log, n_toys=N_TOYS, seed=SEED)
    timings["inference"] = time.time() - t
    timings["total"] = time.time() - t_start
    write_outputs(out_dir, ds, lik, inference, cal, cal_summary, tnp, sel_eff, zx_est, rates, meta, norm, timings, log)
    log(f"done in {timings['total']:.0f} s")
    return 0


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    raise SystemExit(main(sys.argv[1:]))
