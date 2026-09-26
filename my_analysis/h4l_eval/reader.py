"""The pass over the large input files runs in ROOT (C++, src/h4l_reader.cpp), compiled at run time against the ROOT of
the environment; independent chunks of files run in parallel processes.  Its compact outputs (control pairs, tag-and-
probe pairs, loose-lepton event records, Z + 1 lepton rows, file bookkeeping) are small ROOT files read here."""

from __future__ import annotations

import os
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

SRC = Path(__file__).resolve().parent.parent / "src" / "h4l_reader.cpp"


def compile_reader(build_dir: Path, log) -> Path:
    """Compile the ROOT reader with the environment's compiler against its ROOT (the C++ standard of root-config wins; the
    ROOT library directory as rpath so that the program runs without LD_LIBRARY_PATH)."""
    import shutil
    build_dir.mkdir(parents=True, exist_ok=True)
    exe = build_dir / "h4l_reader"
    candidates = [os.environ.get("CXX", ""), "x86_64-conda-linux-gnu-c++", "c++", "g++", "clang++"]
    cxx = next((c for c in candidates if c and shutil.which(c)), None)
    if cxx is None:
        raise RuntimeError("no C++ compiler found (CXX, c++, g++, clang++)")
    flags = subprocess.run(["root-config", "--cflags", "--libs"], check=True, capture_output=True, text=True).stdout.split()
    libdir = subprocess.run(["root-config", "--libdir"], check=True, capture_output=True, text=True).stdout.strip()
    cmd = [cxx, "-O2", "-std=c++17", str(SRC), "-o", str(exe)] + flags + [f"-Wl,-rpath,{libdir}"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"compilation failed: {shlex.join(cmd)}\n{res.stderr}")
    log(f"compiled {exe} ({cxx})")
    return exe


def chunks(tasks: list, mc_chunk: int = 40) -> list:
    """Data files one per chunk; MC files of one sample in chunks of mc_chunk files."""
    out, current = [], []
    for t in tasks:
        if not t["is_mc"]:
            out.append([t])
            continue
        if current and (current[0]["sample"] != t["sample"] or len(current) >= mc_chunk):
            out.append(current)
            current = []
        current.append(t)
    if current:
        out.append(current)
    return out


def run_chunks(exe: Path, work: Path, groups: list, workers: int, log) -> list:
    work.mkdir(parents=True, exist_ok=True)

    def run(i):
        lst = work / f"chunk_{i:05d}.txt"
        out = work / f"chunk_{i:05d}.root"
        with lst.open("w", encoding="utf-8") as stream:
            for t in groups[i]:
                stream.write(f"{t['path']}\t{'mc' if t['is_mc'] else 'data'}\t{t['key']}\t{t['want']}\n")
        res = subprocess.run([str(exe), str(out), str(lst)], capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"h4l_reader failed on chunk {i}: {res.stderr[-2000:]}")
        return read_chunk(out)

    outputs = [None] * len(groups)
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(run, i): i for i in range(len(groups))}
        for fut in futures:
            outputs[futures[fut]] = fut.result()
            done += 1
            if done % 50 == 0:
                log(f"  {done}/{len(groups)} chunks")
    return outputs


def read_chunk(path: Path) -> dict:
    out = {}
    with uproot.open(path) as f:
        files = f["Files"].arrays(library="np")
        out["files"] = {"path": [str(p) for p in files["path"]], "entries": files["entries"], "n_preselection": files["n_preselection"]}
        for tree in ("Calib", "TnP", "ZL"):
            t = f[tree]
            if t.num_entries:
                out[tree] = t.arrays(library="np")
        t = f["Rec"]
        if t.num_entries:
            a = t.arrays()
            rec = {"n_lep": ak.to_numpy(a["n_lep"]).astype(np.int64), "e_met": ak.to_numpy(a["e_met"]),
                   "e_entry": ak.to_numpy(a["e_entry"]), "e_file": ak.to_numpy(a["e_file"])}
            for name in a.fields:
                if name.startswith("l_"):
                    rec[name] = ak.to_numpy(ak.flatten(a[name]))
            out["Rec"] = rec
    return out


def collect(tasks_by_chunk: list, outputs: list) -> dict:
    """Merge the chunk outputs (already read) per sample in chunk order: {sample: {files, entries, n_preselection, outputs}}."""
    names = {"Calib": "calib", "TnP": "tnp", "Rec": "records", "ZL": "zl"}
    samples = {}
    for group, c in zip(tasks_by_chunk, outputs):
        sample = group[0]["sample"]
        s = samples.setdefault(sample, {"files": [], "entries": 0, "n_preselection": 0, "outputs": {}})
        s["files"] += c["files"]["path"]
        s["entries"] += int(np.sum(c["files"]["entries"]))
        s["n_preselection"] += int(np.sum(c["files"]["n_preselection"]))
        for tree, key in names.items():
            if tree in c:
                s["outputs"].setdefault(key, []).append(c[tree])
    for s in samples.values():
        merged = {}
        for key, parts in s["outputs"].items():
            fields = parts[0].keys()
            merged[key] = {k: np.concatenate([p[k] for p in parts]) for k in fields}
        s["outputs"] = merged
    return samples
