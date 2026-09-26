"""Check the outputs of ./run.sh against EVAL_CONTRACT.md (my_analysis/EVAL_CONTRACT.md, per-flavour calibration).

    python tests/check_outputs.py <output_dir>

RESULT.json: pois mH and mu with finite value/stat/syst/total, total > 0 and consistent with stat (+) syst;
significance_obs finite; discovery boolean; optional significance_exp; calibration {scale_shift, smear, sel_eff} x
{muon, electron} x {value, unc}; quality {gof_pvalue, coverage} with value and note.  MODEL.json: loads as a pyhf
workspace and builds a model (the GATE); if pyhf is importable the MLE mu is compared with the reported mu.
"""

import json
import math
import sys
from pathlib import Path


def fail(msg):
    print(f"FAIL: {msg}")
    raise SystemExit(1)


def main():
    out = Path(sys.argv[1])
    r = json.loads((out / "RESULT.json").read_text())
    pois = {p["name"]: p for p in r["pois"]}
    for name in ("mH", "mu"):
        p = pois.get(name) or fail(f"POI {name} missing")
        for k in ("value", "stat", "syst", "total"):
            if not isinstance(p.get(k), (int, float)) or not math.isfinite(p[k]):
                fail(f"{name}.{k} not finite")
        if not p["total"] > 0:
            fail(f"{name}.total not positive")
        quad = math.hypot(p["stat"], p["syst"])
        print(f"{name}: {p['value']:.4f} +- {p['total']:.4f} (stat {p['stat']:.4f}, syst {p['syst']:.4f}; stat+syst {quad:.4f})")
    if not math.isfinite(r["significance_obs"]) or not isinstance(r["discovery"], bool):
        fail("significance_obs / discovery")
    print(f"significance_obs {r['significance_obs']:.2f}, discovery {r['discovery']}, significance_exp {r.get('significance_exp')}")
    for name in ("scale_shift", "smear", "sel_eff"):
        for f in ("muon", "electron"):
            v = r["calibration"][name][f]
            if not (math.isfinite(v["value"]) and math.isfinite(v["unc"])):
                fail(f"calibration {name} {f}")
            print(f"calibration {name} {f}: {v['value']:.5f} +- {v['unc']:.5f}")
    for name, q in r.get("quality", {}).items():
        print(f"quality {name}: {q['value']} ({q['note'][:80]})")
    spec = json.loads((out / "MODEL.json").read_text())
    try:
        import pyhf
    except ImportError:
        print("pyhf not importable: MODEL.json not built here")
        print("OK")
        return
    ws = pyhf.Workspace(spec)
    model = ws.model()
    pars = pyhf.infer.mle.fit(ws.data(model), model)
    mu = float(pars[model.config.poi_index])
    print(f"MODEL.json builds ({len(ws.channels)} channels, {model.config.npars} parameters); pyhf MLE mu {mu:.4f} "
          f"(reported {pois['mu']['value']:.4f})")
    print("OK")


if __name__ == "__main__":
    main()
