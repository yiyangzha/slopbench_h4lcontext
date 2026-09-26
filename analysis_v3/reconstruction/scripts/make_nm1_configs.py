"""The N-1 variants of the final selection configuration (eval_selfreport N-1 expected significances).

    pixi run py -- analysis_v3/reconstruction/scripts/make_nm1_configs.py

Each variant removes one requirement of analysis_v3/reconstruction/config/selection_ul16_v2.json and
keeps the others.  A requirement is relaxed to the floor of the event records (h4l_reco,
reco_ul16_v2.json) where the records carry one, so "without" means "down to the record floor":
    iso     FSR-subtracted isolation < 0.35 -> no isolation requirement;
    sip     SIP < 4 -> SIP < 8 (record floor);
    z1      40 < m_Z1 < 120 GeV -> 30 < m_Z1 < 130 GeV (record floor);
    z2      m_Z2 > 10 GeV (signal region, also in the smart cut) and m_Z2 < 120 GeV -> 3 < m_Z2 < 130 GeV
            (record floor; the Z -> 4l region follows, since it must contain the signal region);
    pt      lepton pT > 5 (muons) / 7 (electrons) GeV -> 3 / 5 GeV (record floor);
    leadpt  leading / subleading lepton pT > 20 / 10 GeV -> 15 / 7 GeV (record floor);
    osmass  every opposite-sign pair m_ll > 4 GeV -> 2 GeV (record floor);
    muid    the muon identification (AN tight: PF, or tracker high-pT) -> none (leptons.muon.require_id);
    elid    the electron identification (mvaFall17V2noIso WPL) -> none (leptons.electron.id "none");
the records keep every loose muon and electron, the identification acting only in h4l_select (the
options exist from the h4l_select of the N-1 selections on).  The trigger, the vertex, |eta|, dxy and
dz, the lepton dR (their record floors are the cuts) and m4l > 70 GeV (no effect inside 118-130 GeV)
have no variant.  The z2 variant's effective floor is 4 GeV (the kept m_ll(OS) > 4 GeV), the Z upper
edges do not act inside the window.  Writes the missing
analysis_v3/reconstruction/config/nm1/selection_ul16_v2_nm1_<cut>.json; an existing file must equal
what the script would write (never overwritten).
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
BASE = REPO / "analysis_v3/reconstruction/config/selection_ul16_v2.json"
RECO = REPO / "analysis_v3/reconstruction/config/reco_ul16_v2.json"
OUT = REPO / "analysis_v3/reconstruction/config/nm1"


def main() -> int:
    base = json.loads(BASE.read_text(encoding="utf-8"))
    reco = json.loads(RECO.read_text(encoding="utf-8"))
    floor_l, floor_c = reco["leptons"], reco["candidates"]
    variants = {}

    def variant(name: str, note: str):
        c = copy.deepcopy(base)
        c["version"] = f"{base['version']}_nm1_{name}"
        c["description"] = f"N-1 variant of {base['version']} (make_nm1_configs.py): {note}"
        variants[name] = c
        return c

    variant("iso", "no isolation requirement")["leptons"]["max_iso_fsr"] = 1.0e9
    variant("sip", f"SIP < {floor_c['max_sip']} (event-record floor) instead of {base['leptons']['max_sip']}")["leptons"]["max_sip"] = \
        floor_c["max_sip"]
    variant("z1", f"{floor_c['z1_mass'][0]} < m_Z1 < {floor_c['z1_mass'][1]} GeV (event-record floor)")["candidate"]["z1_mass"] = \
        list(floor_c["z1_mass"])
    c = variant("z2", f"{floor_c['z2_mass'][0]} < m_Z2 < {floor_c['z2_mass'][1]} GeV in every region (event-record floor)")
    c["candidate"]["z2_high"] = floor_c["z2_mass"][1]
    for region in c["candidate"]["regions"]:
        region["z2_low"] = floor_c["z2_mass"][0]
    c = variant("pt", f"lepton pT > {floor_l['muon_floor_pt']} (muons) / {floor_l['electron_floor_pt']} (electrons) GeV "
                      "(event-record floor)")
    c["leptons"]["muon"]["pt"] = floor_l["muon_floor_pt"]
    c["leptons"]["electron"]["pt"] = floor_l["electron_floor_pt"]
    c = variant("leadpt", f"leading / subleading pT > {floor_c['leading_pt']} / {floor_c['subleading_pt']} GeV (event-record floor)")
    c["candidate"]["leading_pt"] = floor_c["leading_pt"]
    c["candidate"]["subleading_pt"] = floor_c["subleading_pt"]
    variant("osmass", f"m_ll(OS) > {floor_c['min_os_mass']} GeV (event-record floor)")["candidate"]["min_os_mass"] = floor_c["min_os_mass"]
    variant("muid", "no muon identification requirement (loose muons of the records)")["leptons"]["muon"]["require_id"] = False
    variant("elid", "no electron identification requirement (loose electrons of the records)")["leptons"]["electron"]["id"] = "none"
    OUT.mkdir(exist_ok=True)
    for name, config in variants.items():
        path = OUT / f"selection_ul16_v2_nm1_{name}.json"
        text = json.dumps(config, indent=1) + "\n"
        if path.exists():
            if path.read_text(encoding="utf-8") != text:
                raise SystemExit(f"{path} exists with another content")
            print(f"[nm1] {path} (exists, identical)")
            continue
        path.write_text(text, encoding="utf-8")
        print(f"[nm1] {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
