"""The LHCHXSWG YR4 cross sections and H -> ZZ -> 4l branching ratio against m_H (stage 6).

    pixi run py -- analysis_v3/signal_model/scripts/yr4_xsbr.py --version v1

Reads the local copy of the public YR4 spreadsheet
(slopbench_code_fork/results/h4l_seeds/v3/inputs/Higgs_XSBR_YR4_update.xlsx, the one
file under results/ that may be read), verifies its sha256 against the cell map of
benchmark/configs/h4l_pfnano_seed_plan_v3.json ("signal_cross_sections"), and
extracts with the Python standard library (an .xlsx file is a zip archive of XML) the
columns of that map (1-based: A = 1): the 13 TeV ggF, VBF, WH and ZH cross sections
and BR(H -> ZZ -> 4l, l = e, mu).  The effective cross section of a production mode
at m_H is (the documented formula of the seed plan)
    sigma_eff(mode, m_H) = sigma_eff(mode, 125) sigma_YR4(mode, m_H) / sigma_YR4(mode, 125)
                           BR4l(m_H) / BR4l(125),
with VH = WH + ZH.  Writes production_v3/signal_model/yr4/<version>/yr4.json (the
tables and sigma_eff(mode, m_H) on the spreadsheet mass grid).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
SPREADSHEET = Path("/eos/user/y/yiyangz/codex/slopbench_code_fork/results/h4l_seeds/v3/inputs/Higgs_XSBR_YR4_update.xlsx")
SEED_PLAN = Path("/eos/user/y/yiyangz/codex/slopbench_code_fork/benchmark/configs/h4l_pfnano_seed_plan_v3.json")
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
REL_NS = "{http://schemas.openxmlformats.org/package/2006/relationships}"


def column_index(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference).group(0)
    index = 0
    for letter in letters:
        index = index * 26 + (ord(letter) - ord("A") + 1)
    return index


def read_sheets(path: Path, wanted: set[str]) -> dict:
    """{sheet name: {row number: {column (1-based): value}}} for the wanted sheets."""
    with zipfile.ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall("m:si", NS):
                shared.append("".join(t.text or "" for t in item.iter(f"{{{NS['m']}}}t")))
        workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
        rels = ElementTree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {rel.get("Id"): rel.get("Target") for rel in rels.findall(f"{REL_NS}Relationship")}
        out = {}
        for sheet in workbook.find("m:sheets", NS).findall("m:sheet", NS):
            name = sheet.get("name")
            if name not in wanted:
                continue
            target = targets[sheet.get(f"{{{NS['r']}}}id")]
            target = target.lstrip("/")
            target = target if target.startswith("xl/") else "xl/" + target
            root = ElementTree.fromstring(archive.read(target))
            rows = {}
            for row in root.iter(f"{{{NS['m']}}}row"):
                cells = {}
                for cell in row.findall("m:c", NS):
                    value = cell.find("m:v", NS)
                    if value is None or value.text is None:
                        continue
                    text = value.text
                    if cell.get("t") == "s":
                        text = shared[int(text)]
                    cells[column_index(cell.get("r"))] = text
                rows[int(row.get("r"))] = cells
            out[name] = rows
    return out


def column_pairs(rows: dict, mass_column: int, value_column: int) -> list[tuple[float, float]]:
    pairs = []
    for number in sorted(rows):
        cells = rows[number]
        try:
            mass, value = float(cells[mass_column]), float(cells[value_column])
        except (KeyError, ValueError):
            continue
        pairs.append((mass, value))
    return pairs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True)
    args = parser.parse_args()
    out_dir = PRODUCTION / "signal_model" / "yr4" / args.version
    if out_dir.exists():
        raise SystemExit(f"{out_dir} exists")
    plan = json.loads(SEED_PLAN.read_text(encoding="utf-8"))["signal_cross_sections"]
    digest = hashlib.sha256(SPREADSHEET.read_bytes()).hexdigest()
    if digest != plan["source"]["sha256"]:
        raise SystemExit(f"spreadsheet sha256 {digest} differs from the seed plan {plan['source']['sha256']}")
    cells = plan["source"]["cells"]
    sheets = read_sheets(SPREADSHEET, {c["sheet"] for c in cells.values()})
    tables = {}
    for key, cell in cells.items():
        pairs = column_pairs(sheets[cell["sheet"]], cell["mass_column"], cell["value_column"])
        # Keep the monotonic mass grid of the Higgs-mass table (the columns may hold other tables further down).
        grid = []
        for mass, value in pairs:
            if 100.0 <= mass <= 200.0 and (not grid or mass > grid[-1][0]):
                grid.append((mass, value))
        tables[key] = grid
    anchor = plan["anchor_effective_xsec_pb"]
    common = sorted(set(m for m, _ in tables["ggF"]) & set(m for m, _ in tables["VBF"]) & set(m for m, _ in tables["WH"])
                    & set(m for m, _ in tables["ZH"]) & set(m for m, _ in tables["BR4l"]))
    lookup = {key: dict(grid) for key, grid in tables.items()}
    if 125.0 not in common:
        raise SystemExit("125 GeV is not on the common mass grid")
    effective = {}
    for mode, parts in (("GluGluToHToZZ", ["ggF"]), ("VBF_HToZZ", ["VBF"]), ("VHToZZ", ["WH", "ZH"])):
        ref = sum(lookup[p][125.0] for p in parts)
        effective[mode] = [{"mH": m, "sigma_eff_pb": anchor[mode] * sum(lookup[p][m] for p in parts) / ref
                            * lookup["BR4l"][m] / lookup["BR4l"][125.0]} for m in common]
    out_dir.mkdir(parents=True)
    report = {"schema": "h4l_v3_yr4/1", "spreadsheet": str(SPREADSHEET), "sha256": digest, "cell_map": cells,
              "anchor_effective_xsec_pb": anchor, "formula": plan["_description"], "tables": tables, "mass_grid": common,
              "sigma_eff": effective}
    (out_dir / "yr4.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    for key, grid in tables.items():
        at125 = dict(grid).get(125.0)
        print(f"[yr4] {key}: {len(grid)} points {grid[0][0]}-{grid[-1][0]} GeV, value at 125 GeV {at125}")
    for mode, rows in effective.items():
        values = {r["mH"]: r["sigma_eff_pb"] for r in rows}
        print(f"[yr4] sigma_eff {mode}: 120 {values.get(120.0):.6g}, 125 {values.get(125.0):.6g}, 130 {values.get(130.0):.6g} pb")
    print(f"[yr4] {out_dir / 'yr4.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
