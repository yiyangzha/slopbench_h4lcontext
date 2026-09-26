// tnp_extract: the compact tag-and-probe pair ntuple.
//
// For the stage-2 skim outputs of one task (all of one sample) it reads the
// Pairs tree (every opposite-sign same-flavour pair above the skim floors)
// and keeps the pairs with a raw mass in the extraction window and at least
// one tag candidate leg: the event fired a tag path of the pair's flavour,
// the leg is matched to an HLT object carrying a tag filter bit with a
// trigger-object pT above trig_pt_min, and its raw pT exceeds pt_min.  No
// identification, isolation or SIP requirement is applied, so every probe
// denominator of the tag-and-probe stays unbiased.  The final tag and probe
// selection (on the calibrated momenta) is applied later.
//
// The FSR-subtracted isolation of a leg (AN-16-442 section 4) is
//   iso_fsr = iso_chg + max(0, (iso_all - iso_chg) - sum pT_gamma / pT),
// summed over the FSR photons of both legs (muons only) with
// 0.01 < dR(gamma, leg) < 0.3.
//
// Output trees: TnPPairs (one row per kept pair; the bits of tag_candidates
// mark the legs that are tag candidates, 1 = leg 1, 2 = leg 2), Inputs (one
// row per input skim).  Every leg keeps the skim fields the tag-and-probe
// needs, including the electron cut-based ID and MVA inputs and the muon ID
// inputs for an identification re-tuned in stage 4b.  Any read error and any disagreement with the planned entry counts
// aborts the task.
//
// Usage: tnp_extract --task TASK.json --out-json OUT.json --out-root OUT.root

#include "h4l/hash.h"
#include "h4l/io.h"
#include "h4l/kinematics.h"
#include "h4l/root_io.h"

#include <TError.h>
#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <memory>
#include <string>

namespace {
using h4l::json;

struct InputError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

bool g_read_error = false;
std::string g_read_message;
void record_errors(int level, Bool_t abort, const char* location, const char* message) {
  if (level >= kError) {
    g_read_error = true;
    g_read_message = std::string(location ? location : "?") + ": " + (message ? message : "");
  }
  DefaultErrorHandler(level, abort, location, message);
}
void check_read(const std::string& where) {
  if (g_read_error) throw InputError("ROOT error while reading " + where + ": " + g_read_message);
}

std::string hex64(std::uint64_t value) {
  char buffer[17];
  std::snprintf(buffer, sizeof(buffer), "%016llx", static_cast<unsigned long long>(value));
  return buffer;
}

// The float fields of a leg that are copied as they are (the skim names).
// Electron ID inputs (-9 for muons) and muon ID inputs (-9 for electrons) are
// kept for an identification trained or re-tuned in stage 4b.
const char* const kFloatFields[] = {
    "pt", "eta", "eta_sc", "phi", "pt_err", "energy_err", "iso03", "iso03_chg", "iso04", "mini_iso", "sip", "ip3d", "dxy", "dz",
    "dxy_err", "dz_err", "mva_noiso", "mva_iso", "r9", "sieie", "hoe", "e_inv_minus_p_inv", "sc_et_over_pt", "dr03_ecal",
    "dr03_hcal", "dr03_tk", "e_corr", "tunep_rel_pt", "segment_comp", "tk_rel_iso", "mva_lowpt", "soft_mva", "gen_pt",
    "fsr_pt", "fsr_eta", "fsr_phi", "trig_pt"};
constexpr int kFloats = sizeof(kFloatFields) / sizeof(kFloatFields[0]);
// Positions of the fields used here (checked against kFloatFields at start).
enum Field { kPt = 0, kEta = 1, kPhi = 3, kIso03 = 6, kIso03Chg = 7, kFsrPt = 33, kFsrEta = 34, kFsrPhi = 35, kTrigPt = 36 };
void check_fields() {
  const std::pair<Field, const char*> expected[] = {{kPt, "pt"}, {kEta, "eta"}, {kPhi, "phi"}, {kIso03, "iso03"},
                                                    {kIso03Chg, "iso03_chg"}, {kFsrPt, "fsr_pt"}, {kFsrEta, "fsr_eta"},
                                                    {kFsrPhi, "fsr_phi"}, {kTrigPt, "trig_pt"}};
  for (const auto& [index, name] : expected)
    if (std::string(kFloatFields[index]) != name) throw std::runtime_error(std::string("field table out of order at ") + name);
}

struct Leg {
  Short_t index = -1;
  Float_t f[kFloats] = {};
  UShort_t flags = 0;
  Int_t charge = 0, trig_bits = 0, vid_bitmap = 0, tight_charge = -1;
  UChar_t gen_flav = 0, lost_hits = 0, cut_based = 0, n_stations = 0, n_tracker_layers = 0;
  Float_t iso_fsr = 0;
  Char_t charge_out = 0;

  void address(TTree& tree, const std::string& prefix) {
    auto set = [&](const std::string& field, void* where) {
      const std::string name = prefix + field;
      tree.SetBranchStatus(name.c_str(), true);
      if (tree.SetBranchAddress(name.c_str(), where) < 0) throw InputError("cannot address " + name);
    };
    set("index", &index);
    for (int k = 0; k < kFloats; ++k) set(kFloatFields[k], &f[k]);
    set("flags", &flags);
    set("charge", &charge);
    set("trig_bits", &trig_bits);
    set("vid_bitmap", &vid_bitmap);
    set("tight_charge", &tight_charge);
    set("gen_flav", &gen_flav);
    set("lost_hits", &lost_hits);
    set("cut_based", &cut_based);
    set("n_stations", &n_stations);
    set("n_tracker_layers", &n_tracker_layers);
  }
  void book(TTree& tree, const std::string& prefix) {
    auto branch = [&](const std::string& field, void* where, const char* type) {
      tree.Branch((prefix + field).c_str(), where, (prefix + field + "/" + type).c_str());
    };
    branch("index", &index, "S");
    for (int k = 0; k < kFloats; ++k) branch(kFloatFields[k], &f[k], "F");
    branch("iso_fsr", &iso_fsr, "F");
    branch("flags", &flags, "s");
    branch("charge", &charge_out, "B");
    branch("trig_bits", &trig_bits, "I");
    branch("vid_bitmap", &vid_bitmap, "I");
    branch("tight_charge", &tight_charge, "I");
    branch("gen_flav", &gen_flav, "b");
    branch("lost_hits", &lost_hits, "b");
    branch("cut_based", &cut_based, "b");
    branch("n_stations", &n_stations, "b");
    branch("n_tracker_layers", &n_tracker_layers, "b");
  }
};

float fsr_isolation(const Leg& leg, const Leg (&legs)[2]) {
  double photons = 0;
  for (const Leg& owner : legs) {
    if (!(owner.f[kFsrPt] > 0)) continue;
    const double dr = h4l::delta_r(leg.f[kEta], leg.f[kPhi], owner.f[kFsrEta], owner.f[kFsrPhi]);
    if (dr > 0.01 && dr < 0.3) photons += owner.f[kFsrPt];
  }
  const double chg = leg.f[kIso03Chg];
  const double neutral = std::max(0.0, static_cast<double>(leg.f[kIso03]) - chg - photons / leg.f[kPt]);
  return static_cast<float>(chg + neutral);
}

struct TagRule {
  UInt_t path_mask = 0;
  int filter_bits = 0;
  double trig_pt_min = 0, pt_min = 0;
  bool candidate(const Leg& leg, UInt_t trig) const {
    return (trig & path_mask) != 0 && (leg.trig_bits & filter_bits) != 0 && leg.f[kTrigPt] >= trig_pt_min &&
           leg.f[kPt] >= pt_min;
  }
};
}  // namespace

int main(int argc, char** argv) {
  std::string task_path, out_json, out_root;
  for (int index = 1; index < argc; index += 2) {
    const std::string key = argv[index];
    if (index + 1 >= argc) {
      std::cerr << "missing value for " << key << "\n";
      return 64;
    }
    if (key == "--task") task_path = argv[index + 1];
    else if (key == "--out-json") out_json = argv[index + 1];
    else if (key == "--out-root") out_root = argv[index + 1];
    else {
      std::cerr << "unknown argument " << key << "\n";
      return 64;
    }
  }
  if (task_path.empty() || out_json.empty() || out_root.empty()) {
    std::cerr << "usage: tnp_extract --task TASK.json --out-json OUT.json --out-root OUT.root\n";
    return 64;
  }
  gROOT->SetBatch(true);
  SetErrorHandler(record_errors);
  try {
    check_fields();
    const json task = h4l::read_json(task_path);
    const json& extract = task.at("extract_config");
    const bool is_mc = task.at("config").at("kind").get<std::string>() == "mc";
    const double mass_low = extract.at("mass_window").at(0).get<double>();
    const double mass_high = extract.at("mass_window").at(1).get<double>();
    std::map<int, TagRule> rules;
    for (const auto& [name, pdg] : {std::pair{"muon", 13}, std::pair{"electron", 11}}) {
      TagRule rule;
      rule.path_mask = task.at("tag_path_masks").at(name).get<UInt_t>();
      if (rule.path_mask == 0) throw std::runtime_error(std::string("empty tag path mask for ") + name);
      rule.filter_bits = extract.at(name).at("filter_bits").get<int>();
      rule.trig_pt_min = extract.at(name).at("trig_pt_min").get<double>();
      rule.pt_min = extract.at(name).at("pt_min").get<double>();
      rules[pdg] = rule;
    }

    std::unique_ptr<TFile> output(TFile::Open(out_root.c_str(), "RECREATE"));
    if (!output || output->IsZombie()) throw std::runtime_error("cannot create " + out_root);
    output->cd();
    TTree& pairs_out = *new TTree("TnPPairs", "tag-and-probe pairs");
    TTree& inputs_tree = *new TTree("Inputs", "per-input counters");
    ULong64_t file_key = 0;
    Long64_t entry = 0;
    Float_t weight = 1, mass = 0, pair_pt = 0;
    Int_t flavour = 0, npv = 0;
    Char_t flavour_out = 0;
    UInt_t trig = 0;
    UChar_t tag_mask = 0;
    Leg leg[2];
    pairs_out.Branch("file_key", &file_key, "file_key/l");
    pairs_out.Branch("entry", &entry, "entry/L");
    pairs_out.Branch("weight", &weight, "weight/F");
    pairs_out.Branch("flavour", &flavour_out, "flavour/B");
    pairs_out.Branch("mass", &mass, "mass/F");
    pairs_out.Branch("pair_pt", &pair_pt, "pair_pt/F");
    pairs_out.Branch("trig", &trig, "trig/i");
    pairs_out.Branch("npv", &npv, "npv/I");
    pairs_out.Branch("tag_candidates", &tag_mask, "tag_candidates/b");
    leg[0].book(pairs_out, "l1_");
    leg[1].book(pairs_out, "l2_");
    std::string input_path;
    Long64_t input_pairs = 0, input_kept = 0;
    inputs_tree.Branch("path", &input_path);
    inputs_tree.Branch("pairs", &input_pairs, "pairs/L");
    inputs_tree.Branch("kept", &input_kept, "kept/L");

    json inputs_report = json::array();
    long long total_read = 0, total_kept = 0, kept_mm = 0, kept_ee = 0;
    for (const auto& input : task.at("inputs")) {
      input_path = input.at("root").get<std::string>();
      const long long expected = input.at("pairs").get<long long>();
      std::string open_error;
      auto file = h4l::open_input(input_path, 5, open_error);
      if (!file) throw InputError("cannot open " + input_path + ": " + open_error);
      g_read_error = false;
      auto* tree = dynamic_cast<TTree*>(file->Get("Pairs"));
      if (!tree) throw InputError("no Pairs in " + input_path);
      if (tree->GetEntries() != expected)
        throw InputError(input_path + " has " + std::to_string(tree->GetEntries()) + " pairs, the skim report " +
                         std::to_string(expected));
      tree->SetBranchStatus("*", false);
      auto set = [&](const char* name, void* where) {
        tree->SetBranchStatus(name, true);
        if (tree->SetBranchAddress(name, where) < 0) throw InputError(std::string("cannot address ") + name);
      };
      set("file_key", &file_key);
      set("entry", &entry);
      set("weight", &weight);
      set("flavour", &flavour);
      set("mass", &mass);
      set("pair_pt", &pair_pt);
      set("trig", &trig);
      set("npv", &npv);
      leg[0].address(*tree, "l1_");
      leg[1].address(*tree, "l2_");
      input_pairs = tree->GetEntries();
      input_kept = 0;
      for (Long64_t row = 0; row < input_pairs; ++row) {
        if (tree->GetEntry(row) <= 0) throw InputError("cannot read Pairs row " + std::to_string(row) + " of " + input_path);
        check_read(input_path + " row " + std::to_string(row));
        const auto rule = rules.find(flavour);
        if (rule == rules.end()) throw InputError("unknown pair flavour in " + input_path);
        if (mass <= mass_low || mass >= mass_high) continue;
        tag_mask = 0;
        for (int l = 0; l < 2; ++l)
          if (rule->second.candidate(leg[l], trig)) tag_mask |= static_cast<UChar_t>(1 << l);
        if (tag_mask == 0) continue;
        for (auto& l : leg) {
          l.iso_fsr = fsr_isolation(l, leg);
          l.charge_out = static_cast<Char_t>(l.charge);
        }
        if (!is_mc) weight = 1;
        flavour_out = static_cast<Char_t>(flavour);
        if (pairs_out.Fill() <= 0) throw std::runtime_error("TnPPairs Fill failed");
        ++input_kept;
        (flavour == 13 ? kept_mm : kept_ee) += 1;
      }
      check_read(input_path);
      total_read += input_pairs;
      total_kept += input_kept;
      if (inputs_tree.Fill() <= 0) throw std::runtime_error("Inputs Fill failed");
      inputs_report.push_back({{"root", input_path}, {"skim_task_id", input.at("skim_task_id")},
                               {"root_sha256", input.at("root_sha256")}, {"pairs", input_pairs}, {"kept", input_kept}});
      std::cout << "[extract] " << input_path << ": " << input_pairs << " pairs, " << input_kept << " kept\n";
      file.reset();
    }

    output->cd();
    for (TTree* tree : {&pairs_out, &inputs_tree})
      if (tree->Write() <= 0) throw std::runtime_error(std::string("cannot write ") + tree->GetName());
    output->Close();
    output.reset();
    if (g_read_error) throw std::runtime_error("ROOT error while writing the output: " + g_read_message);
    const std::map<std::string, long long> trees = {{"TnPPairs", total_kept},
                                                    {"Inputs", static_cast<long long>(task.at("inputs").size())}};
    h4l::validate_root_output_full(out_root, trees);
    json report = {{"schema", "h4l_v3_tnp_extract_report/1"},
                   {"task_id", task.at("task_id")},
                   {"sample", task.at("config").at("sample")},
                   {"kind", task.at("config").at("kind")},
                   {"extract_config_fnv1a64", hex64(h4l::fnv1a64(extract.dump()))},
                   {"tag_path_masks", task.at("tag_path_masks")},
                   {"frozen_program", task.value("frozen_program", json())},
                   {"inputs", inputs_report},
                   {"inputs_root_sha256_origin", "planned (copied from the skim publication records)"},
                   {"tag_candidates_bits", "1 = leg 1 (l1_), 2 = leg 2 (l2_)"},
                   {"source_file_keys", task.at("source_file_keys")},
                   {"totals", {{"pairs_read", total_read}, {"kept", total_kept}, {"kept_mm", kept_mm}, {"kept_ee", kept_ee}}},
                   {"root_output", {{"trees", trees}, {"objects", json::array()}}},
                   {"finished_utc", h4l::utc_now()}};
    h4l::write_json_file(out_json, report);
  } catch (const InputError& error) {
    std::cerr << "ERROR (input): " << error.what() << "\n";
    return 3;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
