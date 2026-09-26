// calib_extract: the compact Z -> ll pair ntuple of the lepton calibration.
//
// For the stage-2 skim outputs of one task (all of one sample) it reads the
// Pairs tree and keeps the pairs whose event fired the analysis trigger OR,
// whose raw mass lies in the extraction window and whose two legs pass the
// loose extraction preselection (FSR-subtracted isolation, SIP).  The kept
// pairs carry everything the calibration needs to apply its own tight
// selection, to transform the leg momenta, to seed the per-lepton smearing
// (original file key, original entry, original lepton index), to split the
// per-event mass-uncertainty categories and to draw the validation plots.
//
// The FSR-subtracted isolation of a leg (AN-16-442 section 4) is
//   iso_fsr = iso_chg + max(0, (iso_all - iso_chg) - sum pT_gamma / pT),
// summed over the FSR photons of both legs (muons only) with
// 0.01 < dR(gamma, leg) < 0.3.
//
// Output trees:
//   CalibPairs  one row per kept pair;
//   Inputs      one row per input skim (path, pairs read, pairs kept).
// Any read error (GetEntry <= 0 or a ROOT error message while reading) and
// any disagreement with the planned entry counts aborts the task.
//
// Usage: calib_extract --task TASK.json --out-json OUT.json --out-root OUT.root

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

// One leg of the skim Pairs tree (the fields read here) and of CalibPairs.
struct Leg {
  Short_t index = -1;
  Float_t pt = 0, eta = 0, eta_sc = 0, phi = 0, pt_err = 0, iso03 = 0, iso03_chg = 0, sip = 0, dxy = 0, dz = 0,
          mva_noiso = -9, mva_iso = -9, r9 = -9, gen_pt = -1, fsr_pt = -1, fsr_eta = 0, fsr_phi = 0;
  UShort_t flags = 0;
  Int_t charge = 0, trig_bits = 0;
  UChar_t gen_flav = 0;
  // Written only.
  Float_t iso_fsr = 0;
  Char_t charge_out = 0;

  void address(TTree& tree, const std::string& prefix) {
    auto set = [&](const std::string& field, void* where) {
      const std::string name = prefix + field;
      tree.SetBranchStatus(name.c_str(), true);
      if (tree.SetBranchAddress(name.c_str(), where) < 0) throw InputError("cannot address " + name);
    };
    set("index", &index);
    for (auto [field, where] : std::initializer_list<std::pair<const char*, Float_t*>>{
             {"pt", &pt}, {"eta", &eta}, {"eta_sc", &eta_sc}, {"phi", &phi}, {"pt_err", &pt_err}, {"iso03", &iso03},
             {"iso03_chg", &iso03_chg}, {"sip", &sip}, {"dxy", &dxy}, {"dz", &dz}, {"mva_noiso", &mva_noiso},
             {"mva_iso", &mva_iso}, {"r9", &r9}, {"gen_pt", &gen_pt}, {"fsr_pt", &fsr_pt}, {"fsr_eta", &fsr_eta},
             {"fsr_phi", &fsr_phi}})
      set(field, where);
    set("flags", &flags);
    set("charge", &charge);
    set("trig_bits", &trig_bits);
    set("gen_flav", &gen_flav);
  }
  void book(TTree& tree, const std::string& prefix) {
    auto branch = [&](const std::string& field, void* where, const char* type) {
      tree.Branch((prefix + field).c_str(), where, (prefix + field + "/" + type).c_str());
    };
    branch("index", &index, "S");
    for (auto [field, where] : std::initializer_list<std::pair<const char*, Float_t*>>{
             {"pt", &pt}, {"eta", &eta}, {"eta_sc", &eta_sc}, {"phi", &phi}, {"pt_err", &pt_err}, {"iso03", &iso03},
             {"iso03_chg", &iso03_chg}, {"iso_fsr", &iso_fsr}, {"sip", &sip}, {"dxy", &dxy}, {"dz", &dz},
             {"mva_noiso", &mva_noiso}, {"mva_iso", &mva_iso}, {"r9", &r9}, {"gen_pt", &gen_pt}, {"fsr_pt", &fsr_pt},
             {"fsr_eta", &fsr_eta}, {"fsr_phi", &fsr_phi}})
      branch(field, where, "F");
    branch("flags", &flags, "s");
    branch("charge", &charge_out, "B");
    branch("trig_bits", &trig_bits, "I");
    branch("gen_flav", &gen_flav, "b");
  }
};

// FSR-subtracted relative isolation of `leg` with the FSR photons of both legs.
float fsr_isolation(const Leg& leg, const Leg (&legs)[2]) {
  double photons = 0;
  for (const Leg& owner : legs) {
    if (!(owner.fsr_pt > 0)) continue;
    const double dr = h4l::delta_r(leg.eta, leg.phi, owner.fsr_eta, owner.fsr_phi);
    if (dr > 0.01 && dr < 0.3) photons += owner.fsr_pt;
  }
  const double neutral = std::max(0.0, static_cast<double>(leg.iso03) - leg.iso03_chg - photons / leg.pt);
  return static_cast<float>(leg.iso03_chg + neutral);
}
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
    std::cerr << "usage: calib_extract --task TASK.json --out-json OUT.json --out-root OUT.root\n";
    return 64;
  }
  gROOT->SetBatch(true);
  SetErrorHandler(record_errors);
  try {
    const json task = h4l::read_json(task_path);
    const json& extract = task.at("extract_config");
    const bool is_mc = task.at("config").at("kind").get<std::string>() == "mc";
    const UInt_t trigger_mask = task.at("trigger_or_mask").get<UInt_t>();
    if (trigger_mask == 0) throw std::runtime_error("empty trigger OR mask");
    const double mass_low = extract.at("mass_window").at(0).get<double>();
    const double mass_high = extract.at("mass_window").at(1).get<double>();
    const double max_iso = extract.at("max_iso_fsr").get<double>();
    const double max_sip = extract.at("max_sip").get<double>();

    std::unique_ptr<TFile> output(TFile::Open(out_root.c_str(), "RECREATE"));
    if (!output || output->IsZombie()) throw std::runtime_error("cannot create " + out_root);
    output->cd();
    TTree& pairs_out = *new TTree("CalibPairs", "Z->ll calibration pairs");
    TTree& inputs_tree = *new TTree("Inputs", "per-input counters");
    ULong64_t file_key = 0;
    Long64_t entry = 0;
    Float_t weight = 1, mass = 0, mass_fsr = 0, pair_pt = 0, pair_y = 0;
    Int_t flavour = 0, npv = 0;
    Char_t flavour_out = 0;
    UInt_t trig = 0;
    Leg leg[2];
    pairs_out.Branch("file_key", &file_key, "file_key/l");
    pairs_out.Branch("entry", &entry, "entry/L");
    pairs_out.Branch("weight", &weight, "weight/F");
    pairs_out.Branch("flavour", &flavour_out, "flavour/B");
    pairs_out.Branch("mass", &mass, "mass/F");
    pairs_out.Branch("mass_fsr", &mass_fsr, "mass_fsr/F");
    pairs_out.Branch("pair_pt", &pair_pt, "pair_pt/F");
    pairs_out.Branch("pair_y", &pair_y, "pair_y/F");
    pairs_out.Branch("trig", &trig, "trig/i");
    pairs_out.Branch("npv", &npv, "npv/I");
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
      set("mass_fsr", &mass_fsr);
      set("pair_pt", &pair_pt);
      set("pair_y", &pair_y);
      set("trig", &trig);
      set("npv", &npv);
      leg[0].address(*tree, "l1_");
      leg[1].address(*tree, "l2_");
      input_pairs = tree->GetEntries();
      input_kept = 0;
      for (Long64_t row = 0; row < input_pairs; ++row) {
        if (tree->GetEntry(row) <= 0) throw InputError("cannot read Pairs row " + std::to_string(row) + " of " + input_path);
        check_read(input_path + " row " + std::to_string(row));
        if (flavour != 11 && flavour != 13) throw InputError("unknown pair flavour in " + input_path);
        if ((trig & trigger_mask) == 0) continue;
        if (mass <= mass_low || mass >= mass_high) continue;
        bool pass = true;
        for (auto& l : leg) {
          l.iso_fsr = fsr_isolation(l, leg);
          pass = pass && l.iso_fsr < max_iso && l.sip < max_sip;
        }
        if (!pass) continue;
        if (!is_mc) weight = 1;
        flavour_out = static_cast<Char_t>(flavour);
        for (auto& l : leg) l.charge_out = static_cast<Char_t>(l.charge);
        if (pairs_out.Fill() <= 0) throw std::runtime_error("CalibPairs Fill failed");
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
    const std::map<std::string, long long> trees = {{"CalibPairs", total_kept},
                                                    {"Inputs", static_cast<long long>(task.at("inputs").size())}};
    h4l::validate_root_output_full(out_root, trees);

    json report = {{"schema", "h4l_v3_calib_extract_report/2"},
                   {"task_id", task.at("task_id")},
                   {"sample", task.at("config").at("sample")},
                   {"kind", task.at("config").at("kind")},
                   {"extract_config_fnv1a64", hex64(h4l::fnv1a64(extract.dump()))},
                   {"trigger_or_mask", trigger_mask},
                   {"frozen_program", task.value("frozen_program", json())},
                   {"inputs", inputs_report},
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
