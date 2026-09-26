// skim_v3: the one-pass stage-2 skim of the UL16 PFNano v3 analysis.
//
// For the input files of one task (all of one sample) it writes one ROOT
// file with
//   Pairs         every opposite-sign same-flavour pair of leptons above the
//                 skim floors with a mass in the configured window, with both
//                 legs' kinematics, momentum errors, identification flags,
//                 scores and ID inputs, isolation, impact parameters, trigger
//                 match and (MC) generator match: the input of the lepton
//                 calibration and of tag-and-probe;
//   PhotonPairs   (tag electron, Photon) pairs: electron reconstruction probes;
//   TrigObjPairs  (tag muon, merged HLT muon object) pairs: muon
//                 reconstruction probes;
//   Events        a slim NanoAOD copy (configured branches) of every event
//                 with at least N leptons above the floors, with its original
//                 file key and entry in h4l_file_key and h4l_entry;
//   GenTable      (signal MC only) one row per event: the Higgs, its decay
//                 leptons, the dressed prompt leptons, the fiducial selection
//                 of JHEP 11 (2017) 047 Table 4, generator jets, the STXS
//                 stage-0 bin and the VH decay class;
//   Files         per-input counters.
// The floors are below every threshold the analysis may use, so later
// stages can move thresholds without re-skimming.  Every read error (a bad
// TTreeReader status, GetEntry <= 0, or any ROOT error message while events
// are read) aborts the task without output; the task is then retried.
//
// Usage: skim_v3 --task TASK.json --out-json OUT.json --out-root OUT.root

#include "h4l/hash.h"
#include "h4l/io.h"
#include "h4l/kinematics.h"
#include "h4l/root_io.h"
#include "h4l/trigger.h"

#include <TChain.h>
#include <TError.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <fnmatch.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace {
using h4l::json;
using h4l::P4;

struct InputError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// Any ROOT error while events are read (for example a basket that fails to
// decompress) marks the read as failed.
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

template <typename T>
using Array = TTreeReaderArray<T>;
template <typename T>
using Value = TTreeReaderValue<T>;

// ----------------------------------------------------------------------------
// Configuration
struct Floor {
  double pt, abs_eta, abs_dxy, abs_dz;
};

struct Config {
  bool is_mc = false, is_signal = false;
  std::string mode;
  std::vector<std::string> paths;
  double match_dr = 0.1, muon_merge_dr = 0.05;
  h4l::Bending bending;
  Floor muon{}, electron{};
  double mass_low = 40, mass_high = 140;
  double photon_pt = 10, photon_abs_eta = 2.5, photon_min_dr = 0.3, trigobj_min_dr = 0.3;
  double tag_pt_muon = 20, tag_pt_electron = 25;
  double fsr_pt = 2, fsr_abs_eta = 2.4, fsr_rel_iso = 1.8, fsr_dr_et2 = 0.012, fsr_max_dr = 0.5;
  int multilepton_min = 3;
  double stxs_abs_y = 2.5;
  std::vector<std::string> slim_keep;
  json fiducial;
};

Config read_config(const json& task) {
  const json& analysis = task.at("analysis_config");
  const json& skim = analysis.at("skim");
  const json& triggers = analysis.at("triggers");
  Config config;
  config.is_mc = task.at("config").at("kind") == "mc";
  config.is_signal = task.at("config").value("role", "") == "signal";
  config.mode = task.at("config").value("mode", json()).is_string() ? task.at("config").at("mode").get<std::string>() : "";
  if (config.is_signal && config.mode != "ggH" && config.mode != "VBF" && config.mode != "VH")
    throw std::runtime_error("signal task without a known production mode");
  for (const auto& path : triggers.at("paths")) config.paths.push_back(path.get<std::string>());
  if (config.paths.size() > 32) throw std::runtime_error("more than 32 trigger paths do not fit the trigger mask");
  config.match_dr = triggers.at("trigobj_match_dr").get<double>();
  config.muon_merge_dr = triggers.at("muon_object_merge_dr").get<double>();
  const json& position = triggers.at("electron_object_position");
  config.bending = {position.at("field_tesla").get<double>(), position.at("barrel_radius_m").get<double>(),
                    position.at("endcap_z_m").get<double>()};
  auto floor = [](const json& node, const char* eta) {
    return Floor{node.at("pt").get<double>(), node.at(eta).get<double>(), node.at("abs_dxy").get<double>(),
                 node.at("abs_dz").get<double>()};
  };
  config.muon = floor(skim.at("muon_floor"), "abs_eta");
  config.electron = floor(skim.at("electron_floor"), "abs_eta_sc");
  config.mass_low = skim.at("pair_mass_window").at(0).get<double>();
  config.mass_high = skim.at("pair_mass_window").at(1).get<double>();
  config.photon_pt = skim.at("probe_photon").at("pt").get<double>();
  config.photon_abs_eta = skim.at("probe_photon").at("abs_eta").get<double>();
  config.photon_min_dr = skim.at("probe_photon").at("min_dr_to_tag").get<double>();
  config.trigobj_min_dr = skim.at("probe_trigobj").at("min_dr_to_tag").get<double>();
  config.tag_pt_muon = skim.at("tag_min_pt").at("muon").get<double>();
  config.tag_pt_electron = skim.at("tag_min_pt").at("electron").get<double>();
  const json& fsr = skim.at("fsr_an");
  config.fsr_pt = fsr.at("pt").get<double>();
  config.fsr_abs_eta = fsr.at("abs_eta").get<double>();
  config.fsr_rel_iso = fsr.at("rel_iso").get<double>();
  config.fsr_dr_et2 = fsr.at("dr_over_et2").get<double>();
  config.fsr_max_dr = fsr.at("max_dr").get<double>();
  config.multilepton_min = skim.at("multilepton_min_leptons").get<int>();
  config.stxs_abs_y = analysis.at("stxs").at("abs_y_max").get<double>();
  for (const auto& pattern : skim.at("slim_keep")) config.slim_keep.push_back(pattern.get<std::string>());
  config.fiducial = analysis.at("fiducial");
  return config;
}

// ----------------------------------------------------------------------------
// Input: every branch the skim reads.  MC-only and signal-only branches are
// created only for such samples.
struct Mc {
  Value<Float_t> genWeight;
  Value<Float_t> nTrueInt;
  Array<Float_t> gen_pt, gen_eta, gen_phi, gen_mass;
  Array<Int_t> gen_pdg, gen_status, gen_flags, gen_mother;
  Array<Int_t> mu_gen_idx, el_gen_idx;
  Array<UChar_t> mu_gen_flav, el_gen_flav;
  explicit Mc(TTreeReader& r)
      : genWeight(r, "genWeight"), nTrueInt(r, "Pileup_nTrueInt"), gen_pt(r, "GenPart_pt"), gen_eta(r, "GenPart_eta"),
        gen_phi(r, "GenPart_phi"), gen_mass(r, "GenPart_mass"), gen_pdg(r, "GenPart_pdgId"),
        gen_status(r, "GenPart_status"), gen_flags(r, "GenPart_statusFlags"), gen_mother(r, "GenPart_genPartIdxMother"),
        mu_gen_idx(r, "Muon_genPartIdx"), el_gen_idx(r, "Electron_genPartIdx"), mu_gen_flav(r, "Muon_genPartFlav"),
        el_gen_flav(r, "Electron_genPartFlav") {}
};

struct Signal {
  Array<Float_t> cand_pt, cand_eta, cand_phi, cand_mass;
  Array<Int_t> cand_pdg;
  Array<Float_t> jet_pt, jet_eta, jet_phi, jet_mass;
  explicit Signal(TTreeReader& r)
      : cand_pt(r, "GenCands_pt"), cand_eta(r, "GenCands_eta"), cand_phi(r, "GenCands_phi"), cand_mass(r, "GenCands_mass"),
        cand_pdg(r, "GenCands_pdgId"), jet_pt(r, "GenJet_pt"), jet_eta(r, "GenJet_eta"), jet_phi(r, "GenJet_phi"),
        jet_mass(r, "GenJet_mass") {}
};

#define H4L_F(name, branch) Array<Float_t> name{reader, branch}
#define H4L_I(name, branch) Array<Int_t> name{reader, branch}
#define H4L_B(name, branch) Array<Bool_t> name{reader, branch}
#define H4L_U(name, branch) Array<UChar_t> name{reader, branch}

struct Reader {
  TTreeReader reader;
  // Muons.
  H4L_F(mu_pt, "Muon_pt"); H4L_F(mu_eta, "Muon_eta"); H4L_F(mu_phi, "Muon_phi");
  H4L_F(mu_dxy, "Muon_dxy"); H4L_F(mu_dz, "Muon_dz"); H4L_F(mu_sip, "Muon_sip3d"); H4L_F(mu_ip3d, "Muon_ip3d");
  H4L_F(mu_dxy_err, "Muon_dxyErr"); H4L_F(mu_dz_err, "Muon_dzErr");
  H4L_F(mu_iso03, "Muon_pfRelIso03_all"); H4L_F(mu_iso03_chg, "Muon_pfRelIso03_chg"); H4L_F(mu_iso04, "Muon_pfRelIso04_all");
  H4L_F(mu_mini_iso, "Muon_miniPFRelIso_all"); H4L_F(mu_tk_iso, "Muon_tkRelIso"); H4L_F(mu_pterr, "Muon_ptErr");
  H4L_F(mu_tunep, "Muon_tunepRelPt"); H4L_F(mu_segment, "Muon_segmentComp"); H4L_F(mu_mva_lowpt, "Muon_mvaLowPt");
  H4L_F(mu_soft_mva, "Muon_softMva");
  H4L_I(mu_charge, "Muon_charge"); H4L_I(mu_nstations, "Muon_nStations"); H4L_I(mu_nlayers, "Muon_nTrackerLayers");
  H4L_B(mu_global, "Muon_isGlobal"); H4L_B(mu_tracker, "Muon_isTracker"); H4L_B(mu_standalone, "Muon_isStandalone");
  H4L_B(mu_pf, "Muon_isPFcand"); H4L_B(mu_loose, "Muon_looseId"); H4L_B(mu_medium, "Muon_mediumId");
  H4L_B(mu_tight, "Muon_tightId"); H4L_B(mu_soft, "Muon_softId"); H4L_B(mu_trigloose, "Muon_triggerIdLoose");
  H4L_B(mu_intime, "Muon_inTimeMuon"); H4L_B(mu_medium_prompt, "Muon_mediumPromptId"); H4L_B(mu_high_purity, "Muon_highPurity");
  H4L_B(mu_soft_mva_id, "Muon_softMvaId");
  H4L_U(mu_highpt, "Muon_highPtId");
  // Electrons.
  H4L_F(el_pt, "Electron_pt"); H4L_F(el_eta, "Electron_eta"); H4L_F(el_phi, "Electron_phi");
  H4L_F(el_deta_sc, "Electron_deltaEtaSC"); H4L_F(el_dxy, "Electron_dxy"); H4L_F(el_dz, "Electron_dz");
  H4L_F(el_sip, "Electron_sip3d"); H4L_F(el_ip3d, "Electron_ip3d"); H4L_F(el_dxy_err, "Electron_dxyErr");
  H4L_F(el_dz_err, "Electron_dzErr"); H4L_F(el_iso03, "Electron_pfRelIso03_all"); H4L_F(el_iso03_chg, "Electron_pfRelIso03_chg");
  H4L_F(el_mini_iso, "Electron_miniPFRelIso_all"); H4L_F(el_energy_err, "Electron_energyErr");
  H4L_F(el_mva_noiso, "Electron_mvaFall17V2noIso"); H4L_F(el_mva_iso, "Electron_mvaFall17V2Iso");
  H4L_F(el_sieie, "Electron_sieie"); H4L_F(el_hoe, "Electron_hoe"); H4L_F(el_r9, "Electron_r9");
  H4L_F(el_einv, "Electron_eInvMinusPInv"); H4L_F(el_sc_et, "Electron_scEtOverPt");
  H4L_F(el_dr03_ecal, "Electron_dr03EcalRecHitSumEt"); H4L_F(el_dr03_hcal, "Electron_dr03HcalDepth1TowerSumEt");
  H4L_F(el_dr03_tk, "Electron_dr03TkSumPt"); H4L_F(el_ecorr, "Electron_eCorr");
  H4L_I(el_charge, "Electron_charge"); H4L_I(el_cutbased, "Electron_cutBased"); H4L_I(el_vid, "Electron_vidNestedWPBitmap");
  H4L_I(el_tight_charge, "Electron_tightCharge");
  H4L_B(el_noiso_wpl, "Electron_mvaFall17V2noIso_WPL"); H4L_B(el_noiso_wp90, "Electron_mvaFall17V2noIso_WP90");
  H4L_B(el_noiso_wp80, "Electron_mvaFall17V2noIso_WP80"); H4L_B(el_iso_wpl, "Electron_mvaFall17V2Iso_WPL");
  H4L_B(el_iso_wp90, "Electron_mvaFall17V2Iso_WP90"); H4L_B(el_iso_wp80, "Electron_mvaFall17V2Iso_WP80");
  H4L_B(el_convveto, "Electron_convVeto");
  H4L_U(el_losthits, "Electron_lostHits");
  // FSR photons, photons, trigger objects.
  H4L_F(fsr_pt, "FsrPhoton_pt"); H4L_F(fsr_eta, "FsrPhoton_eta"); H4L_F(fsr_phi, "FsrPhoton_phi");
  H4L_F(fsr_iso, "FsrPhoton_relIso03"); H4L_F(fsr_dret2, "FsrPhoton_dROverEt2"); H4L_I(fsr_muon, "FsrPhoton_muonIdx");
  H4L_F(ph_pt, "Photon_pt"); H4L_F(ph_eta, "Photon_eta"); H4L_F(ph_phi, "Photon_phi"); H4L_F(ph_hoe, "Photon_hoe");
  H4L_F(ph_r9, "Photon_r9"); H4L_F(ph_sieie, "Photon_sieie"); H4L_F(ph_iso, "Photon_pfRelIso03_all");
  H4L_I(ph_el_idx, "Photon_electronIdx"); H4L_B(ph_pixel, "Photon_pixelSeed"); H4L_B(ph_eb, "Photon_isScEtaEB");
  H4L_B(ph_ee, "Photon_isScEtaEE");
  H4L_F(to_pt, "TrigObj_pt"); H4L_F(to_eta, "TrigObj_eta"); H4L_F(to_phi, "TrigObj_phi");
  H4L_I(to_id, "TrigObj_id"); H4L_I(to_bits, "TrigObj_filterBits");
  Value<Int_t> npv{reader, "PV_npvs"}, npv_good{reader, "PV_npvsGood"};
  Value<Float_t> rho{reader, "fixedGridRhoFastjetAll"}, met{reader, "MET_pt"}, met_phi{reader, "MET_phi"};
  std::vector<std::unique_ptr<Value<Bool_t>>> paths;
  std::unique_ptr<Mc> mc;
  std::unique_ptr<Signal> signal;

  Reader(TTree* tree, const Config& config) : reader(tree) {
    for (const auto& path : config.paths) paths.push_back(std::make_unique<Value<Bool_t>>(reader, path.c_str()));
    if (config.is_mc) mc = std::make_unique<Mc>(reader);
    if (config.is_signal) signal = std::make_unique<Signal>(reader);
  }
};
#undef H4L_F
#undef H4L_I
#undef H4L_B
#undef H4L_U

// ----------------------------------------------------------------------------
// Output rows.
struct Leg {
  Short_t index = -1;
  Float_t pt = 0, eta = 0, phi = 0, eta_sc = 0, dxy = 0, dz = 0, sip = 0, ip3d = 0, dxy_err = 0, dz_err = 0, iso03 = 0,
          iso03_chg = 0, iso04 = -1, mini_iso = 0, pt_err = 0, energy_err = -1, mva_noiso = -9, mva_iso = -9;
  // Electron ID inputs (-9 for muons).
  Float_t sieie = -9, hoe = -9, r9 = -9, e_inv_minus_p_inv = -9, sc_et_over_pt = -9, dr03_ecal = -9, dr03_hcal = -9,
          dr03_tk = -9, e_corr = -9;
  Int_t vid_bitmap = 0, tight_charge = -1;
  // Muon ID inputs (-9 for electrons).
  Float_t tunep_rel_pt = -9, segment_comp = -9, tk_rel_iso = -9, mva_lowpt = -9, soft_mva = -9;
  UChar_t n_tracker_layers = 0;
  Int_t charge = 0;
  UShort_t flags = 0;  // muon or electron flag bits, see flag_doc()
  UChar_t cut_based = 0, lost_hits = 0, n_stations = 0;
  Float_t fsr_pt = -1, fsr_eta = 0, fsr_phi = 0, fsr_dret2 = -1, fsr_iso = -1;
  Int_t trig_bits = 0, trig_objects = 0;
  Float_t trig_pt = -1, trig_dr = -1;
  Short_t gen_idx = -1;
  UChar_t gen_flav = 0;
  Float_t gen_pt = -1;

  void book(TTree& tree, const std::string& prefix) {
    auto branch = [&](const char* field, auto* address) { tree.Branch((prefix + field).c_str(), address); };
    branch("index", &index);
    for (auto [field, address] : std::vector<std::pair<const char*, Float_t*>>{
             {"pt", &pt}, {"eta", &eta}, {"phi", &phi}, {"eta_sc", &eta_sc}, {"dxy", &dxy}, {"dz", &dz}, {"sip", &sip},
             {"ip3d", &ip3d}, {"dxy_err", &dxy_err}, {"dz_err", &dz_err}, {"iso03", &iso03}, {"iso03_chg", &iso03_chg},
             {"iso04", &iso04}, {"mini_iso", &mini_iso}, {"pt_err", &pt_err}, {"energy_err", &energy_err},
             {"mva_noiso", &mva_noiso}, {"mva_iso", &mva_iso}, {"sieie", &sieie}, {"hoe", &hoe}, {"r9", &r9},
             {"e_inv_minus_p_inv", &e_inv_minus_p_inv}, {"sc_et_over_pt", &sc_et_over_pt}, {"dr03_ecal", &dr03_ecal},
             {"dr03_hcal", &dr03_hcal}, {"dr03_tk", &dr03_tk}, {"e_corr", &e_corr}, {"tunep_rel_pt", &tunep_rel_pt},
             {"segment_comp", &segment_comp}, {"tk_rel_iso", &tk_rel_iso}, {"mva_lowpt", &mva_lowpt},
             {"soft_mva", &soft_mva}, {"fsr_pt", &fsr_pt}, {"fsr_eta", &fsr_eta}, {"fsr_phi", &fsr_phi},
             {"fsr_dret2", &fsr_dret2}, {"fsr_iso", &fsr_iso}, {"trig_pt", &trig_pt}, {"trig_dr", &trig_dr},
             {"gen_pt", &gen_pt}})
      branch(field, address);
    branch("vid_bitmap", &vid_bitmap);
    branch("tight_charge", &tight_charge);
    branch("n_tracker_layers", &n_tracker_layers);
    branch("charge", &charge);
    branch("flags", &flags);
    branch("cut_based", &cut_based);
    branch("lost_hits", &lost_hits);
    branch("n_stations", &n_stations);
    branch("trig_bits", &trig_bits);
    branch("trig_objects", &trig_objects);
    branch("gen_idx", &gen_idx);
    branch("gen_flav", &gen_flav);
  }
};

json flag_doc() {
  return {{"muon", {"isGlobal", "isTracker", "isStandalone", "isPFcand", "looseId", "mediumId", "tightId", "softId",
                    "highPtId>=1", "highPtId==2", "nStations>0", "triggerIdLoose", "inTimeMuon", "mediumPromptId",
                    "highPurity", "softMvaId"}},
          {"electron", {"mvaFall17V2noIso_WPL", "mvaFall17V2noIso_WP90", "mvaFall17V2noIso_WP80", "mvaFall17V2Iso_WPL",
                        "mvaFall17V2Iso_WP90", "mvaFall17V2Iso_WP80", "convVeto"}}};
}

struct EventFields {
  ULong64_t file_key = 0;
  Long64_t entry = 0;
  Float_t weight = 1, pu_true = -1, rho = 0, met = 0, met_phi = 0;
  Int_t npv = 0, npv_good = 0;
  UInt_t trig = 0;
  UChar_t n_mu = 0, n_el = 0;
  void book(TTree& tree) {
    tree.Branch("file_key", &file_key);
    tree.Branch("entry", &entry);
    tree.Branch("weight", &weight);
    tree.Branch("pu_true", &pu_true);
    tree.Branch("rho", &rho);
    tree.Branch("met", &met);
    tree.Branch("met_phi", &met_phi);
    tree.Branch("npv", &npv);
    tree.Branch("npv_good", &npv_good);
    tree.Branch("trig", &trig);
    tree.Branch("n_mu", &n_mu);
    tree.Branch("n_el", &n_el);
  }
};

// ----------------------------------------------------------------------------
// Per-entry helpers.
struct Skimmer {
  const Config& config;
  Reader& in;
  std::vector<h4l::TrigObject> muon_objects, electron_objects;

  Skimmer(const Config& c, Reader& r) : config(c), in(r) {}

  void prepare_trigger_objects() {
    muon_objects = h4l::trigger_objects(in.to_pt, in.to_eta, in.to_phi, in.to_id, in.to_bits, 13, config.muon_merge_dr);
    electron_objects = h4l::trigger_objects(in.to_pt, in.to_eta, in.to_phi, in.to_id, in.to_bits, 11, 0.0);
  }

  bool muon_floor(std::size_t i) const {
    return in.mu_pt[i] > config.muon.pt && std::fabs(in.mu_eta[i]) < config.muon.abs_eta &&
           std::fabs(in.mu_dxy[i]) < config.muon.abs_dxy && std::fabs(in.mu_dz[i]) < config.muon.abs_dz;
  }
  // The AN loose electron uses abs(eta) < 2.5; keep abs(eta) or abs(eta_SC).
  bool electron_floor(std::size_t i) const {
    const double eta_sc = in.el_eta[i] + in.el_deta_sc[i];
    return in.el_pt[i] > config.electron.pt &&
           (std::fabs(in.el_eta[i]) < config.electron.abs_eta || std::fabs(eta_sc) < config.electron.abs_eta) &&
           std::fabs(in.el_dxy[i]) < config.electron.abs_dxy && std::fabs(in.el_dz[i]) < config.electron.abs_dz;
  }
  P4 muon_p4(std::size_t i) const { return P4(in.mu_pt[i], in.mu_eta[i], in.mu_phi[i], h4l::kMuonMass); }
  P4 electron_p4(std::size_t i) const { return P4(in.el_pt[i], in.el_eta[i], in.el_phi[i], h4l::kElectronMass); }

  // The FsrPhoton of this muon that passes the AN cuts with the lowest dR/ET^2.
  bool an_fsr(std::size_t muon, std::size_t& photon) const {
    bool found = false;
    double best = 0;
    for (std::size_t k = 0; k < in.fsr_pt.GetSize(); ++k) {
      if (in.fsr_muon[k] != static_cast<int>(muon)) continue;
      const double dr = h4l::delta_r(in.mu_eta[muon], in.mu_phi[muon], in.fsr_eta[k], in.fsr_phi[k]);
      if (in.fsr_pt[k] <= config.fsr_pt || std::fabs(in.fsr_eta[k]) >= config.fsr_abs_eta ||
          in.fsr_iso[k] >= config.fsr_rel_iso || in.fsr_dret2[k] >= config.fsr_dr_et2 || dr >= config.fsr_max_dr)
        continue;
      if (!found || in.fsr_dret2[k] < best) {
        best = in.fsr_dret2[k];
        photon = k;
        found = true;
      }
    }
    return found;
  }

  void set_match(Leg& leg, const h4l::TrigMatch& match) const {
    leg.trig_bits = match.bits;
    leg.trig_objects = match.objects;
    leg.trig_pt = match.pt;
    leg.trig_dr = match.dr;
  }

  void fill_muon(Leg& leg, std::size_t i) const {
    leg = Leg{};
    leg.index = static_cast<Short_t>(i);
    leg.pt = in.mu_pt[i];
    leg.eta = in.mu_eta[i];
    leg.phi = in.mu_phi[i];
    leg.eta_sc = in.mu_eta[i];
    leg.dxy = in.mu_dxy[i];
    leg.dz = in.mu_dz[i];
    leg.sip = in.mu_sip[i];
    leg.ip3d = in.mu_ip3d[i];
    leg.dxy_err = in.mu_dxy_err[i];
    leg.dz_err = in.mu_dz_err[i];
    leg.iso03 = in.mu_iso03[i];
    leg.iso03_chg = in.mu_iso03_chg[i];
    leg.iso04 = in.mu_iso04[i];
    leg.mini_iso = in.mu_mini_iso[i];
    leg.pt_err = in.mu_pterr[i];
    leg.tunep_rel_pt = in.mu_tunep[i];
    leg.segment_comp = in.mu_segment[i];
    leg.tk_rel_iso = in.mu_tk_iso[i];
    leg.mva_lowpt = in.mu_mva_lowpt[i];
    leg.soft_mva = in.mu_soft_mva[i];
    leg.n_tracker_layers = static_cast<UChar_t>(std::clamp(in.mu_nlayers[i], 0, 255));
    leg.charge = in.mu_charge[i];
    leg.n_stations = static_cast<UChar_t>(std::clamp(in.mu_nstations[i], 0, 255));
    const bool bits[] = {in.mu_global[i],  in.mu_tracker[i],       in.mu_standalone[i], in.mu_pf[i],
                         in.mu_loose[i],   in.mu_medium[i],        in.mu_tight[i],      in.mu_soft[i],
                         in.mu_highpt[i] >= 1, in.mu_highpt[i] == 2, in.mu_nstations[i] > 0, in.mu_trigloose[i],
                         in.mu_intime[i],  in.mu_medium_prompt[i], in.mu_high_purity[i], in.mu_soft_mva_id[i]};
    for (std::size_t b = 0; b < std::size(bits); ++b) leg.flags |= static_cast<UShort_t>(bits[b]) << b;
    std::size_t photon = 0;
    if (an_fsr(i, photon)) {
      leg.fsr_pt = in.fsr_pt[photon];
      leg.fsr_eta = in.fsr_eta[photon];
      leg.fsr_phi = in.fsr_phi[photon];
      leg.fsr_dret2 = in.fsr_dret2[photon];
      leg.fsr_iso = in.fsr_iso[photon];
    }
    set_match(leg, h4l::match_trigger(muon_objects, leg.eta, leg.phi, config.match_dr));
    if (in.mc) {
      leg.gen_idx = static_cast<Short_t>(in.mc->mu_gen_idx[i]);
      leg.gen_flav = in.mc->mu_gen_flav[i];
      if (leg.gen_idx >= 0 && static_cast<std::size_t>(leg.gen_idx) < in.mc->gen_pt.GetSize()) leg.gen_pt = in.mc->gen_pt[leg.gen_idx];
    }
  }

  void fill_electron(Leg& leg, std::size_t i) const {
    leg = Leg{};
    leg.index = static_cast<Short_t>(i);
    leg.pt = in.el_pt[i];
    leg.eta = in.el_eta[i];
    leg.phi = in.el_phi[i];
    leg.eta_sc = in.el_eta[i] + in.el_deta_sc[i];
    leg.dxy = in.el_dxy[i];
    leg.dz = in.el_dz[i];
    leg.sip = in.el_sip[i];
    leg.ip3d = in.el_ip3d[i];
    leg.dxy_err = in.el_dxy_err[i];
    leg.dz_err = in.el_dz_err[i];
    leg.iso03 = in.el_iso03[i];
    leg.iso03_chg = in.el_iso03_chg[i];
    leg.mini_iso = in.el_mini_iso[i];
    leg.energy_err = in.el_energy_err[i];
    leg.pt_err = in.el_energy_err[i] / std::cosh(in.el_eta[i]);
    leg.mva_noiso = in.el_mva_noiso[i];
    leg.mva_iso = in.el_mva_iso[i];
    leg.sieie = in.el_sieie[i];
    leg.hoe = in.el_hoe[i];
    leg.r9 = in.el_r9[i];
    leg.e_inv_minus_p_inv = in.el_einv[i];
    leg.sc_et_over_pt = in.el_sc_et[i];
    leg.dr03_ecal = in.el_dr03_ecal[i];
    leg.dr03_hcal = in.el_dr03_hcal[i];
    leg.dr03_tk = in.el_dr03_tk[i];
    leg.e_corr = in.el_ecorr[i];
    leg.vid_bitmap = in.el_vid[i];
    leg.tight_charge = in.el_tight_charge[i];
    leg.charge = in.el_charge[i];
    leg.cut_based = static_cast<UChar_t>(std::clamp(in.el_cutbased[i], 0, 255));
    leg.lost_hits = in.el_losthits[i];
    const bool bits[] = {in.el_noiso_wpl[i], in.el_noiso_wp90[i], in.el_noiso_wp80[i], in.el_iso_wpl[i],
                         in.el_iso_wp90[i],  in.el_iso_wp80[i],   in.el_convveto[i]};
    for (std::size_t b = 0; b < std::size(bits); ++b) leg.flags |= static_cast<UShort_t>(bits[b]) << b;
    const auto [eta_hlt, phi_hlt] = h4l::electron_hlt_position(leg.eta_sc, leg.phi, leg.charge, leg.pt, config.bending);
    set_match(leg, h4l::match_trigger(electron_objects, eta_hlt, phi_hlt, config.match_dr));
    if (in.mc) {
      leg.gen_idx = static_cast<Short_t>(in.mc->el_gen_idx[i]);
      leg.gen_flav = in.mc->el_gen_flav[i];
      if (leg.gen_idx >= 0 && static_cast<std::size_t>(leg.gen_idx) < in.mc->gen_pt.GetSize()) leg.gen_pt = in.mc->gen_pt[leg.gen_idx];
    }
  }

  P4 leg_p4_with_fsr(const Leg& leg, double mass) const {
    P4 p4(leg.pt, leg.eta, leg.phi, mass);
    if (leg.fsr_pt > 0) p4 += P4(leg.fsr_pt, leg.fsr_eta, leg.fsr_phi, 0.0);
    return p4;
  }
};

// ----------------------------------------------------------------------------
// Generator table of the signal samples.
constexpr int kMaxDressed = 8;

struct GenRow {
  Float_t weight = 1;
  Float_t h_pt = -1, h_eta = 0, h_phi = 0, h_mass = -1, h_y = 0;
  Int_t n_hlep = 0, final_state_true = -1, stage0 = 0, vh_class = 0;
  std::array<Float_t, 4> hlep_pt{}, hlep_eta{}, hlep_phi{};
  std::array<Int_t, 4> hlep_pdg{};
  Float_t m4l_bare = -1, m4l_dressed = -1;
  Int_t n_dressed = 0, n_selected = 0;
  std::array<Float_t, kMaxDressed> d_pt{}, d_eta{}, d_phi{}, d_mass{}, d_iso{};
  std::array<Int_t, kMaxDressed> d_pdg{};
  std::array<Bool_t, kMaxDressed> d_from_h{}, d_selected{};
  Bool_t fid_pass = false;
  std::array<Int_t, 4> fid_idx{-1, -1, -1, -1};
  Float_t fid_m4l = -1, fid_mz1 = -1, fid_mz2 = -1, fid_pt4l = -1, fid_y4l = 0;
  Int_t fid_final_state = -1, fid_njets = 0;
  Float_t fid_jet1_pt = -1;
};

struct Dressed {
  P4 p4;
  int pdg = 0;
  bool from_h = false;
  double iso = 0;
};

int final_state_code(int n_electrons, int n_muons) {
  if (n_muons == 4) return 0;                      // 4mu
  if (n_electrons == 4) return 1;                  // 4e
  if (n_muons == 2 && n_electrons == 2) return 2;  // 2e2mu
  return -1;
}

bool has_ancestor(const Mc& mc, int index, int pdg) {
  int mother = mc.gen_mother[index];
  for (int guard = 0; mother >= 0 && guard < 200; ++guard) {
    if (mc.gen_pdg[mother] == pdg) return true;
    mother = mc.gen_mother[mother];
  }
  return false;
}

// The V of a VH event: hard-process W/Z that does not descend from the H.
// Its decay products (first generation after the last copy) classify the
// event: 1 W->l nu, 2 Z->ll, 3 Z->nu nu, 4 W->qq, 5 Z->qq, 0 unknown.
int vh_decay_class(const Mc& mc) {
  const std::size_t n = mc.gen_pt.GetSize();
  for (std::size_t i = 0; i < n; ++i) {
    const int pdg = std::abs(mc.gen_pdg[i]);
    if ((pdg != 23 && pdg != 24) || !(mc.gen_flags[i] >> 7 & 1) || has_ancestor(mc, static_cast<int>(i), 25)) continue;
    // Walk to the last copy of this V.
    int last = static_cast<int>(i);
    for (bool moved = true; moved;) {
      moved = false;
      for (std::size_t j = 0; j < n; ++j)
        if (mc.gen_mother[j] == last && mc.gen_pdg[j] == mc.gen_pdg[i]) {
          last = static_cast<int>(j);
          moved = true;
          break;
        }
    }
    bool quark = false, charged = false, neutrino = false;
    for (std::size_t j = 0; j < n; ++j) {
      if (mc.gen_mother[j] != last) continue;
      const int d = std::abs(mc.gen_pdg[j]);
      if (d >= 1 && d <= 5) quark = true;
      if (d == 11 || d == 13 || d == 15) charged = true;
      if (d == 12 || d == 14 || d == 16) neutrino = true;
    }
    if (pdg == 24) return quark ? 4 : (charged ? 1 : 0);
    return quark ? 5 : (charged ? 2 : (neutrino ? 3 : 0));
  }
  return 0;
}

void fill_gen(GenRow& row, Reader& in, const Config& config) {
  row = GenRow{};
  const json& fid = config.fiducial;
  Mc& mc = *in.mc;
  Signal& signal = *in.signal;
  row.weight = *mc.genWeight;
  const std::size_t n = mc.gen_pt.GetSize();
  // The Higgs: the last copy with pdgId 25.
  for (std::size_t i = 0; i < n; ++i) {
    if (mc.gen_pdg[i] == 25 && (mc.gen_flags[i] >> 13 & 1)) {
      const P4 h(mc.gen_pt[i], mc.gen_eta[i], mc.gen_phi[i], mc.gen_mass[i]);
      row.h_pt = h.Pt();
      row.h_eta = h.Eta();
      row.h_phi = h.Phi();
      row.h_mass = h.M();
      row.h_y = h.Rapidity();
    }
  }
  // STXS stage 0 from the production mode and the generator record.
  const bool forward = std::fabs(row.h_y) > config.stxs_abs_y;
  if (config.mode == "ggH") row.stage0 = forward ? 10 : 11;
  if (config.mode == "VBF") row.stage0 = forward ? 20 : 21;
  if (config.mode == "VH") {
    row.vh_class = vh_decay_class(mc);
    if (row.vh_class == 4 || row.vh_class == 5) row.stage0 = forward ? 22 : 23;
    if (row.vh_class == 1) row.stage0 = forward ? 30 : 31;
    if (row.vh_class == 2 || row.vh_class == 3) row.stage0 = forward ? 40 : 41;
  }

  // Prompt stable electrons and muons; the H-decay leptons also come from
  // the hard process (this excludes gamma* -> ll in the shower).
  std::vector<std::size_t> prompt;
  std::vector<bool> from_h;
  for (std::size_t i = 0; i < n; ++i) {
    const int pdg = std::abs(mc.gen_pdg[i]);
    if ((pdg != 11 && pdg != 13) || mc.gen_status[i] != 1 || !(mc.gen_flags[i] & 1)) continue;
    prompt.push_back(i);
    from_h.push_back((mc.gen_flags[i] >> 8 & 1) && has_ancestor(mc, static_cast<int>(i), 25));
  }
  int ne = 0, nm = 0;
  P4 bare;
  for (std::size_t k = 0; k < prompt.size(); ++k) {
    if (!from_h[k]) continue;
    const std::size_t i = prompt[k];
    (std::abs(mc.gen_pdg[i]) == 11 ? ne : nm) += 1;
    bare += P4(mc.gen_pt[i], mc.gen_eta[i], mc.gen_phi[i], mc.gen_mass[i]);
    if (row.n_hlep < 4) {
      row.hlep_pt[row.n_hlep] = mc.gen_pt[i];
      row.hlep_eta[row.n_hlep] = mc.gen_eta[i];
      row.hlep_phi[row.n_hlep] = mc.gen_phi[i];
      row.hlep_pdg[row.n_hlep] = mc.gen_pdg[i];
    }
    ++row.n_hlep;
  }
  if (row.n_hlep == 4) {
    row.final_state_true = final_state_code(ne, nm);
    row.m4l_bare = bare.M();
  }

  // Dressing: every GenCands photon goes to the closest prompt lepton within
  // the dressing cone.
  const double dress_dr = fid.at("dressing_dr").get<double>();
  std::vector<Dressed> leptons(prompt.size());
  std::vector<P4> bare_p4(prompt.size());
  for (std::size_t k = 0; k < prompt.size(); ++k) {
    const std::size_t i = prompt[k];
    bare_p4[k] = P4(mc.gen_pt[i], mc.gen_eta[i], mc.gen_phi[i], mc.gen_mass[i]);
    leptons[k] = {bare_p4[k], mc.gen_pdg[i], from_h[k], 0.0};
  }
  const std::size_t n_cands = signal.cand_pt.GetSize();
  std::vector<bool> dressing(n_cands, false);
  for (std::size_t c = 0; c < n_cands; ++c) {
    if (signal.cand_pdg[c] != 22) continue;
    double best = dress_dr;
    int owner = -1;
    for (std::size_t k = 0; k < bare_p4.size(); ++k) {
      const double dr = h4l::delta_r(bare_p4[k].Eta(), bare_p4[k].Phi(), signal.cand_eta[c], signal.cand_phi[c]);
      if (dr < best) {
        best = dr;
        owner = static_cast<int>(k);
      }
    }
    if (owner >= 0) {
      leptons[owner].p4 += P4(signal.cand_pt[c], signal.cand_eta[c], signal.cand_phi[c], 0.0);
      dressing[c] = true;
    }
  }
  if (row.n_hlep == 4) {
    P4 dressed_sum;
    for (const auto& lepton : leptons)
      if (lepton.from_h) dressed_sum += lepton.p4;
    row.m4l_dressed = dressed_sum.M();
  }

  // Kinematics and isolation (stable particles except e, mu, neutrinos and
  // the dressing photons, within the isolation cone).
  const double iso_dr = fid.at("isolation_dr").get<double>();
  const double iso_max = fid.at("isolation_max_rel").get<double>();
  std::vector<std::size_t> selected;
  for (std::size_t k = 0; k < leptons.size(); ++k) {
    double sum = 0;
    for (std::size_t c = 0; c < n_cands; ++c) {
      const int pdg = std::abs(signal.cand_pdg[c]);
      if (pdg == 11 || pdg == 13 || pdg == 12 || pdg == 14 || pdg == 16 || dressing[c]) continue;
      if (h4l::delta_r(leptons[k].p4.Eta(), leptons[k].p4.Phi(), signal.cand_eta[c], signal.cand_phi[c]) < iso_dr)
        sum += signal.cand_pt[c];
    }
    leptons[k].iso = leptons[k].p4.Pt() > 0 ? sum / leptons[k].p4.Pt() : 1e9;
    const bool electron = std::abs(leptons[k].pdg) == 11;
    const json& cuts = fid.at(electron ? "electron" : "muon");
    if (leptons[k].p4.Pt() > cuts.at("pt").get<double>() && std::fabs(leptons[k].p4.Eta()) < cuts.at("abs_eta").get<double>() &&
        leptons[k].iso < iso_max)
      selected.push_back(k);
  }
  row.n_dressed = static_cast<Int_t>(leptons.size());
  row.n_selected = static_cast<Int_t>(selected.size());
  for (std::size_t k = 0; k < leptons.size() && k < static_cast<std::size_t>(kMaxDressed); ++k) {
    row.d_pt[k] = leptons[k].p4.Pt();
    row.d_eta[k] = leptons[k].p4.Eta();
    row.d_phi[k] = leptons[k].p4.Phi();
    row.d_mass[k] = leptons[k].p4.M();
    row.d_iso[k] = leptons[k].iso;
    row.d_pdg[k] = leptons[k].pdg;
    row.d_from_h[k] = leptons[k].from_h;
    row.d_selected[k] = std::find(selected.begin(), selected.end(), k) != selected.end();
  }

  // Candidate: Z1 closest to m_Z, then the Z2 with the highest sum pT that
  // satisfies every requirement together with Z1.
  const double zmass = fid.at("z_mass").get<double>();
  const double z1_low = fid.at("z1_mass").at(0).get<double>(), z1_high = fid.at("z1_mass").at(1).get<double>();
  const double z2_low = fid.at("z2_mass").at(0).get<double>(), z2_high = fid.at("z2_mass").at(1).get<double>();
  const double lead = fid.at("leading_pt").get<double>(), sublead = fid.at("subleading_pt").get<double>();
  const double min_dr = fid.at("min_dr_leptons").get<double>(), min_os = fid.at("min_os_mass").get<double>();
  const double m4l_low = fid.at("m4l").at(0).get<double>(), m4l_high = fid.at("m4l").at(1).get<double>();
  int z1a = -1, z1b = -1;
  double z1_distance = 1e9;
  for (std::size_t a = 0; a < selected.size(); ++a)
    for (std::size_t b = a + 1; b < selected.size(); ++b) {
      const Dressed& x = leptons[selected[a]];
      const Dressed& y = leptons[selected[b]];
      if (x.pdg + y.pdg != 0) continue;  // opposite-sign same-flavour
      const double m = (x.p4 + y.p4).M();
      if (m <= z1_low || m >= z1_high) continue;
      if (std::fabs(m - zmass) < z1_distance) {
        z1_distance = std::fabs(m - zmass);
        z1a = static_cast<int>(selected[a]);
        z1b = static_cast<int>(selected[b]);
      }
    }
  if (z1a < 0) return;
  int z2a = -1, z2b = -1;
  double best_sum = -1;
  for (std::size_t a = 0; a < selected.size(); ++a)
    for (std::size_t b = a + 1; b < selected.size(); ++b) {
      const int i = static_cast<int>(selected[a]), j = static_cast<int>(selected[b]);
      if (i == z1a || i == z1b || j == z1a || j == z1b) continue;
      const Dressed& x = leptons[i];
      const Dressed& y = leptons[j];
      if (x.pdg + y.pdg != 0) continue;
      const double m = (x.p4 + y.p4).M();
      if (m <= z2_low || m >= z2_high) continue;
      const std::array<int, 4> four{z1a, z1b, i, j};
      std::vector<double> pts;
      bool ok = true;
      for (int p = 0; p < 4 && ok; ++p) {
        pts.push_back(leptons[four[p]].p4.Pt());
        for (int q = p + 1; q < 4 && ok; ++q) {
          const Dressed& u = leptons[four[p]];
          const Dressed& v = leptons[four[q]];
          if (h4l::delta_r(u.p4.Eta(), u.p4.Phi(), v.p4.Eta(), v.p4.Phi()) <= min_dr) ok = false;
          if ((u.pdg > 0) != (v.pdg > 0) && (u.p4 + v.p4).M() <= min_os) ok = false;  // every OS pair
        }
      }
      std::sort(pts.rbegin(), pts.rend());
      if (!ok || pts[0] <= lead || pts[1] <= sublead) continue;
      const double sum = x.p4.Pt() + y.p4.Pt();
      if (sum > best_sum) {
        best_sum = sum;
        z2a = i;
        z2b = j;
      }
    }
  if (z2a < 0) return;
  const P4 z1 = leptons[z1a].p4 + leptons[z1b].p4;
  const P4 z2 = leptons[z2a].p4 + leptons[z2b].p4;
  const P4 zz = z1 + z2;
  row.fid_idx = {z1a, z1b, z2a, z2b};
  row.fid_m4l = zz.M();
  row.fid_mz1 = z1.M();
  row.fid_mz2 = z2.M();
  row.fid_pt4l = zz.Pt();
  row.fid_y4l = zz.Rapidity();
  const int electrons = (std::abs(leptons[z1a].pdg) == 11) * 2 + (std::abs(leptons[z2a].pdg) == 11) * 2;
  row.fid_final_state = final_state_code(electrons, 4 - electrons);
  row.fid_pass = row.fid_m4l > m4l_low && row.fid_m4l < m4l_high;
  // Generator jets, cleaned from the selected dressed leptons.
  const json& jets = fid.at("jets");
  for (std::size_t j = 0; j < signal.jet_pt.GetSize(); ++j) {
    if (signal.jet_pt[j] <= jets.at("pt").get<double>() || std::fabs(signal.jet_eta[j]) >= jets.at("abs_eta").get<double>())
      continue;
    bool clean = true;
    for (std::size_t k : selected)
      if (h4l::delta_r(leptons[k].p4.Eta(), leptons[k].p4.Phi(), signal.jet_eta[j], signal.jet_phi[j]) <
          jets.at("min_dr_leptons").get<double>())
        clean = false;
    if (!clean) continue;
    ++row.fid_njets;
    row.fid_jet1_pt = std::max<Float_t>(row.fid_jet1_pt, signal.jet_pt[j]);
  }
}

void book_gen(TTree& tree, GenRow& row, EventFields& event) {
  tree.Branch("file_key", &event.file_key);
  tree.Branch("entry", &event.entry);
  tree.Branch("weight", &row.weight);
  for (auto [field, address] : std::vector<std::pair<const char*, Float_t*>>{
           {"h_pt", &row.h_pt}, {"h_eta", &row.h_eta}, {"h_phi", &row.h_phi}, {"h_mass", &row.h_mass}, {"h_y", &row.h_y},
           {"m4l_bare", &row.m4l_bare}, {"m4l_dressed", &row.m4l_dressed}, {"fid_m4l", &row.fid_m4l},
           {"fid_mz1", &row.fid_mz1}, {"fid_mz2", &row.fid_mz2}, {"fid_pt4l", &row.fid_pt4l}, {"fid_y4l", &row.fid_y4l},
           {"fid_jet1_pt", &row.fid_jet1_pt}})
    tree.Branch(field, address);
  for (auto [field, address] : std::vector<std::pair<const char*, Int_t*>>{
           {"n_hlep", &row.n_hlep}, {"final_state_true", &row.final_state_true}, {"stage0", &row.stage0},
           {"vh_class", &row.vh_class}, {"n_dressed", &row.n_dressed}, {"n_selected", &row.n_selected},
           {"fid_final_state", &row.fid_final_state}, {"fid_njets", &row.fid_njets}})
    tree.Branch(field, address);
  tree.Branch("fid_pass", &row.fid_pass);
  tree.Branch("fid_idx", row.fid_idx.data(), "fid_idx[4]/I");
  tree.Branch("hlep_pt", row.hlep_pt.data(), "hlep_pt[4]/F");
  tree.Branch("hlep_eta", row.hlep_eta.data(), "hlep_eta[4]/F");
  tree.Branch("hlep_phi", row.hlep_phi.data(), "hlep_phi[4]/F");
  tree.Branch("hlep_pdg", row.hlep_pdg.data(), "hlep_pdg[4]/I");
  const std::string m = std::to_string(kMaxDressed);
  tree.Branch("d_pt", row.d_pt.data(), ("d_pt[" + m + "]/F").c_str());
  tree.Branch("d_eta", row.d_eta.data(), ("d_eta[" + m + "]/F").c_str());
  tree.Branch("d_phi", row.d_phi.data(), ("d_phi[" + m + "]/F").c_str());
  tree.Branch("d_mass", row.d_mass.data(), ("d_mass[" + m + "]/F").c_str());
  tree.Branch("d_iso", row.d_iso.data(), ("d_iso[" + m + "]/F").c_str());
  tree.Branch("d_pdg", row.d_pdg.data(), ("d_pdg[" + m + "]/I").c_str());
  tree.Branch("d_from_h", row.d_from_h.data(), ("d_from_h[" + m + "]/O").c_str());
  tree.Branch("d_selected", row.d_selected.data(), ("d_selected[" + m + "]/O").c_str());
}

// ----------------------------------------------------------------------------
struct FileCounters {
  std::string file_key, path;
  long long entries = 0, pairs = 0, photon_pairs = 0, trigobj_pairs = 0, slim = 0;
  double sum_weight = 0;
};

std::string hex64(std::uint64_t value) {
  char buffer[17];
  std::snprintf(buffer, sizeof(buffer), "%016llx", static_cast<unsigned long long>(value));
  return buffer;
}

}  // namespace

int main(int argc, char** argv) {
  std::string task_path, out_json, out_root;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--task" && index + 1 < argc) task_path = argv[++index];
    else if (argument == "--out-json" && index + 1 < argc) out_json = argv[++index];
    else if (argument == "--out-root" && index + 1 < argc) out_root = argv[++index];
    else {
      std::cerr << "usage: skim_v3 --task TASK.json --out-json OUT.json --out-root OUT.root\n";
      return 2;
    }
  }
  if (task_path.empty() || out_json.empty() || out_root.empty()) {
    std::cerr << "usage: skim_v3 --task TASK.json --out-json OUT.json --out-root OUT.root\n";
    return 2;
  }
  gROOT->SetBatch(true);
  SetErrorHandler(record_errors);
  try {
    const json task = h4l::read_json(task_path);
    const Config config = read_config(task);
    const json& inputs = task.at("inputs");

    auto output = std::make_unique<TFile>(out_root.c_str(), "RECREATE", "", 505);
    if (!output || output->IsZombie()) throw std::runtime_error("cannot create " + out_root);
    // The trees belong to the output file, which deletes them when it closes.
    output->cd();
    TTree& pairs = *new TTree("Pairs", "opposite-sign same-flavour lepton pairs");
    TTree& photon_pairs = *new TTree("PhotonPairs", "tag electron + Photon probe pairs");
    TTree& trigobj_pairs = *new TTree("TrigObjPairs", "tag muon + merged HLT muon-object probe pairs");
    TTree& gen_table = *new TTree("GenTable", "signal generator table");
    TTree& files_tree = *new TTree("Files", "per-input counters");

    std::vector<FileCounters> counters;
    std::vector<std::vector<Long64_t>> slim_entries;
    std::vector<std::string> paths;
    std::vector<ULong64_t> keys;

    // Branch buffers shared by every file.
    EventFields event;
    Int_t pair_flavour = 0;
    Float_t pair_mass = 0, pair_mass_fsr = 0, pair_pt = 0, pair_y = 0;
    Leg leg1, leg2, tag, probe_lepton;
    Float_t probe_pt = 0, probe_eta = 0, probe_phi = 0, probe_hoe = 0, probe_r9 = 0, probe_sieie = 0, probe_iso = 0,
            probe_mass = 0, probe_match_dr = -1;
    Int_t probe_el_idx = -1, probe_bits = 0, probe_matched = 0, probe_members = 0;
    UChar_t probe_flags = 0;
    GenRow gen_row;
    json skipped = json::array();

    event.book(pairs);
    pairs.Branch("flavour", &pair_flavour);
    pairs.Branch("mass", &pair_mass);
    pairs.Branch("mass_fsr", &pair_mass_fsr);
    pairs.Branch("pair_pt", &pair_pt);
    pairs.Branch("pair_y", &pair_y);
    leg1.book(pairs, "l1_");
    leg2.book(pairs, "l2_");

    event.book(photon_pairs);
    tag.book(photon_pairs, "t_");
    for (auto [field, address] : std::vector<std::pair<const char*, Float_t*>>{
             {"p_pt", &probe_pt}, {"p_eta", &probe_eta}, {"p_phi", &probe_phi}, {"p_hoe", &probe_hoe}, {"p_r9", &probe_r9},
             {"p_sieie", &probe_sieie}, {"p_iso", &probe_iso}, {"mass", &probe_mass}, {"p_match_dr", &probe_match_dr}})
      photon_pairs.Branch(field, address);
    photon_pairs.Branch("p_electron_idx", &probe_el_idx);
    photon_pairs.Branch("p_flags", &probe_flags);  // bit 0 pixelSeed, 1 isScEtaEB, 2 isScEtaEE
    photon_pairs.Branch("p_matched", &probe_matched);
    probe_lepton.book(photon_pairs, "e_");

    event.book(trigobj_pairs);
    tag.book(trigobj_pairs, "t_");
    trigobj_pairs.Branch("p_pt", &probe_pt);
    trigobj_pairs.Branch("p_eta", &probe_eta);
    trigobj_pairs.Branch("p_phi", &probe_phi);
    trigobj_pairs.Branch("p_bits", &probe_bits);
    trigobj_pairs.Branch("p_members", &probe_members);
    trigobj_pairs.Branch("mass", &probe_mass);
    trigobj_pairs.Branch("p_match_dr", &probe_match_dr);
    trigobj_pairs.Branch("p_matched", &probe_matched);
    probe_lepton.book(trigobj_pairs, "m_");

    if (config.is_signal) book_gen(gen_table, gen_row, event);
    long long stage0_unknown = 0;

    for (const auto& input : inputs) {
      const std::string path = input.at("path").get<std::string>();
      const std::string key_hex = input.at("file_key").get<std::string>();
      const ULong64_t key = h4l::parse_file_key(key_hex);
      g_read_error = false;
      std::string open_error;
      auto file = h4l::open_input(path, 5, open_error);
      if (!file) {
        // An unreadable MC file is recorded and skipped (its coverage is
        // lost); a pseudo-data shard can never be skipped.
        if (!config.is_mc) throw InputError("cannot open " + path + ": " + open_error);
        skipped.push_back({{"file_key", key_hex}, {"path", path}, {"reason", open_error}});
        std::cerr << "WARN: skipped unreadable MC input " << path << ": " << open_error << "\n";
        g_read_error = false;
        continue;
      }
      auto* tree = dynamic_cast<TTree*>(file->Get("Events"));
      if (!tree) throw InputError("no Events in " + path);
      const long long expected_entries = input.value("events_entries", -1LL);
      if (expected_entries >= 0 && tree->GetEntries() != expected_entries)
        throw InputError(path + " has " + std::to_string(tree->GetEntries()) + " entries, the manifest " +
                         std::to_string(expected_entries));
      Reader in(tree, config);
      Skimmer skim(config, in);
      FileCounters counter;
      counter.file_key = key_hex;
      counter.path = path;
      std::vector<Long64_t> selected;
      const Long64_t entries = tree->GetEntries();
      for (Long64_t entry = 0; entry < entries; ++entry) {
        const auto status = in.reader.SetEntry(entry);
        if (status != TTreeReader::kEntryValid) throw InputError("read error at entry " + std::to_string(entry) + " of " + path);
        event = EventFields{};
        event.file_key = key;
        event.entry = entry;
        event.weight = in.mc ? *in.mc->genWeight : 1.f;
        event.pu_true = in.mc ? *in.mc->nTrueInt : -1.f;
        event.rho = *in.rho;
        event.met = *in.met;
        event.met_phi = *in.met_phi;
        event.npv = *in.npv;
        event.npv_good = *in.npv_good;
        for (std::size_t p = 0; p < in.paths.size(); ++p) event.trig |= static_cast<UInt_t>(**in.paths[p]) << p;
        counter.sum_weight += event.weight;
        skim.prepare_trigger_objects();

        std::vector<std::size_t> muons, electrons;
        for (std::size_t i = 0; i < in.mu_pt.GetSize(); ++i)
          if (skim.muon_floor(i)) muons.push_back(i);
        for (std::size_t i = 0; i < in.el_pt.GetSize(); ++i)
          if (skim.electron_floor(i)) electrons.push_back(i);
        event.n_mu = static_cast<UChar_t>(std::min<std::size_t>(muons.size(), 255));
        event.n_el = static_cast<UChar_t>(std::min<std::size_t>(electrons.size(), 255));

        // Pairs.
        for (int flavour : {13, 11}) {
          const auto& list = flavour == 13 ? muons : electrons;
          for (std::size_t a = 0; a < list.size(); ++a)
            for (std::size_t b = a + 1; b < list.size(); ++b) {
              const std::size_t i = list[a], j = list[b];
              const bool muon = flavour == 13;
              const int qi = muon ? in.mu_charge[i] : in.el_charge[i];
              const int qj = muon ? in.mu_charge[j] : in.el_charge[j];
              if (qi * qj >= 0) continue;
              const P4 pi = muon ? skim.muon_p4(i) : skim.electron_p4(i);
              const P4 pj = muon ? skim.muon_p4(j) : skim.electron_p4(j);
              const P4 sum = pi + pj;
              if (sum.M() <= config.mass_low || sum.M() >= config.mass_high) continue;
              const bool first_leads = pi.Pt() >= pj.Pt();
              if (muon) {
                skim.fill_muon(leg1, first_leads ? i : j);
                skim.fill_muon(leg2, first_leads ? j : i);
              } else {
                skim.fill_electron(leg1, first_leads ? i : j);
                skim.fill_electron(leg2, first_leads ? j : i);
              }
              const double mass = muon ? h4l::kMuonMass : h4l::kElectronMass;
              pair_flavour = flavour;
              pair_mass = sum.M();
              pair_mass_fsr = (skim.leg_p4_with_fsr(leg1, mass) + skim.leg_p4_with_fsr(leg2, mass)).M();
              pair_pt = sum.Pt();
              pair_y = sum.Rapidity();
              if (pairs.Fill() <= 0) throw std::runtime_error("Pairs Fill failed");
              ++counter.pairs;
            }
        }

        // Electron reconstruction probes: tag electron + Photon.
        for (std::size_t t : electrons) {
          if (in.el_pt[t] < config.tag_pt_electron) continue;
          const P4 tp = skim.electron_p4(t);
          for (std::size_t k = 0; k < in.ph_pt.GetSize(); ++k) {
            if (in.ph_pt[k] <= config.photon_pt || std::fabs(in.ph_eta[k]) >= config.photon_abs_eta) continue;
            if (h4l::delta_r(in.el_eta[t], in.el_phi[t], in.ph_eta[k], in.ph_phi[k]) <= config.photon_min_dr) continue;
            const P4 pp(in.ph_pt[k], in.ph_eta[k], in.ph_phi[k], 0.0);
            const double m = (tp + pp).M();
            if (m <= config.mass_low || m >= config.mass_high) continue;
            skim.fill_electron(tag, t);
            probe_pt = in.ph_pt[k];
            probe_eta = in.ph_eta[k];
            probe_phi = in.ph_phi[k];
            probe_hoe = in.ph_hoe[k];
            probe_r9 = in.ph_r9[k];
            probe_sieie = in.ph_sieie[k];
            probe_iso = in.ph_iso[k];
            probe_mass = m;
            probe_el_idx = in.ph_el_idx[k];
            probe_flags = static_cast<UChar_t>(in.ph_pixel[k] | (in.ph_eb[k] << 1) | (in.ph_ee[k] << 2));
            probe_matched = 0;
            probe_match_dr = -1;
            probe_lepton = Leg{};
            if (probe_el_idx >= 0 && static_cast<std::size_t>(probe_el_idx) < in.el_pt.GetSize()) {
              skim.fill_electron(probe_lepton, static_cast<std::size_t>(probe_el_idx));
              probe_matched = 1;
              probe_match_dr = static_cast<Float_t>(h4l::delta_r(in.ph_eta[k], in.ph_phi[k], in.el_eta[probe_el_idx], in.el_phi[probe_el_idx]));
            }
            if (photon_pairs.Fill() <= 0) throw std::runtime_error("PhotonPairs Fill failed");
            ++counter.photon_pairs;
          }
        }

        // Muon reconstruction probes: tag muon + merged HLT muon object.
        for (std::size_t t : muons) {
          if (in.mu_pt[t] < config.tag_pt_muon) continue;
          const P4 tp = skim.muon_p4(t);
          for (const auto& object : skim.muon_objects) {
            if (h4l::delta_r(in.mu_eta[t], in.mu_phi[t], object.eta, object.phi) <= config.trigobj_min_dr) continue;
            const P4 pp(object.pt, object.eta, object.phi, h4l::kMuonMass);
            const double m = (tp + pp).M();
            if (m <= config.mass_low || m >= config.mass_high) continue;
            skim.fill_muon(tag, t);
            probe_pt = object.pt;
            probe_eta = object.eta;
            probe_phi = object.phi;
            probe_bits = object.bits;
            probe_members = object.members;
            probe_mass = m;
            probe_matched = 0;
            probe_match_dr = -1;
            probe_lepton = Leg{};
            double best = config.match_dr;
            int best_index = -1;
            for (std::size_t i = 0; i < in.mu_pt.GetSize(); ++i) {
              const double dr = h4l::delta_r(in.mu_eta[i], in.mu_phi[i], object.eta, object.phi);
              if (dr < best) {
                best = dr;
                best_index = static_cast<int>(i);
              }
            }
            if (best_index >= 0) {
              skim.fill_muon(probe_lepton, static_cast<std::size_t>(best_index));
              probe_matched = 1;
              probe_match_dr = static_cast<Float_t>(best);
            }
            if (trigobj_pairs.Fill() <= 0) throw std::runtime_error("TrigObjPairs Fill failed");
            ++counter.trigobj_pairs;
          }
        }

        if (static_cast<int>(muons.size() + electrons.size()) >= config.multilepton_min) selected.push_back(entry);
        if (config.is_signal) {
          fill_gen(gen_row, in, config);
          stage0_unknown += gen_row.stage0 == 0;
          if (gen_table.Fill() <= 0) throw std::runtime_error("GenTable Fill failed");
        }
        check_read(path + " entry " + std::to_string(entry));
      }
      counter.entries = entries;
      counter.slim = static_cast<long long>(selected.size());
      counters.push_back(counter);
      slim_entries.push_back(selected);
      paths.push_back(path);
      keys.push_back(key);
      std::cout << "[skim] " << path << ": " << entries << " entries, " << counter.pairs << " pairs, " << selected.size()
                << " slim\n";
    }

    // Slim NanoAOD: an explicit copy of the configured branches of the
    // selected entries, with the original file key and entry.
    TChain chain("Events");
    for (const auto& path : paths) chain.Add(path.c_str());
    chain.SetBranchStatus("*", false);
    std::vector<std::string> patterns = config.slim_keep;
    patterns.insert(patterns.end(), config.paths.begin(), config.paths.end());
    std::set<std::string> kept;
    if (chain.LoadTree(0) < 0) throw InputError("cannot load the first slim input");
    TObjArray* branches = chain.GetTree()->GetListOfBranches();
    for (int b = 0; b < branches->GetEntriesFast(); ++b) {
      const std::string name = branches->At(b)->GetName();
      for (const auto& pattern : patterns) {
        if (fnmatch(pattern.c_str(), name.c_str(), 0) == 0) {
          chain.SetBranchStatus(name.c_str(), true);
          kept.insert(name);
          break;
        }
      }
    }
    output->cd();
    TTree* slim = chain.CloneTree(0);
    if (!slim) throw std::runtime_error("CloneTree failed");
    ULong64_t slim_key = 0;
    Long64_t slim_entry = 0;
    slim->Branch("h4l_file_key", &slim_key);
    slim->Branch("h4l_entry", &slim_entry);
    Long64_t slim_total = 0, offset = 0, processed_entries = 0;
    for (const auto& c : counters) processed_entries += c.entries;
    g_read_error = false;
    // Every processed file must be in the chain with its first-pass entries
    // (a file that fails to open would be dropped silently by TChain).
    if (chain.GetEntries() != processed_entries) throw InputError("the slim chain does not hold the processed entries");
    check_read("the slim chain");
    for (std::size_t f = 0; f < paths.size(); ++f) {
      for (Long64_t entry : slim_entries[f]) {
        const Long64_t global = offset + entry;
        if (chain.LoadTree(global) != entry || chain.GetTreeNumber() != static_cast<int>(f))
          throw InputError("slim chain entry " + std::to_string(global) + " is not entry " + std::to_string(entry) + " of " + paths[f]);
        if (chain.GetEntry(global) <= 0) throw InputError("slim copy cannot read entry " + std::to_string(entry) + " of " + paths[f]);
        check_read(paths[f] + " (slim copy) entry " + std::to_string(entry));
        slim_key = keys[f];
        slim_entry = entry;
        if (slim->Fill() <= 0) throw std::runtime_error("slim Fill failed");
        ++slim_total;
      }
      offset += counters[f].entries;
    }

    // Files tree.
    std::string file_key_text, file_path_text;
    Long64_t f_entries = 0, f_pairs = 0, f_photon = 0, f_trigobj = 0, f_slim = 0;
    Double_t f_weight = 0;
    files_tree.Branch("file_key", &file_key_text);
    files_tree.Branch("path", &file_path_text);
    files_tree.Branch("entries", &f_entries);
    files_tree.Branch("pairs", &f_pairs);
    files_tree.Branch("photon_pairs", &f_photon);
    files_tree.Branch("trigobj_pairs", &f_trigobj);
    files_tree.Branch("slim", &f_slim);
    files_tree.Branch("sum_weight", &f_weight);
    json files_json = json::array();
    long long total_pairs = 0, total_photon = 0, total_trigobj = 0, total_entries = 0;
    for (const auto& c : counters) {
      file_key_text = c.file_key;
      file_path_text = c.path;
      f_entries = c.entries;
      f_pairs = c.pairs;
      f_photon = c.photon_pairs;
      f_trigobj = c.trigobj_pairs;
      f_slim = c.slim;
      f_weight = c.sum_weight;
      if (files_tree.Fill() <= 0) throw std::runtime_error("Files Fill failed");
      files_json.push_back({{"file_key", c.file_key}, {"path", c.path}, {"entries", c.entries}, {"pairs", c.pairs},
                            {"photon_pairs", c.photon_pairs}, {"trigobj_pairs", c.trigobj_pairs}, {"slim", c.slim},
                            {"sum_weight", c.sum_weight}});
      total_pairs += c.pairs;
      total_photon += c.photon_pairs;
      total_trigobj += c.trigobj_pairs;
      total_entries += c.entries;
    }

    output->cd();
    for (TTree* tree : {&pairs, &photon_pairs, &trigobj_pairs, &files_tree, slim})
      if (tree->Write() <= 0) throw std::runtime_error(std::string("cannot write ") + tree->GetName());
    if (config.is_signal && gen_table.Write() <= 0) throw std::runtime_error("cannot write GenTable");
    const long long gen_rows = config.is_signal ? gen_table.GetEntries() : 0;
    output->Close();
    output.reset();

    std::map<std::string, long long> trees = {{"Pairs", total_pairs},        {"PhotonPairs", total_photon},
                                              {"TrigObjPairs", total_trigobj}, {"Events", slim_total},
                                              {"Files", static_cast<long long>(counters.size())}};
    if (config.is_signal) trees["GenTable"] = total_entries;
    if (config.is_signal && gen_rows != total_entries) throw std::runtime_error("GenTable rows differ from the entries");
    h4l::validate_root_output_full(out_root, trees);
    // The slim entries carry exactly the planned (file key, entry) sequence.
    {
      std::unique_ptr<TFile> check(TFile::Open(out_root.c_str(), "READ"));
      auto* events = dynamic_cast<TTree*>(check->Get("Events"));
      events->SetBranchStatus("*", false);
      events->SetBranchStatus("h4l_file_key", true);
      events->SetBranchStatus("h4l_entry", true);
      ULong64_t read_key = 0;
      Long64_t read_entry = 0;
      if (events->SetBranchAddress("h4l_file_key", &read_key) < 0 || events->SetBranchAddress("h4l_entry", &read_entry) < 0)
        throw std::runtime_error("cannot address h4l_file_key/h4l_entry");
      Long64_t position = 0;
      for (std::size_t f = 0; f < paths.size(); ++f)
        for (Long64_t entry : slim_entries[f]) {
          if (events->GetEntry(position) <= 0 || read_key != keys[f] || read_entry != entry)
            throw std::runtime_error("slim entry " + std::to_string(position) + " does not carry its planned origin");
          ++position;
        }
    }

    json report = {{"schema", "h4l_v3_skim_report/2"},
                   {"task_id", task.at("task_id")},
                   {"sample", task.at("config").at("sample")},
                   {"kind", task.at("config").at("kind")},
                   {"analysis_config_version", task.at("analysis_config").at("version")},
                   {"analysis_config_fnv1a64", hex64(h4l::fnv1a64(task.at("analysis_config").dump()))},
                   {"frozen_program", task.value("frozen_program", json())},
                   {"flag_bits", flag_doc()},
                   {"trigger_bits", config.paths},
                   {"slim_branches_kept", kept.size()},
                   {"stage0_unknown", stage0_unknown},
                   {"skipped_files", skipped},
                   {"files", files_json},
                   {"totals", {{"entries", total_entries}, {"pairs", total_pairs}, {"photon_pairs", total_photon},
                               {"trigobj_pairs", total_trigobj}, {"slim", slim_total}}},
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
