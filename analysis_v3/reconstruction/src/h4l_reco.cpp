// h4l_reco: 4-lepton event records of the H -> ZZ* -> 4l analysis
// (AN-16-442 sections 3-4) from the stage-2 multilepton slim NanoAOD.
//
// Usage: h4l_reco --task TASK.json --out-json OUT.json --out-root OUT.root
//
// For the skim outputs of one task (all of one sample) it reads the slim
// Events tree and writes one row of Events4l for every event with at least
// three loose leptons and an opposite-sign same-flavour pair with
// 30 < m(ll gamma) < 120 GeV:
//   * the event: original file key and entry, genWeight, the skim trigger
//     bits and the analysis trigger OR, MET, good vertices;
//   * every loose lepton, ordered by decreasing calibrated pT (at most 12):
//       muons: pT > floor, |eta| < 2.4, |dxy| < 0.5, |dz| < 1, global or
//       tracker with a matched station; ghost cleaning (substitute for the
//       AN segment sharing, not in NanoAOD): of two muons within dR < 0.02
//       the PF one, else the higher-pT one, is kept;
//       electrons: pT > floor, |eta| < 2.5, |dxy| < 0.5, |dz| < 1; the AN
//       cross cleaning against selected muons is applied downstream with
//       l_overlap (the closest loose muon within dR < 0.05);
//     with the calibrated pT (data pT exp(-u), MC pT (1 + r N), N the
//     stream-0 deviate of (file key, entry, collection, original index);
//     h4l/calibration.h; the eta variable of each payload is checked), the
//     raw pT, N, the relative pT error, the isolation (all and charged) and
//     the isolation with the AN FSR photons of all loose muons inside
//     0.01 < dR < 0.3 removed from its neutral part (the downstream
//     selection can rebuild other definitions from the stored parts), SIP,
//     dxy, dz, identification flags (muons in the stage-2 layout of
//     h4l/calibration.h) and scores, the AN FSR photon (muons only: pT > 2,
//     |eta| < 2.4, relIso < 1.8, dR/ET^2 < 0.012, dR < 0.5, the lowest
//     dR/ET^2), the matched HLT objects (bits and pT, the skim matching), the
//     generator match, and for electrons the loose muons within dR < 0.05
//     (closest index and bitmask);
//   * every jet with pT above the floor (20 GeV, so that jet energy
//     variations can cross the analysis threshold of 30 GeV), |eta| < 4.7 and
//     the tight jet ID (UL NanoAOD has no loose jet ID), with CSVv2 and
//     DeepJet scores;
//   * every ZZ candidate under the loosest thresholds the analysis scans:
//     Z1 an opposite-sign same-flavour pair with 30 < m(ll gamma) < 130 GeV,
//     Z2 another same-flavour pair (opposite or same sign; same sign for the
//     SS control region) with 3 < m(ll gamma) < 130 GeV, four distinct
//     leptons with SIP < 8, dR > 0.02 between all, the leading and
//     subleading pT above 15 and 7 GeV, every opposite-sign pair above 2 GeV
//     (without FSR) and m4l > 65 GeV (floors below the analysis thresholds
//     so that scale variations stay exact).  Both orderings of two opposite-sign
//     pairs are stored; c_z1_closer marks the one with Z1 closer to m_Z (the
//     signal-region definition; an exact tie goes to the pair with the lower
//     first index).  For 4e and 4mu with an opposite-sign Z2,
//     the alternative pairing (Za closer to m_Z, Zb) is stored for the smart
//     cut;
//   * MC: the LHE scale and PDF weights and the parton-shower weights of
//     events with a candidate; signal: the generator-table fields.
// The final selection (the AN or the optimized thresholds, the tight
// lepton, the cross cleaning, the smart cut, the candidate choice by
// D_bkg^kin), the jet cleaning, the categories, D_mass, the Z1 refit, the
// weights and the variations are applied downstream.  Any read error or
// disagreement with the planned slim-event counts aborts the task; the
// output is fully re-read before the JSON report is written.

#include "h4l/calibration.h"
#include "h4l/hash.h"
#include "h4l/io.h"
#include "h4l/kinematics.h"
#include "h4l/root_io.h"
#include "h4l/trigger.h"

#include <TError.h>
#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {
using h4l::json;
using h4l::P4;

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

template <typename T>
using Array = TTreeReaderArray<T>;
template <typename T>
using Value = TTreeReaderValue<T>;

constexpr int kMaxLeptons = 12, kMaxJets = 24, kMaxCandidates = 1024, kMaxScale = 16, kMaxPdf = 128, kMaxPs = 8;

struct Reader {
  TTreeReader reader;
  Value<ULong64_t> key;
  Value<Long64_t> entry;
  Value<Float_t> met_pt, met_phi;
  Value<Int_t> npv_good;
  Array<Float_t> mu_pt, mu_eta, mu_phi, mu_dxy, mu_dz, mu_sip, mu_iso, mu_iso_chg, mu_pterr;
  Array<Int_t> mu_charge, mu_nstations;
  Array<Bool_t> mu_global, mu_tracker, mu_standalone, mu_pf, mu_loose, mu_medium, mu_tight, mu_soft, mu_trigloose, mu_intime,
      mu_medium_prompt, mu_high_purity, mu_soft_mva_id;
  Array<UChar_t> mu_highpt;
  Value<Float_t> rho;
  Array<Float_t> el_pt, el_eta, el_phi, el_deta_sc, el_dxy, el_dz, el_sip, el_iso, el_iso_chg, el_energy_err, el_mva_noiso,
      el_mva_iso;
  Array<Int_t> el_charge, el_cutbased;
  Array<Bool_t> el_wp90_noiso, el_wp80_noiso, el_wpl_noiso, el_wp90_iso, el_wp80_iso, el_wpl_iso, el_convveto;
  Array<UChar_t> el_losthits;
  Array<Float_t> fsr_pt, fsr_eta, fsr_phi, fsr_iso, fsr_dret2;
  Array<Int_t> fsr_muon;
  Array<Float_t> jet_pt, jet_eta, jet_phi, jet_mass, jet_csv, jet_deepjet;
  Array<Int_t> jet_id, jet_puid;
  Array<Float_t> to_pt, to_eta, to_phi;
  Array<Int_t> to_id, to_bits;
  std::vector<std::unique_ptr<Value<Bool_t>>> paths;
  // MC only.
  std::unique_ptr<Value<Float_t>> gen_weight, pu_true;
  std::unique_ptr<Array<UChar_t>> mu_gen_flav, el_gen_flav;
  std::unique_ptr<Array<Int_t>> mu_gen_idx, el_gen_idx, jet_hadron_flavour, jet_gen_idx;
  std::unique_ptr<Array<Float_t>> gen_pt, lhe_scale, lhe_pdf, ps;

  Reader(TTree* tree, bool mc, const std::vector<std::string>& trigger_paths)
      : reader(tree),
        key(reader, "h4l_file_key"), entry(reader, "h4l_entry"), met_pt(reader, "MET_pt"), met_phi(reader, "MET_phi"),
        npv_good(reader, "PV_npvsGood"),
        mu_pt(reader, "Muon_pt"), mu_eta(reader, "Muon_eta"), mu_phi(reader, "Muon_phi"), mu_dxy(reader, "Muon_dxy"),
        mu_dz(reader, "Muon_dz"), mu_sip(reader, "Muon_sip3d"), mu_iso(reader, "Muon_pfRelIso03_all"),
        mu_iso_chg(reader, "Muon_pfRelIso03_chg"), mu_pterr(reader, "Muon_ptErr"), mu_charge(reader, "Muon_charge"),
        mu_nstations(reader, "Muon_nStations"), mu_global(reader, "Muon_isGlobal"), mu_tracker(reader, "Muon_isTracker"),
        mu_standalone(reader, "Muon_isStandalone"), mu_pf(reader, "Muon_isPFcand"), mu_loose(reader, "Muon_looseId"),
        mu_medium(reader, "Muon_mediumId"), mu_tight(reader, "Muon_tightId"), mu_soft(reader, "Muon_softId"),
        mu_trigloose(reader, "Muon_triggerIdLoose"), mu_intime(reader, "Muon_inTimeMuon"),
        mu_medium_prompt(reader, "Muon_mediumPromptId"), mu_high_purity(reader, "Muon_highPurity"),
        mu_soft_mva_id(reader, "Muon_softMvaId"), mu_highpt(reader, "Muon_highPtId"), rho(reader, "fixedGridRhoFastjetAll"),
        el_pt(reader, "Electron_pt"), el_eta(reader, "Electron_eta"), el_phi(reader, "Electron_phi"),
        el_deta_sc(reader, "Electron_deltaEtaSC"), el_dxy(reader, "Electron_dxy"), el_dz(reader, "Electron_dz"),
        el_sip(reader, "Electron_sip3d"), el_iso(reader, "Electron_pfRelIso03_all"),
        el_iso_chg(reader, "Electron_pfRelIso03_chg"), el_energy_err(reader, "Electron_energyErr"),
        el_mva_noiso(reader, "Electron_mvaFall17V2noIso"), el_mva_iso(reader, "Electron_mvaFall17V2Iso"),
        el_charge(reader, "Electron_charge"), el_cutbased(reader, "Electron_cutBased"),
        el_wp90_noiso(reader, "Electron_mvaFall17V2noIso_WP90"), el_wp80_noiso(reader, "Electron_mvaFall17V2noIso_WP80"),
        el_wpl_noiso(reader, "Electron_mvaFall17V2noIso_WPL"), el_wp90_iso(reader, "Electron_mvaFall17V2Iso_WP90"),
        el_wp80_iso(reader, "Electron_mvaFall17V2Iso_WP80"), el_wpl_iso(reader, "Electron_mvaFall17V2Iso_WPL"),
        el_convveto(reader, "Electron_convVeto"), el_losthits(reader, "Electron_lostHits"),
        fsr_pt(reader, "FsrPhoton_pt"), fsr_eta(reader, "FsrPhoton_eta"), fsr_phi(reader, "FsrPhoton_phi"),
        fsr_iso(reader, "FsrPhoton_relIso03"), fsr_dret2(reader, "FsrPhoton_dROverEt2"), fsr_muon(reader, "FsrPhoton_muonIdx"),
        jet_pt(reader, "Jet_pt"), jet_eta(reader, "Jet_eta"), jet_phi(reader, "Jet_phi"), jet_mass(reader, "Jet_mass"),
        jet_csv(reader, "Jet_btagCSVV2"), jet_deepjet(reader, "Jet_btagDeepFlavB"), jet_id(reader, "Jet_jetId"),
        jet_puid(reader, "Jet_puId"), to_pt(reader, "TrigObj_pt"), to_eta(reader, "TrigObj_eta"), to_phi(reader, "TrigObj_phi"),
        to_id(reader, "TrigObj_id"), to_bits(reader, "TrigObj_filterBits") {
    for (const auto& path : trigger_paths) paths.push_back(std::make_unique<Value<Bool_t>>(reader, path.c_str()));
    if (mc) {
      gen_weight = std::make_unique<Value<Float_t>>(reader, "genWeight");
      pu_true = std::make_unique<Value<Float_t>>(reader, "Pileup_nTrueInt");
      mu_gen_flav = std::make_unique<Array<UChar_t>>(reader, "Muon_genPartFlav");
      el_gen_flav = std::make_unique<Array<UChar_t>>(reader, "Electron_genPartFlav");
      mu_gen_idx = std::make_unique<Array<Int_t>>(reader, "Muon_genPartIdx");
      el_gen_idx = std::make_unique<Array<Int_t>>(reader, "Electron_genPartIdx");
      gen_pt = std::make_unique<Array<Float_t>>(reader, "GenPart_pt");
      jet_hadron_flavour = std::make_unique<Array<Int_t>>(reader, "Jet_hadronFlavour");
      jet_gen_idx = std::make_unique<Array<Int_t>>(reader, "Jet_genJetIdx");
      // Theory weights exist only in some samples (none in the private signal samples).
      if (tree->GetBranch("LHEScaleWeight")) lhe_scale = std::make_unique<Array<Float_t>>(reader, "LHEScaleWeight");
      if (tree->GetBranch("LHEPdfWeight")) lhe_pdf = std::make_unique<Array<Float_t>>(reader, "LHEPdfWeight");
      if (tree->GetBranch("PSWeight")) ps = std::make_unique<Array<Float_t>>(reader, "PSWeight");
    }
  }
};

// One loose lepton.
struct Lepton {
  int pdg = 0, charge = 0, index = 0, overlap = -1, cut_based = 0, lost_hits = 0, gen_flav = -1, trig_bits = 0;
  double pt = 0, pt_raw = 0, eta = 0, eta_sc = 0, phi = 0, normal = 0, rel_err = 0, iso = 0, iso_chg = 0, iso_fsr = 0, sip = 0,
         dxy = 0, dz = 0, mva_noiso = -9, mva_iso = -9, gen_pt = -1, trig_pt = -1;
  unsigned flags = 0;  // muons: the stage-2 layout (0 isGlobal, 1 isTracker, 2 isStandalone, 3 isPFcand, 4 looseId,
                       // 5 mediumId, 6 tightId, 7 softId, 8 highPtId >= 1, 9 highPtId == 2, 10 nStations > 0,
                       // 11 triggerIdLoose, 12 inTimeMuon, 13 mediumPromptId, 14 highPurity, 15 softMvaId);
                       // electrons: 0/1/2 noIso WPL/90/80, 3/4/5 Iso WPL/90/80, 6 convVeto
  unsigned overlap_mask = 0;  // electrons: bit k = loose muon k within the overlap cone
  double fsr_pt = -1, fsr_eta = 0, fsr_phi = 0, fsr_dret2 = -1;
  bool muon() const { return std::abs(pdg) == 13; }
  P4 p4() const { return P4(pt, eta, phi, muon() ? h4l::kMuonMass : h4l::kElectronMass); }
  P4 p4_fsr() const {
    P4 p = p4();
    if (fsr_pt > 0) p += P4(fsr_pt, fsr_eta, fsr_phi, 0.0);
    return p;
  }
};

struct Candidate {
  std::array<int, 4> leptons{};  // Z1 legs, Z2 legs
  double m_z1 = 0, m_z2 = 0, m4l = 0, m_za = -1, m_zb = -1;
  bool z2_same_sign = false, z1_closer = true;
};

struct Config {
  double mu_floor_pt = 3, el_floor_pt = 5, mu_abs_eta = 2.4, el_abs_eta = 2.5, abs_dxy = 0.5, abs_dz = 1.0, ghost_dr = 0.02,
         overlap_dr = 0.05;
  double fsr_pt = 2, fsr_eta = 2.4, fsr_iso = 1.8, fsr_dret2 = 0.012, fsr_dr = 0.5, iso_veto_dr = 0.01, iso_dr = 0.3;
  int min_leptons = 3;
  double event_z_low = 30, event_z_high = 120;
  double cand_sip = 8, z1_low = 30, z1_high = 120, z2_low = 4, z2_high = 120, lead_pt = 15, sublead_pt = 7, min_dr = 0.02,
         min_os_mass = 2, min_m4l = 70;
  double jet_pt = 30, jet_eta = 4.7;
  int jet_id_mask = 2;
  double match_dr = 0.1, muon_merge_dr = 0.05;
  h4l::Bending bending;
};

Config read_config(const json& reco, const json& triggers) {
  Config c;
  const json& l = reco.at("leptons");
  c.mu_floor_pt = l.at("muon_floor_pt").get<double>();
  c.el_floor_pt = l.at("electron_floor_pt").get<double>();
  c.mu_abs_eta = l.at("muon_abs_eta").get<double>();
  c.el_abs_eta = l.at("electron_abs_eta").get<double>();
  c.abs_dxy = l.at("abs_dxy").get<double>();
  c.abs_dz = l.at("abs_dz").get<double>();
  c.ghost_dr = l.at("ghost_dr").get<double>();
  c.overlap_dr = l.at("overlap_dr").get<double>();
  const json& f = reco.at("fsr");
  c.fsr_pt = f.at("pt").get<double>();
  c.fsr_eta = f.at("abs_eta").get<double>();
  c.fsr_iso = f.at("rel_iso").get<double>();
  c.fsr_dret2 = f.at("dr_over_et2").get<double>();
  c.fsr_dr = f.at("max_dr").get<double>();
  c.iso_veto_dr = f.at("iso_veto_dr").get<double>();
  c.iso_dr = f.at("iso_dr").get<double>();
  const json& e = reco.at("event");
  c.min_leptons = e.at("min_leptons").get<int>();
  c.event_z_low = e.at("z_mass").at(0).get<double>();
  c.event_z_high = e.at("z_mass").at(1).get<double>();
  const json& k = reco.at("candidates");
  c.cand_sip = k.at("max_sip").get<double>();
  c.z1_low = k.at("z1_mass").at(0).get<double>();
  c.z1_high = k.at("z1_mass").at(1).get<double>();
  c.z2_low = k.at("z2_mass").at(0).get<double>();
  c.z2_high = k.at("z2_mass").at(1).get<double>();
  c.lead_pt = k.at("leading_pt").get<double>();
  c.sublead_pt = k.at("subleading_pt").get<double>();
  c.min_dr = k.at("min_dr").get<double>();
  c.min_os_mass = k.at("min_os_mass").get<double>();
  c.min_m4l = k.at("min_m4l").get<double>();
  const json& j = reco.at("jets");
  c.jet_pt = j.at("pt").get<double>();
  c.jet_eta = j.at("abs_eta").get<double>();
  c.jet_id_mask = j.at("jet_id_mask").get<int>();
  c.match_dr = triggers.at("trigobj_match_dr").get<double>();
  c.muon_merge_dr = triggers.at("muon_object_merge_dr").get<double>();
  const json& position = triggers.at("electron_object_position");
  c.bending = {position.at("field_tesla").get<double>(), position.at("barrel_radius_m").get<double>(),
               position.at("endcap_z_m").get<double>()};
  return c;
}

struct GenFields {
  Int_t stage0 = 0, vh_class = 0, final_state_true = -1, fid_final_state = -1, fid_njets = -1;
  Bool_t fid_pass = false;
  Float_t h_pt = -1, h_eta = 0, h_phi = 0, h_mass = -1, h_y = 0, m4l_dressed = -1, fid_m4l = -1, fid_mz1 = -1, fid_mz2 = -1,
          fid_pt4l = -1, fid_y4l = 0, fid_jet1_pt = -1;
  std::array<Float_t, 4> hlep_pt{}, hlep_eta{}, hlep_phi{};
  std::array<Int_t, 4> hlep_pdg{};
};

// Output buffers of one Events4l row.
struct Row {
  ULong64_t file_key = 0;
  Long64_t entry = 0;
  Float_t weight = 1, met = 0, met_phi = 0, rho = 0, pu_true = -1;
  Int_t npv = 0, nlep_all = 0, njet_all = 0;
  UInt_t trig_bits = 0;
  Bool_t trigger = false;
  Int_t nlep = 0;
  Int_t l_pdg[kMaxLeptons], l_index[kMaxLeptons], l_overlap[kMaxLeptons], l_cut_based[kMaxLeptons], l_lost_hits[kMaxLeptons],
      l_gen_flav[kMaxLeptons], l_trig_bits[kMaxLeptons];
  UShort_t l_flags[kMaxLeptons], l_overlap_mask[kMaxLeptons];
  Float_t l_pt[kMaxLeptons], l_pt_raw[kMaxLeptons], l_eta[kMaxLeptons], l_eta_sc[kMaxLeptons], l_phi[kMaxLeptons],
      l_normal[kMaxLeptons], l_rel_err[kMaxLeptons], l_iso[kMaxLeptons], l_iso_chg[kMaxLeptons], l_iso_fsr[kMaxLeptons],
      l_sip[kMaxLeptons],
      l_dxy[kMaxLeptons], l_dz[kMaxLeptons], l_mva_noiso[kMaxLeptons], l_mva_iso[kMaxLeptons], l_gen_pt[kMaxLeptons],
      l_trig_pt[kMaxLeptons], l_fsr_pt[kMaxLeptons], l_fsr_eta[kMaxLeptons], l_fsr_phi[kMaxLeptons], l_fsr_dret2[kMaxLeptons];
  Int_t njet = 0;
  Float_t j_pt[kMaxJets], j_eta[kMaxJets], j_phi[kMaxJets], j_mass[kMaxJets], j_csv[kMaxJets], j_deepjet[kMaxJets];
  Int_t j_id[kMaxJets], j_puid[kMaxJets], j_flavour[kMaxJets], j_gen[kMaxJets];
  Int_t ncand = 0, ncand_overflow = 0;
  Char_t c_l1[kMaxCandidates], c_l2[kMaxCandidates], c_l3[kMaxCandidates], c_l4[kMaxCandidates];
  Bool_t c_ss[kMaxCandidates], c_z1_closer[kMaxCandidates];
  Float_t c_mz1[kMaxCandidates], c_mz2[kMaxCandidates], c_m4l[kMaxCandidates], c_mza[kMaxCandidates], c_mzb[kMaxCandidates];
  Int_t n_scale = 0, n_pdf = 0, n_ps = 0;
  Float_t w_scale[kMaxScale], w_pdf[kMaxPdf], w_ps[kMaxPs];
  Bool_t g_present = false;
  GenFields g;

  void book(TTree& t) {
    t.Branch("file_key", &file_key, "file_key/l");
    t.Branch("entry", &entry, "entry/L");
    t.Branch("weight", &weight, "weight/F");
    t.Branch("met", &met, "met/F");
    t.Branch("met_phi", &met_phi, "met_phi/F");
    t.Branch("npv", &npv, "npv/I");
    t.Branch("rho", &rho, "rho/F");
    t.Branch("pu_true", &pu_true, "pu_true/F");
    t.Branch("nlep_all", &nlep_all, "nlep_all/I");
    t.Branch("njet_all", &njet_all, "njet_all/I");
    t.Branch("trig_bits", &trig_bits, "trig_bits/i");
    t.Branch("trigger", &trigger, "trigger/O");
    t.Branch("nlep", &nlep, "nlep/I");
    for (auto [name, a] : std::initializer_list<std::pair<const char*, Int_t*>>{
             {"l_pdg", l_pdg}, {"l_index", l_index}, {"l_overlap", l_overlap}, {"l_cut_based", l_cut_based},
             {"l_lost_hits", l_lost_hits}, {"l_gen_flav", l_gen_flav}, {"l_trig_bits", l_trig_bits}})
      t.Branch(name, a, (std::string(name) + "[nlep]/I").c_str());
    t.Branch("l_flags", l_flags, "l_flags[nlep]/s");
    t.Branch("l_overlap_mask", l_overlap_mask, "l_overlap_mask[nlep]/s");
    for (auto [name, a] : std::initializer_list<std::pair<const char*, Float_t*>>{
             {"l_pt", l_pt}, {"l_pt_raw", l_pt_raw}, {"l_eta", l_eta}, {"l_eta_sc", l_eta_sc}, {"l_phi", l_phi},
             {"l_normal", l_normal}, {"l_rel_err", l_rel_err}, {"l_iso", l_iso}, {"l_iso_chg", l_iso_chg},
             {"l_iso_fsr", l_iso_fsr}, {"l_sip", l_sip},
             {"l_dxy", l_dxy}, {"l_dz", l_dz}, {"l_mva_noiso", l_mva_noiso}, {"l_mva_iso", l_mva_iso}, {"l_gen_pt", l_gen_pt},
             {"l_trig_pt", l_trig_pt}, {"l_fsr_pt", l_fsr_pt}, {"l_fsr_eta", l_fsr_eta}, {"l_fsr_phi", l_fsr_phi},
             {"l_fsr_dret2", l_fsr_dret2}})
      t.Branch(name, a, (std::string(name) + "[nlep]/F").c_str());
    t.Branch("njet", &njet, "njet/I");
    for (auto [name, a] : std::initializer_list<std::pair<const char*, Float_t*>>{
             {"j_pt", j_pt}, {"j_eta", j_eta}, {"j_phi", j_phi}, {"j_mass", j_mass}, {"j_csv", j_csv}, {"j_deepjet", j_deepjet}})
      t.Branch(name, a, (std::string(name) + "[njet]/F").c_str());
    for (auto [name, a] : std::initializer_list<std::pair<const char*, Int_t*>>{
             {"j_id", j_id}, {"j_puid", j_puid}, {"j_flavour", j_flavour}, {"j_gen", j_gen}})
      t.Branch(name, a, (std::string(name) + "[njet]/I").c_str());
    t.Branch("ncand", &ncand, "ncand/I");
    t.Branch("ncand_overflow", &ncand_overflow, "ncand_overflow/I");
    for (auto [name, a] : std::initializer_list<std::pair<const char*, Char_t*>>{
             {"c_l1", c_l1}, {"c_l2", c_l2}, {"c_l3", c_l3}, {"c_l4", c_l4}})
      t.Branch(name, a, (std::string(name) + "[ncand]/B").c_str());
    t.Branch("c_ss", c_ss, "c_ss[ncand]/O");
    t.Branch("c_z1_closer", c_z1_closer, "c_z1_closer[ncand]/O");
    for (auto [name, a] : std::initializer_list<std::pair<const char*, Float_t*>>{
             {"c_mz1", c_mz1}, {"c_mz2", c_mz2}, {"c_m4l", c_m4l}, {"c_mza", c_mza}, {"c_mzb", c_mzb}})
      t.Branch(name, a, (std::string(name) + "[ncand]/F").c_str());
    t.Branch("n_scale", &n_scale, "n_scale/I");
    t.Branch("w_scale", w_scale, "w_scale[n_scale]/F");
    t.Branch("n_pdf", &n_pdf, "n_pdf/I");
    t.Branch("w_pdf", w_pdf, "w_pdf[n_pdf]/F");
    t.Branch("n_ps", &n_ps, "n_ps/I");
    t.Branch("w_ps", w_ps, "w_ps[n_ps]/F");
    t.Branch("g_present", &g_present, "g_present/O");
    for (auto [name, a] : std::initializer_list<std::pair<const char*, Int_t*>>{
             {"g_stage0", &g.stage0}, {"g_vh_class", &g.vh_class}, {"g_final_state_true", &g.final_state_true},
             {"g_fid_final_state", &g.fid_final_state}, {"g_fid_njets", &g.fid_njets}})
      t.Branch(name, a, (std::string(name) + "/I").c_str());
    t.Branch("g_fid_pass", &g.fid_pass, "g_fid_pass/O");
    for (auto [name, a] : std::initializer_list<std::pair<const char*, Float_t*>>{
             {"g_h_pt", &g.h_pt}, {"g_h_eta", &g.h_eta}, {"g_h_phi", &g.h_phi}, {"g_h_mass", &g.h_mass}, {"g_h_y", &g.h_y},
             {"g_m4l_dressed", &g.m4l_dressed}, {"g_fid_m4l", &g.fid_m4l}, {"g_fid_mz1", &g.fid_mz1},
             {"g_fid_mz2", &g.fid_mz2}, {"g_fid_pt4l", &g.fid_pt4l}, {"g_fid_y4l", &g.fid_y4l},
             {"g_fid_jet1_pt", &g.fid_jet1_pt}})
      t.Branch(name, a, (std::string(name) + "/F").c_str());
    t.Branch("g_hlep_pt", g.hlep_pt.data(), "g_hlep_pt[4]/F");
    t.Branch("g_hlep_eta", g.hlep_eta.data(), "g_hlep_eta[4]/F");
    t.Branch("g_hlep_phi", g.hlep_phi.data(), "g_hlep_phi[4]/F");
    t.Branch("g_hlep_pdg", g.hlep_pdg.data(), "g_hlep_pdg[4]/I");
  }
};

// The generator table of a signal skim, keyed by (file key, entry).
std::map<std::pair<ULong64_t, Long64_t>, GenFields> read_gen_table(TFile& file, const std::string& path, long long expected) {
  auto* tree = dynamic_cast<TTree*>(file.Get("GenTable"));
  if (!tree) throw InputError("no GenTable in " + path);
  if (tree->GetEntries() != expected)
    throw InputError(path + ": GenTable has " + std::to_string(tree->GetEntries()) + " rows, the skim report " +
                     std::to_string(expected));
  ULong64_t key = 0;
  Long64_t entry = 0;
  GenFields g;
  auto set = [&](const char* name, void* where) {
    if (tree->SetBranchAddress(name, where) < 0) throw InputError(std::string("cannot address GenTable ") + name);
  };
  set("file_key", &key);
  set("entry", &entry);
  for (auto [name, a] : std::initializer_list<std::pair<const char*, void*>>{
           {"stage0", &g.stage0}, {"vh_class", &g.vh_class}, {"final_state_true", &g.final_state_true},
           {"fid_final_state", &g.fid_final_state}, {"fid_njets", &g.fid_njets}, {"fid_pass", &g.fid_pass},
           {"h_pt", &g.h_pt}, {"h_eta", &g.h_eta}, {"h_phi", &g.h_phi}, {"h_mass", &g.h_mass}, {"h_y", &g.h_y},
           {"m4l_dressed", &g.m4l_dressed}, {"fid_m4l", &g.fid_m4l}, {"fid_mz1", &g.fid_mz1}, {"fid_mz2", &g.fid_mz2},
           {"fid_pt4l", &g.fid_pt4l}, {"fid_y4l", &g.fid_y4l}, {"fid_jet1_pt", &g.fid_jet1_pt},
           {"hlep_pt", g.hlep_pt.data()}, {"hlep_eta", g.hlep_eta.data()}, {"hlep_phi", g.hlep_phi.data()},
           {"hlep_pdg", g.hlep_pdg.data()}})
    set(name, a);
  std::map<std::pair<ULong64_t, Long64_t>, GenFields> table;
  for (Long64_t r = 0; r < tree->GetEntries(); ++r) {
    if (tree->GetEntry(r) <= 0) throw InputError("cannot read GenTable row of " + path);
    if (!table.emplace(std::make_pair(key, entry), g).second) throw InputError("duplicate GenTable row in " + path);
  }
  check_read(path + " GenTable");
  tree->ResetBranchAddresses();
  return table;
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
    std::cerr << "usage: h4l_reco --task TASK.json --out-json OUT.json --out-root OUT.root\n";
    return 64;
  }
  gROOT->SetBatch(true);
  SetErrorHandler(record_errors);
  try {
    const json task = h4l::read_json(task_path);
    const json& triggers = task.at("triggers");
    const Config cfg = read_config(task.at("reco_config"), triggers);
    const bool is_mc = task.at("config").at("kind").get<std::string>() == "mc";
    const bool is_signal = task.at("config").at("role").get<std::string>() == "signal";
    const auto muon_payload = h4l::FactorizedPayload::from_json(task.at("calibration_payload").at("muon"));
    const auto electron_payload = h4l::FactorizedPayload::from_json(task.at("calibration_payload").at("electron"));
    if (muon_payload.model.eta_sc || !electron_payload.model.eta_sc)
      throw std::runtime_error("the payloads must use |eta| for muons and |eta_SC| for electrons");
    // The skim trigger bits are the order of triggers.paths; the analysis OR is triggers.analysis_or.
    std::vector<std::string> trigger_paths;
    for (const auto& name : triggers.at("paths")) trigger_paths.push_back(name.get<std::string>());
    if (trigger_paths.size() > 32) throw std::runtime_error("more than 32 trigger paths");
    UInt_t analysis_mask = 0;
    for (const auto& name : triggers.at("analysis_or")) {
      const auto found = std::find(trigger_paths.begin(), trigger_paths.end(), name.get<std::string>());
      if (found == trigger_paths.end()) throw std::runtime_error("analysis path not in the skim paths: " + name.get<std::string>());
      analysis_mask |= 1U << (found - trigger_paths.begin());
    }

    std::set<ULong64_t> planned_keys;
    for (const auto& key : task.at("source_file_keys")) planned_keys.insert(h4l::parse_file_key(key.get<std::string>()));

    std::unique_ptr<TFile> output(TFile::Open(out_root.c_str(), "RECREATE"));
    if (!output || output->IsZombie()) throw std::runtime_error("cannot create " + out_root);
    output->cd();
    TTree& events_out = *new TTree("Events4l", "4-lepton event records");
    TTree& inputs_tree = *new TTree("Inputs", "per-input counters");
    Row row;
    row.book(events_out);
    std::string input_path;
    Long64_t input_events = 0, input_kept = 0, input_with_candidate = 0, input_candidates = 0;
    inputs_tree.Branch("path", &input_path);
    inputs_tree.Branch("events", &input_events, "events/L");
    inputs_tree.Branch("kept", &input_kept, "kept/L");
    inputs_tree.Branch("with_candidate", &input_with_candidate, "with_candidate/L");
    inputs_tree.Branch("candidates", &input_candidates, "candidates/L");

    json inputs_report = json::array();
    long long total_events = 0, total_kept = 0, total_with_candidate = 0, total_candidates = 0, candidate_overflow = 0,
              lepton_overflow = 0, jet_overflow = 0;
    for (const auto& input : task.at("inputs")) {
      input_path = input.at("root").get<std::string>();
      const long long expected = input.at("events").get<long long>();
      std::string open_error;
      auto file = h4l::open_input(input_path, 5, open_error);
      if (!file) throw InputError("cannot open " + input_path + ": " + open_error);
      g_read_error = false;
      auto* tree = dynamic_cast<TTree*>(file->Get("Events"));
      if (!tree) throw InputError("no Events in " + input_path);
      if (tree->GetEntries() != expected)
        throw InputError(input_path + " has " + std::to_string(tree->GetEntries()) + " slim events, the skim report " +
                         std::to_string(expected));
      std::map<std::pair<ULong64_t, Long64_t>, GenFields> gen;
      if (is_signal) gen = read_gen_table(*file, input_path, input.at("gen_rows").get<long long>());
      // Theory-weight branches: required for the roles that use them, recorded for all.
      json theory = json::object();
      for (const char* name : {"LHEScaleWeight", "LHEPdfWeight", "PSWeight"}) theory[name] = {{"present", tree->GetBranch(name) != nullptr}};
      if (is_mc) {
        for (const auto& name : task.at("required_theory_branches"))
          if (!tree->GetBranch(name.get<std::string>().c_str()))
            throw InputError("required theory branch " + name.get<std::string>() + " missing in " + input_path);
      }
      int length_min[3] = {1 << 30, 1 << 30, 1 << 30}, length_max[3] = {-1, -1, -1};
      input_events = tree->GetEntries();
      input_kept = input_with_candidate = input_candidates = 0;
      {
      Reader in(tree, is_mc, trigger_paths);
      for (Long64_t e = 0; e < input_events; ++e) {
        if (in.reader.SetEntry(e) != TTreeReader::kEntryValid)
          throw InputError("read error at entry " + std::to_string(e) + " of " + input_path);
        const ULong64_t key = *in.key;
        const Long64_t entry = *in.entry;
        if (!planned_keys.count(key)) throw InputError("event of a file key outside the plan in " + input_path);
        const auto muon_objects =
            h4l::trigger_objects(in.to_pt, in.to_eta, in.to_phi, in.to_id, in.to_bits, 13, cfg.muon_merge_dr);
        const auto electron_objects = h4l::trigger_objects(in.to_pt, in.to_eta, in.to_phi, in.to_id, in.to_bits, 11, 0.0);
        // Loose muons (calibrated), ghost-cleaned.
        std::vector<Lepton> muons;
        for (std::size_t i = 0; i < in.mu_pt.GetSize(); ++i) {
          const double abs_eta = std::fabs(in.mu_eta[i]);
          if (abs_eta >= cfg.mu_abs_eta || std::fabs(in.mu_dxy[i]) >= cfg.abs_dxy || std::fabs(in.mu_dz[i]) >= cfg.abs_dz) continue;
          if (!(in.mu_global[i] || (in.mu_tracker[i] && in.mu_nstations[i] > 0))) continue;
          Lepton l;
          l.pdg = -13 * in.mu_charge[i];
          l.charge = in.mu_charge[i];
          l.index = static_cast<int>(i);
          l.pt_raw = in.mu_pt[i];
          l.normal = is_mc ? h4l::object_normal(key, static_cast<std::uint64_t>(entry), h4l::Collection::Muon, i, 0) : 0.0;
          l.pt = is_mc ? muon_payload.smeared_mc_pt(l.pt_raw, abs_eta, l.normal) : muon_payload.corrected_data_pt(l.pt_raw, abs_eta);
          if (l.pt <= cfg.mu_floor_pt) continue;
          l.eta = in.mu_eta[i];
          l.eta_sc = in.mu_eta[i];
          l.phi = in.mu_phi[i];
          l.rel_err = in.mu_pterr[i] / in.mu_pt[i];
          l.iso = in.mu_iso[i];
          l.iso_chg = in.mu_iso_chg[i];
          l.sip = in.mu_sip[i];
          l.dxy = in.mu_dxy[i];
          l.dz = in.mu_dz[i];
          const bool bits[] = {in.mu_global[i],        in.mu_tracker[i],     in.mu_standalone[i],    in.mu_pf[i],
                               in.mu_loose[i],         in.mu_medium[i],      in.mu_tight[i],         in.mu_soft[i],
                               in.mu_highpt[i] >= 1,   in.mu_highpt[i] == 2, in.mu_nstations[i] > 0, in.mu_trigloose[i],
                               in.mu_intime[i],        in.mu_medium_prompt[i], in.mu_high_purity[i], in.mu_soft_mva_id[i]};
          for (std::size_t b = 0; b < std::size(bits); ++b) l.flags |= static_cast<unsigned>(bits[b]) << b;
          const auto match = h4l::match_trigger(muon_objects, l.eta, l.phi, cfg.match_dr);
          l.trig_bits = match.bits;
          l.trig_pt = match.pt;
          if (is_mc) {
            l.gen_flav = (*in.mu_gen_flav)[i];
            const int g = (*in.mu_gen_idx)[i];
            if (g >= 0 && static_cast<std::size_t>(g) < in.gen_pt->GetSize()) l.gen_pt = (*in.gen_pt)[g];
          }
          // AN FSR photon: the lowest dR/ET^2 among the photons of this muon.
          for (std::size_t k = 0; k < in.fsr_pt.GetSize(); ++k) {
            if (in.fsr_muon[k] != static_cast<int>(i)) continue;
            const double dr = h4l::delta_r(l.eta, l.phi, in.fsr_eta[k], in.fsr_phi[k]);
            if (in.fsr_pt[k] <= cfg.fsr_pt || std::fabs(in.fsr_eta[k]) >= cfg.fsr_eta || in.fsr_iso[k] >= cfg.fsr_iso ||
                in.fsr_dret2[k] >= cfg.fsr_dret2 || dr >= cfg.fsr_dr)
              continue;
            if (l.fsr_dret2 < 0 || in.fsr_dret2[k] < l.fsr_dret2) {
              l.fsr_pt = in.fsr_pt[k];
              l.fsr_eta = in.fsr_eta[k];
              l.fsr_phi = in.fsr_phi[k];
              l.fsr_dret2 = in.fsr_dret2[k];
            }
          }
          muons.push_back(l);
        }
        std::vector<bool> ghost(muons.size(), false);
        for (std::size_t a = 0; a < muons.size(); ++a)
          for (std::size_t b = a + 1; b < muons.size(); ++b) {
            if (ghost[a] || ghost[b]) continue;
            if (h4l::delta_r(muons[a].eta, muons[a].phi, muons[b].eta, muons[b].phi) >= cfg.ghost_dr) continue;
            const bool a_pf = muons[a].flags >> 3 & 1U, b_pf = muons[b].flags >> 3 & 1U;
            const bool keep_a = a_pf != b_pf ? a_pf : muons[a].pt >= muons[b].pt;
            ghost[keep_a ? b : a] = true;
          }
        std::vector<Lepton> leptons;
        for (std::size_t a = 0; a < muons.size(); ++a)
          if (!ghost[a]) leptons.push_back(muons[a]);
        // Loose electrons (calibrated).
        for (std::size_t i = 0; i < in.el_pt.GetSize(); ++i) {
          if (std::fabs(in.el_eta[i]) >= cfg.el_abs_eta || std::fabs(in.el_dxy[i]) >= cfg.abs_dxy ||
              std::fabs(in.el_dz[i]) >= cfg.abs_dz)
            continue;
          Lepton l;
          l.pdg = -11 * in.el_charge[i];
          l.charge = in.el_charge[i];
          l.index = static_cast<int>(i);
          l.pt_raw = in.el_pt[i];
          l.eta = in.el_eta[i];
          l.eta_sc = in.el_eta[i] + in.el_deta_sc[i];
          const double abs_sc = std::fabs(l.eta_sc);
          l.normal = is_mc ? h4l::object_normal(key, static_cast<std::uint64_t>(entry), h4l::Collection::Electron, i, 0) : 0.0;
          l.pt = is_mc ? electron_payload.smeared_mc_pt(l.pt_raw, abs_sc, l.normal)
                       : electron_payload.corrected_data_pt(l.pt_raw, abs_sc);
          if (l.pt <= cfg.el_floor_pt) continue;
          l.phi = in.el_phi[i];
          l.rel_err = in.el_energy_err[i] / (in.el_pt[i] * std::cosh(in.el_eta[i]));
          l.iso = in.el_iso[i];
          l.iso_chg = in.el_iso_chg[i];
          l.sip = in.el_sip[i];
          l.dxy = in.el_dxy[i];
          l.dz = in.el_dz[i];
          l.mva_noiso = in.el_mva_noiso[i];
          l.mva_iso = in.el_mva_iso[i];
          const bool bits[] = {in.el_wpl_noiso[i], in.el_wp90_noiso[i], in.el_wp80_noiso[i], in.el_wpl_iso[i],
                               in.el_wp90_iso[i],  in.el_wp80_iso[i],   in.el_convveto[i]};
          for (std::size_t b = 0; b < std::size(bits); ++b) l.flags |= static_cast<unsigned>(bits[b]) << b;
          l.cut_based = in.el_cutbased[i];
          l.lost_hits = in.el_losthits[i];
          const auto [eta_hlt, phi_hlt] = h4l::electron_hlt_position(l.eta_sc, l.phi, l.charge, l.pt_raw, cfg.bending);
          const auto match = h4l::match_trigger(electron_objects, eta_hlt, phi_hlt, cfg.match_dr);
          l.trig_bits = match.bits;
          l.trig_pt = match.pt;
          if (is_mc) {
            l.gen_flav = (*in.el_gen_flav)[i];
            const int g = (*in.el_gen_idx)[i];
            if (g >= 0 && static_cast<std::size_t>(g) < in.gen_pt->GetSize()) l.gen_pt = (*in.gen_pt)[g];
          }
          leptons.push_back(l);
        }
        if (static_cast<int>(leptons.size()) < cfg.min_leptons) continue;
        // FSR-subtracted isolation with the FSR photons of every loose muon.
        for (auto& l : leptons) {
          double photons = 0;
          for (const auto& owner : leptons) {
            if (!(owner.fsr_pt > 0)) continue;
            const double dr = h4l::delta_r(l.eta, l.phi, owner.fsr_eta, owner.fsr_phi);
            if (dr > cfg.iso_veto_dr && dr < cfg.iso_dr) photons += owner.fsr_pt;
          }
          l.iso_fsr = l.iso_chg + std::max(0.0, (l.iso - l.iso_chg) - photons / l.pt_raw);
        }
        std::stable_sort(leptons.begin(), leptons.end(), [](const Lepton& a, const Lepton& b) { return a.pt > b.pt; });
        row.nlep_all = static_cast<Int_t>(leptons.size());
        if (static_cast<int>(leptons.size()) > kMaxLeptons) {
          ++lepton_overflow;
          leptons.resize(kMaxLeptons);
        }
        // Electron-muon overlap (the cross cleaning is decided downstream).
        for (auto& l : leptons) {
          if (l.muon()) continue;
          double best = cfg.overlap_dr;
          for (std::size_t m = 0; m < leptons.size(); ++m) {
            if (!leptons[m].muon()) continue;
            const double dr = h4l::delta_r(l.eta, l.phi, leptons[m].eta, leptons[m].phi);
            if (dr < cfg.overlap_dr) l.overlap_mask |= 1U << m;
            if (dr < best) {
              best = dr;
              l.overlap = static_cast<int>(m);
            }
          }
        }
        // Same-flavour pairs.
        struct Pair {
          int a, b;
          double mass;
          bool os;
        };
        std::vector<Pair> pairs;
        bool z_like = false;
        const int n = static_cast<int>(leptons.size());
        for (int a = 0; a < n; ++a)
          for (int b = a + 1; b < n; ++b) {
            if (std::abs(leptons[a].pdg) != std::abs(leptons[b].pdg)) continue;
            const bool os = leptons[a].charge * leptons[b].charge < 0;
            const double mass = (leptons[a].p4_fsr() + leptons[b].p4_fsr()).M();
            if (os && mass > cfg.event_z_low && mass < cfg.event_z_high) z_like = true;
            if (leptons[a].sip >= cfg.cand_sip || leptons[b].sip >= cfg.cand_sip) continue;
            if (mass <= cfg.z2_low || mass >= cfg.z2_high) continue;
            pairs.push_back({a, b, mass, os});
          }
        if (!z_like) continue;
        std::vector<Candidate> candidates;
        int overflow = 0;
        for (const auto& z1 : pairs) {
          if (!z1.os || z1.mass <= cfg.z1_low || z1.mass >= cfg.z1_high) continue;
          for (const auto& z2 : pairs) {
            if (z2.a == z1.a || z2.a == z1.b || z2.b == z1.a || z2.b == z1.b) continue;
            const int idx[4] = {z1.a, z1.b, z2.a, z2.b};
            double pts[4];
            bool pass = true;
            for (int u = 0; u < 4 && pass; ++u) {
              pts[u] = leptons[idx[u]].pt;
              for (int v = u + 1; v < 4 && pass; ++v) {
                const Lepton &x = leptons[idx[u]], &y = leptons[idx[v]];
                if (h4l::delta_r(x.eta, x.phi, y.eta, y.phi) <= cfg.min_dr) pass = false;
                if (x.charge * y.charge < 0 && (x.p4() + y.p4()).M() <= cfg.min_os_mass) pass = false;
              }
            }
            if (!pass) continue;
            std::sort(pts, pts + 4, std::greater<double>());
            if (pts[0] <= cfg.lead_pt || pts[1] <= cfg.sublead_pt) continue;
            const double m4l = (leptons[z1.a].p4_fsr() + leptons[z1.b].p4_fsr() + leptons[z2.a].p4_fsr() +
                                leptons[z2.b].p4_fsr()).M();
            if (m4l <= cfg.min_m4l) continue;
            Candidate c;
            c.leptons = {z1.a, z1.b, z2.a, z2.b};
            c.m_z1 = z1.mass;
            c.m_z2 = z2.mass;
            c.m4l = m4l;
            c.z2_same_sign = !z2.os;
            const double d1 = std::fabs(z1.mass - h4l::kZMass), d2 = std::fabs(z2.mass - h4l::kZMass);
            c.z1_closer = !z2.os || d1 < d2 || (d1 == d2 && z1.a < z2.a);
            // Alternative pairing of 4e / 4mu with two opposite-sign pairs (smart cut): each Z1 lepton
            // with the Z2 lepton of opposite charge.
            if (z2.os && std::abs(leptons[z1.a].pdg) == std::abs(leptons[z2.a].pdg)) {
              const int partner_a = leptons[z2.a].charge != leptons[z1.a].charge ? z2.a : z2.b;
              const int partner_b = partner_a == z2.a ? z2.b : z2.a;
              const double m1 = (leptons[z1.a].p4_fsr() + leptons[partner_a].p4_fsr()).M();
              const double m2 = (leptons[z1.b].p4_fsr() + leptons[partner_b].p4_fsr()).M();
              const bool first_closer = std::fabs(m1 - h4l::kZMass) < std::fabs(m2 - h4l::kZMass);
              c.m_za = first_closer ? m1 : m2;
              c.m_zb = first_closer ? m2 : m1;
            }
            if (static_cast<int>(candidates.size()) < kMaxCandidates) candidates.push_back(c);
            else ++overflow;
          }
        }
        // Fill the row.
        row.file_key = key;
        row.entry = entry;
        row.weight = is_mc ? **in.gen_weight : 1.f;
        row.met = *in.met_pt;
        row.met_phi = *in.met_phi;
        row.npv = *in.npv_good;
        row.rho = *in.rho;
        row.pu_true = is_mc ? **in.pu_true : -1.f;
        row.trig_bits = 0;
        for (std::size_t p = 0; p < in.paths.size(); ++p)
          if (**in.paths[p]) row.trig_bits |= 1U << p;
        row.trigger = (row.trig_bits & analysis_mask) != 0;
        row.nlep = n;
        for (int k = 0; k < n; ++k) {
          const Lepton& l = leptons[k];
          row.l_pdg[k] = l.pdg;
          row.l_index[k] = l.index;
          row.l_overlap[k] = l.overlap;
          row.l_cut_based[k] = l.cut_based;
          row.l_lost_hits[k] = l.lost_hits;
          row.l_gen_flav[k] = l.gen_flav;
          row.l_trig_bits[k] = l.trig_bits;
          row.l_flags[k] = static_cast<UShort_t>(l.flags);
          row.l_overlap_mask[k] = static_cast<UShort_t>(l.overlap_mask);
          row.l_iso_chg[k] = static_cast<Float_t>(l.iso_chg);
          row.l_pt[k] = static_cast<Float_t>(l.pt);
          row.l_pt_raw[k] = static_cast<Float_t>(l.pt_raw);
          row.l_eta[k] = static_cast<Float_t>(l.eta);
          row.l_eta_sc[k] = static_cast<Float_t>(l.eta_sc);
          row.l_phi[k] = static_cast<Float_t>(l.phi);
          row.l_normal[k] = static_cast<Float_t>(l.normal);
          row.l_rel_err[k] = static_cast<Float_t>(l.rel_err);
          row.l_iso[k] = static_cast<Float_t>(l.iso);
          row.l_iso_fsr[k] = static_cast<Float_t>(l.iso_fsr);
          row.l_sip[k] = static_cast<Float_t>(l.sip);
          row.l_dxy[k] = static_cast<Float_t>(l.dxy);
          row.l_dz[k] = static_cast<Float_t>(l.dz);
          row.l_mva_noiso[k] = static_cast<Float_t>(l.mva_noiso);
          row.l_mva_iso[k] = static_cast<Float_t>(l.mva_iso);
          row.l_gen_pt[k] = static_cast<Float_t>(l.gen_pt);
          row.l_trig_pt[k] = static_cast<Float_t>(l.trig_pt);
          row.l_fsr_pt[k] = static_cast<Float_t>(l.fsr_pt);
          row.l_fsr_eta[k] = static_cast<Float_t>(l.fsr_eta);
          row.l_fsr_phi[k] = static_cast<Float_t>(l.fsr_phi);
          row.l_fsr_dret2[k] = static_cast<Float_t>(l.fsr_dret2);
        }
        row.njet = 0;
        row.njet_all = 0;
        for (std::size_t j = 0; j < in.jet_pt.GetSize(); ++j) {
          if (in.jet_pt[j] <= cfg.jet_pt || std::fabs(in.jet_eta[j]) >= cfg.jet_eta || !(in.jet_id[j] & cfg.jet_id_mask)) continue;
          ++row.njet_all;
          if (row.njet >= kMaxJets) {
            ++jet_overflow;
            continue;
          }
          const int k = row.njet++;
          row.j_pt[k] = in.jet_pt[j];
          row.j_eta[k] = in.jet_eta[j];
          row.j_phi[k] = in.jet_phi[j];
          row.j_mass[k] = in.jet_mass[j];
          row.j_csv[k] = in.jet_csv[j];
          row.j_deepjet[k] = in.jet_deepjet[j];
          row.j_id[k] = in.jet_id[j];
          row.j_puid[k] = in.jet_puid[j];
          row.j_flavour[k] = is_mc ? (*in.jet_hadron_flavour)[j] : -1;
          row.j_gen[k] = is_mc ? (*in.jet_gen_idx)[j] : -1;
        }
        row.ncand = static_cast<Int_t>(candidates.size());
        row.ncand_overflow = overflow;
        for (int k = 0; k < row.ncand; ++k) {
          const Candidate& c = candidates[k];
          row.c_l1[k] = static_cast<Char_t>(c.leptons[0]);
          row.c_l2[k] = static_cast<Char_t>(c.leptons[1]);
          row.c_l3[k] = static_cast<Char_t>(c.leptons[2]);
          row.c_l4[k] = static_cast<Char_t>(c.leptons[3]);
          row.c_ss[k] = c.z2_same_sign;
          row.c_z1_closer[k] = c.z1_closer;
          row.c_mz1[k] = static_cast<Float_t>(c.m_z1);
          row.c_mz2[k] = static_cast<Float_t>(c.m_z2);
          row.c_m4l[k] = static_cast<Float_t>(c.m4l);
          row.c_mza[k] = static_cast<Float_t>(c.m_za);
          row.c_mzb[k] = static_cast<Float_t>(c.m_zb);
        }
        // Theory weights of MC events with a candidate.
        row.n_scale = row.n_pdf = row.n_ps = 0;
        if (is_mc && row.ncand > 0) {
          auto copy = [&](const std::unique_ptr<Array<Float_t>>& from, Float_t* to, Int_t& count, int capacity) {
            if (!from) return;
            const int size = static_cast<int>(from->GetSize());
            if (size > capacity) throw std::runtime_error("theory weight vector longer than its buffer in " + input_path);
            count = size;
            for (int k = 0; k < size; ++k) to[k] = (*from)[k];
          };
          copy(in.lhe_scale, row.w_scale, row.n_scale, kMaxScale);
          copy(in.lhe_pdf, row.w_pdf, row.n_pdf, kMaxPdf);
          copy(in.ps, row.w_ps, row.n_ps, kMaxPs);
          const int lengths[3] = {row.n_scale, row.n_pdf, row.n_ps};
          for (int w = 0; w < 3; ++w) {
            length_min[w] = std::min(length_min[w], lengths[w]);
            length_max[w] = std::max(length_max[w], lengths[w]);
          }
        }
        row.g_present = false;
        row.g = GenFields{};
        if (is_signal) {
          const auto found = gen.find({key, entry});
          if (found == gen.end()) throw InputError("event without a generator-table row in " + input_path);
          row.g_present = true;
          row.g = found->second;
        }
        check_read(input_path + " entry " + std::to_string(e));
        output->cd();
        if (events_out.Fill() <= 0) throw std::runtime_error("Events4l Fill failed");
        ++input_kept;
        if (row.ncand > 0) ++input_with_candidate;
        input_candidates += row.ncand;
        candidate_overflow += overflow;
      }
      }  // the reader is destroyed before its file
      check_read(input_path);
      {
        int w = 0;
        for (const char* name : {"LHEScaleWeight", "LHEPdfWeight", "PSWeight"}) {
          theory[name]["copied_length_min"] = length_max[w] < 0 ? json() : json(length_min[w]);
          theory[name]["copied_length_max"] = length_max[w] < 0 ? json() : json(length_max[w]);
          ++w;
        }
      }
      total_events += input_events;
      total_kept += input_kept;
      total_with_candidate += input_with_candidate;
      total_candidates += input_candidates;
      if (inputs_tree.Fill() <= 0) throw std::runtime_error("Inputs Fill failed");
      inputs_report.push_back({{"root", input_path}, {"skim_task_id", input.at("skim_task_id")},
                               {"root_sha256", input.at("root_sha256")}, {"events", input_events},
                               {"gen_rows", input.value("gen_rows", json())}, {"kept", input_kept},
                               {"with_candidate", input_with_candidate}, {"candidates", input_candidates},
                               {"theory_branches", theory}});
      std::cout << "[reco] " << input_path << ": " << input_events << " slim events, " << input_kept << " kept, "
                << input_with_candidate << " with a candidate, " << input_candidates << " candidates\n";
      file.reset();
      check_read(input_path + " (close)");
    }
    output->cd();
    for (TTree* tree : {&events_out, &inputs_tree})
      if (tree->Write("", TObject::kOverwrite) <= 0) throw std::runtime_error(std::string("cannot write ") + tree->GetName());
    output->Close();
    output.reset();
    if (g_read_error) throw std::runtime_error("ROOT error while writing: " + g_read_message);
    const std::map<std::string, long long> trees = {{"Events4l", total_kept},
                                                    {"Inputs", static_cast<long long>(task.at("inputs").size())}};
    h4l::validate_root_output_full(out_root, trees);
    if (g_read_error) throw std::runtime_error("ROOT error while validating the output: " + g_read_message);
    json report = {{"schema", "h4l_v3_reco_report/1"},
                   {"task_id", task.at("task_id")},
                   {"sample", task.at("config").at("sample")},
                   {"kind", task.at("config").at("kind")},
                   {"role", task.at("config").at("role")},
                   {"reco_config_version", task.at("reco_config").at("version")},
                   {"reco_config_fnv1a64", hex64(h4l::fnv1a64(task.at("reco_config").dump()))},
                   {"triggers_fnv1a64", hex64(h4l::fnv1a64(task.at("triggers").dump()))},
                   {"inputs_root_sha256_origin", "planned"},
                   {"calibration_payload_fnv1a64", hex64(h4l::fnv1a64(task.at("calibration_payload").dump()))},
                   {"trigger_paths", trigger_paths},
                   {"analysis_trigger_mask", analysis_mask},
                   {"frozen_program", task.value("frozen_program", json())},
                   {"inputs", inputs_report},
                   {"source_file_keys", task.at("source_file_keys")},
                   {"totals", {{"events", total_events}, {"kept", total_kept}, {"with_candidate", total_with_candidate},
                               {"candidates", total_candidates}, {"candidate_overflow", candidate_overflow},
                               {"lepton_overflow", lepton_overflow}, {"jet_overflow", jet_overflow}}},
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
