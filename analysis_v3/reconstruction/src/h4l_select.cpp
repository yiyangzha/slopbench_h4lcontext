// h4l_select: the final H -> ZZ* -> 4l selection (AN-16-442 sections 4-6)
// on the 4-lepton event records of h4l_reco.
//
// Usage: h4l_select --task TASK.json --out-json OUT.json --out-root OUT.root
//        h4l_select --mela-selftest OUT.json --sqrts TEV --mh GEV
//
// The task embeds the selection configuration (selection_ul16_v2.json), the
// lepton calibration payload (applied to the stored raw pT of every lepton:
// data pT exp(-u), MC pT (1 + r N) with the stored per-lepton deviate, so a
// new calibration needs no new event records; with the payload of the event
// records the stored pT is reproduced, which is checked), the lambda payload (per-event mass
// uncertainty), the refit inputs (true Z1 line shape, FSR photon resolution)
// and the MELA block: the sha256 of every MELA library, dictionary and data
// file (verified before MELA is initialized), the data files MELA opens for
// update (copied into the empty working directory first, so that the shared
// installation is only read) and the reference values of a fixed MELA
// self-test (six candidates: 4mu, 4e, 2e2mu without jets, one 1-jet event, a
// VBF-like and a VH-like 2-jet event), which every task repeats after the
// initialization and must
// reproduce to 1e-6 (a MELA constant file that failed to load would give
// valid-looking but wrong discriminants).  The --mela-selftest mode writes
// those reference values.  Per event (trigger OR and a good vertex
// required):
//   * leptons: the AN loose lepton (pT > 5 / 7 GeV; a tracker-only muon
//     within dR < 0.05 of a same-charge PF muon is a ghost; a lepton without
//     finite kinematics, momentum error, isolation and SIP is not loose), the AN
//     FSR-subtracted isolation (3.3: photons of loose muons passing SIP are
//     removed from the isolation of loose leptons passing SIP; veto
//     dR > 0.01 for muons, |eta_SC| < 1.479 or dR > 0.08 for electrons), the
//     selected lepton (loose + tight ID + SIP < 4 + isolation < 0.35; tight
//     muon = PF, or tracker high-pT above 200 GeV; the N-1 variants drop it with
//     leptons.muon.require_id = false or leptons.electron.id = "none"), and the cross cleaning
//     (an electron within dR < 0.05 of any selected muon is removed);
//   * SR and SRZ4l (Z2 > 12 or 4 GeV): candidates of four selected leptons,
//     Z1 the pair closer to m_Z, 40 < m_Z1 < 120, z2 < m_Z2 < 120, pT 20 / 10,
//     m(OS) > 4 without FSR for all opposite-sign pairs, dR > 0.02, m4l > 70,
//     smart cut; of the candidates built from the same four leptons the one
//     with Z1 closest to m_Z, then the highest D_bkg^kin (MELA);
//   * CR (Z+X), for every signal region (region field; events of that signal
//     region excluded): 2P2F and 3P1F (Z1 of selected leptons, Z2 an
//     opposite-sign pair of loose leptons with SIP < 4 of which two or one
//     fail the tight ID and isolation, the region's kinematic cuts, best by
//     D_bkg^kin per type) and SS (Z2 a same-sign pair of loose leptons with
//     SIP < 4, m4l > 70 GeV stored, the AN cut at 100 GeV applied downstream;
//     best by Z1 closest to m_Z then the largest Z2 scalar pT sum);
//   * ZL (fake rates): a Z1 of selected leptons with 40 < m_Z1 < 120 GeV
//     (closest to m_Z; the |m_Z1 - m_Z| < 7 GeV window of the fake rates and
//     the windows of the SS method are applied downstream) and pT 20 / 10
//     GeV, exactly one additional loose lepton (no ID, isolation or SIP
//     requirement; its SIP is stored for the loose + SIP denominator) with
//     m(probe, opposite-sign Z1 lepton) > 4 GeV; the 3-lepton mass for the
//     conversion correction.
// For every selected candidate: masses (with FSR), the per-event mass
// uncertainty (AN 5.3: quadrature of the m4l changes for each lepton moved by
// its lambda-corrected momentum error, FSR photons with the simulation
// parametrization), the Z1 kinematic refit (AN 5.4: Gaussian pT constraints
// times the true Z1 line shape (cubic interpolation), photons fixed; the
// covariance from the numerical Hessian of the likelihood), MELA
// probabilities and discriminants with the jets cleaned for that candidate
// (pT > 30, |eta| < 4.7, dR > 0.4 from the selected leptons, the candidate's
// leptons and the selected FSR photons; for a same-sign Z2 one charge is
// flipped for MELA), b tags, additional leptons, MET, the ID flags of the
// candidate leptons (tighter working points downstream) and, for MC, the
// responses of m4l to the lepton scale and smearing of
// each flavour.  Normalization and scale factors are applied downstream.

#include "h4l/calibration.h"
#include "h4l/hash.h"
#include "h4l/io.h"
#include "h4l/kinematics.h"
#include "h4l/root_io.h"
#include "h4l/sha256.h"

#include "MelaWrapper.h"

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <TError.h>
#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <unistd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <limits>
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

// Canonical text of a JSON value for the configuration hashes: sorted keys, no spaces, floating-point
// numbers as %.17g (the shortest round-trip form is ambiguous at ties, e.g. 66.587677001953125, where
// nlohmann and Python choose different last digits; scan_select.py applies the same rule).
std::string canonical(const json& value) {
  switch (value.type()) {
    case json::value_t::object: {
      std::string out = "{";
      bool first = true;
      for (const auto& [key, item] : value.items()) {
        out += (first ? "" : ",") + json(key).dump() + ":" + canonical(item);
        first = false;
      }
      return out + "}";
    }
    case json::value_t::array: {
      std::string out = "[";
      for (std::size_t i = 0; i < value.size(); ++i) out += (i ? "," : "") + canonical(value[i]);
      return out + "]";
    }
    case json::value_t::number_float: {
      char buffer[40];
      std::snprintf(buffer, sizeof(buffer), "%.17g", value.get<double>());
      return buffer;
    }
    default:
      return value.dump();
  }
}
std::string config_hash(const json& value) { return hex64(h4l::fnv1a64(canonical(value))); }

// The fixed MELA self-test: six candidates (dressed leptons Z1, Z1, Z2, Z2; jets) whose probabilities,
// constants and discriminants must be reproduced by every task.
json mela_selftest() {
  struct Case {
    const char* name;
    h4lmela::Lepton l[4];
    std::vector<h4lmela::Jet> jets;
  };
  const std::vector<Case> cases = {
      {"4mu_0jet", {{48.3, 0.41, 0.52, 0.1057, 13}, {31.7, -0.83, 2.91, 0.1057, -13}, {22.4, 1.21, -1.73, 0.1057, 13},
                    {12.9, -0.12, -2.64, 0.1057, -13}}, {}},
      {"4e_0jet", {{52.1, -0.35, 1.12, 0.000511, 11}, {27.8, 1.02, -2.18, 0.000511, -11}, {19.6, -1.44, -0.41, 0.000511, 11},
                   {14.2, 0.77, 2.37, 0.000511, -11}}, {}},
      {"2e2mu_0jet", {{44.9, 0.18, -0.95, 0.1057, 13}, {35.2, -1.31, 2.22, 0.1057, -13}, {17.3, 0.92, 0.61, 0.000511, 11},
                      {11.8, -0.54, -2.83, 0.000511, -11}}, {}},
      {"2e2mu_1jet", {{61.4, 0.64, 0.33, 0.000511, 11}, {29.5, -0.27, -2.71, 0.000511, -11}, {24.1, 1.58, 1.95, 0.1057, 13},
                      {10.6, -1.05, -1.22, 0.1057, -13}}, {{87.2, -2.13, 2.64, 9.8}}},
      {"4mu_2jet", {{39.7, 1.12, -0.21, 0.1057, 13}, {33.5, 0.08, 2.47, 0.1057, -13}, {21.9, -0.66, 1.36, 0.1057, 13},
                    {15.4, 1.73, -2.05, 0.1057, -13}}, {{104.6, 2.71, 0.88, 12.3}, {63.8, -2.46, -2.31, 8.7}}},
      {"2e2mu_2jet_vh", {{55.3, -0.21, 2.02, 0.1057, 13}, {30.8, 0.95, -0.47, 0.1057, -13}, {20.7, -1.12, -2.36, 0.000511, 11},
                         {13.1, 0.36, 0.84, 0.000511, -11}}, {{60.0, 0.5, 1.0, 7.5}, {45.0, -0.3, 2.6, 6.2}}}};
  json out = json::object();
  for (const auto& c : cases) {
    h4lmela::MelaInput in;
    for (int i = 0; i < 4; ++i) in.leptons[i] = c.l[i];
    in.jets = c.jets;
    h4lmela::MelaOutput mo;
    const bool ok = h4lmela::compute(in, mo, true);
    out[c.name] = {{"ok", ok}, {"nJets", mo.nJets}, {"p_sig_ggH", mo.p_sig_ggH}, {"p_bkg_qqZZ", mo.p_bkg_qqZZ},
                   {"p_JJVBF_sig", mo.p_JJVBF_sig}, {"p_JJQCD_sig", mo.p_JJQCD_sig}, {"p_HadWH_sig", mo.p_HadWH_sig},
                   {"p_HadZH_sig", mo.p_HadZH_sig}, {"p_JQCD_sig", mo.p_JQCD_sig}, {"p_JVBF_sig", mo.p_JVBF_sig},
                   {"pAux_JVBF_sig", mo.pAux_JVBF_sig}, {"p_HadWH_mavjj", mo.p_HadWH_mavjj},
                   {"p_HadZH_mavjj", mo.p_HadZH_mavjj}, {"c_bkg_kin", mo.c_bkg_kin}, {"c_2jet", mo.c_2jet},
                   {"c_1jet", mo.c_1jet}, {"c_WH", mo.c_WH}, {"c_ZH", mo.c_ZH}, {"D_bkg_kin", mo.D_bkg_kin},
                   {"D_2jet", mo.D_2jet}, {"D_1jet", mo.D_1jet}, {"D_WH", mo.D_WH}, {"D_ZH", mo.D_ZH}};
  }
  return out;
}

// Compare a self-test with its reference (relative 1e-6, absolute 1e-300 for vanishing values).
std::string selftest_mismatch(const json& got, const json& reference) {
  if (!reference.is_object() || reference.empty()) return "no reference values in the task";
  for (const auto& [name, values] : reference.items()) {
    if (!got.contains(name)) return name + " missing";
    for (const auto& [key, value] : values.items()) {
      const json& g = got.at(name).at(key);
      if (value.is_boolean() || value.is_number_integer()) {
        if (g != value) return name + "." + key + " differs";
        continue;
      }
      const double a = g.get<double>(), b = value.get<double>();
      if (!std::isfinite(a) || std::fabs(a - b) > 1e-6 * std::max(std::fabs(a), std::fabs(b)) + 1e-300)
        return name + "." + key + ": " + std::to_string(a) + " against the reference " + std::to_string(b);
    }
  }
  return "";
}

constexpr double kZ = h4l::kZMass;
template <typename T>
using Array = TTreeReaderArray<T>;
template <typename T>
using Value = TTreeReaderValue<T>;

double upper(const json& range) {
  return range.at(1).is_null() ? std::numeric_limits<double>::infinity() : range.at(1).get<double>();
}

// ---------------------------------------------------------------- configuration
struct Config {
  double mu_pt = 5, el_pt = 7, mu_eta = 2.4, el_eta = 2.5, max_sip = 4, max_iso = 0.35, muon_high_pt = 200,
         cross_dr = 0.05, ghost_dr = 0.05, ss_min_m4l = 100;
  int electron_id_bit = 1;  // -1: no electron identification requirement (N-1 variant)
  bool muon_id = true;       // false: no muon identification requirement (N-1 variant)
  double z1_low = 40, z1_high = 120, z2_high = 120, lead_pt = 20, sublead_pt = 10, min_dr = 0.02, min_os_mass = 4,
         min_m4l = 70;
  std::vector<std::pair<std::string, double>> regions;  // name, z2 lower bound
  double zl_z1_low = 40, zl_z1_high = 120, zl_lead = 20, zl_sublead = 10, zl_probe_os_mass = 4;
  double jet_pt = 30, jet_eta = 4.7, jet_clean_dr = 0.4, btag_threshold = 0.2489, btag_eta = 2.4;
  int min_good_pv = 1;
};

Config read_config(const json& s) {
  Config c;
  const json& l = s.at("leptons");
  c.mu_pt = l.at("muon").at("pt").get<double>();
  c.mu_eta = l.at("muon").at("abs_eta").get<double>();
  c.muon_high_pt = l.at("muon").at("high_pt_alternative_above").get<double>();
  c.el_pt = l.at("electron").at("pt").get<double>();
  c.el_eta = l.at("electron").at("abs_eta").get<double>();
  const std::string id = l.at("electron").at("id").get<std::string>();
  const std::map<std::string, int> bits = {{"mvaFall17V2noIso_WPL", 0}, {"mvaFall17V2noIso_WP90", 1},
                                           {"mvaFall17V2noIso_WP80", 2}, {"mvaFall17V2Iso_WPL", 3},
                                           {"mvaFall17V2Iso_WP90", 4},   {"mvaFall17V2Iso_WP80", 5}};
  if (id == "none") {
    c.electron_id_bit = -1;
  } else {
    if (!bits.count(id)) throw std::runtime_error("unknown electron id " + id);
    c.electron_id_bit = bits.at(id);
  }
  c.muon_id = l.at("muon").value("require_id", true);
  c.max_sip = l.at("max_sip").get<double>();
  c.max_iso = l.at("max_iso_fsr").get<double>();
  c.cross_dr = l.at("cross_clean_dr").get<double>();
  c.ghost_dr = l.at("ghost_tracker_muon_dr").get<double>();
  const json& k = s.at("candidate");
  c.z1_low = k.at("z1_mass").at(0).get<double>();
  c.z1_high = k.at("z1_mass").at(1).get<double>();
  c.z2_high = k.at("z2_high").get<double>();
  c.lead_pt = k.at("leading_pt").get<double>();
  c.sublead_pt = k.at("subleading_pt").get<double>();
  c.min_dr = k.at("min_dr").get<double>();
  c.min_os_mass = k.at("min_os_mass").get<double>();
  c.min_m4l = k.at("min_m4l").get<double>();
  c.ss_min_m4l = k.at("ss_min_m4l").get<double>();
  for (const auto& r : k.at("regions")) c.regions.push_back({r.at("name").get<std::string>(), r.at("z2_low").get<double>()});
  if (c.regions.empty() || c.regions.front().first != "SR") throw std::runtime_error("the first region must be SR");
  const json& z = s.at("zl");
  c.zl_z1_low = z.at("z1_mass").at(0).get<double>();
  c.zl_z1_high = z.at("z1_mass").at(1).get<double>();
  c.zl_lead = z.at("leading_pt").get<double>();
  c.zl_sublead = z.at("subleading_pt").get<double>();
  c.zl_probe_os_mass = z.at("min_probe_os_mass").get<double>();
  const json& j = s.at("jets");
  c.jet_pt = j.at("pt").get<double>();
  c.jet_eta = j.at("abs_eta").get<double>();
  c.jet_clean_dr = j.at("clean_dr").get<double>();
  c.btag_threshold = j.at("btag_threshold").get<double>();
  c.btag_eta = j.at("btag_abs_eta").get<double>();
  c.min_good_pv = s.at("event").at("min_good_pv").get<int>();
  return c;
}

// Per-lepton momentum-error correction (lambda payload).
struct LambdaRegion {
  double eta_low = 0, eta_high = 0, err_low = 0, err_high = 0, value = 1;
};
std::vector<LambdaRegion> read_lambda(const json& flavour, const std::string& role) {
  std::vector<LambdaRegion> out;
  for (const auto& entry : flavour.at("regions")) {
    const json& region = entry.at("region");
    LambdaRegion r;
    r.eta_low = region.at("abs_eta").at(0).get<double>();
    r.eta_high = upper(region.at("abs_eta"));
    if (region.contains("rel_err")) {
      r.err_low = region.at("rel_err").at(0).get<double>();
      r.err_high = upper(region.at("rel_err"));
    } else {
      r.err_high = std::numeric_limits<double>::infinity();
    }
    r.value = entry.at(role).at("lambda").get<double>();
    out.push_back(r);
  }
  return out;
}
double lambda_of(const std::vector<LambdaRegion>& regions, double abs_eta, double rel_err) {
  for (const auto& r : regions)
    if (abs_eta >= r.eta_low && abs_eta < r.eta_high && rel_err >= r.err_low && rel_err < r.err_high) return r.value;
  throw std::runtime_error("lepton outside every lambda region");
}

// Refit inputs: the true Z1 line shape and the FSR photon resolution.
struct Lineshape {
  double lo = 0, width = 0;
  std::vector<double> log_density;
  // Catmull-Rom cubic interpolation of the log density (continuous first derivative, so that the
  // curvature the refit covariance needs is defined); steep linear fall-off outside the table.
  double log_at(double m) const {
    const double x = (m - lo) / width - 0.5;
    const int n = static_cast<int>(log_density.size());
    if (x <= 0) return log_density.front() - 50.0 * std::max(0.0, -x * width);
    if (x >= n - 1) return log_density.back() - 50.0 * std::max(0.0, (x - (n - 1)) * width);
    const int i = static_cast<int>(x);
    const double t = x - i;
    const double p0 = log_density[std::max(i - 1, 0)], p1 = log_density[i], p2 = log_density[i + 1],
                 p3 = log_density[std::min(i + 2, n - 1)];
    return 0.5 * (2 * p1 + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t * t + (-p0 + 3 * p1 - 3 * p2 + p3) * t * t * t);
  }
};

// ---------------------------------------------------------------- input rows

struct Reader {
  TTreeReader reader;
  Value<ULong64_t> file_key;
  Value<Long64_t> entry;
  Value<Float_t> weight, met, met_phi;
  Value<Int_t> npv, nlep, njet, ncand, nlep_all, njet_all, ncand_overflow;
  Value<UInt_t> trig_bits;
  Value<Bool_t> trigger;
  Array<Int_t> l_pdg, l_index, l_lost_hits, l_gen_flav, l_trig_bits;
  Array<UShort_t> l_flags, l_overlap_mask;
  Array<Float_t> l_pt, l_pt_raw, l_eta, l_eta_sc, l_phi, l_normal, l_rel_err, l_iso, l_iso_chg, l_sip, l_mva_noiso, l_gen_pt,
      l_trig_pt, l_fsr_pt, l_fsr_eta, l_fsr_phi;
  Array<Float_t> j_pt, j_eta, j_phi, j_mass, j_deepjet;
  Array<Int_t> j_flavour;
  Array<Char_t> c_l1, c_l2, c_l3, c_l4;
  Array<Bool_t> c_ss;
  Array<Float_t> c_mz1, c_mz2, c_m4l;
  Value<Int_t> n_scale, n_pdf, n_ps;
  Array<Float_t> w_scale, w_pdf, w_ps;
  Value<Bool_t> g_present, g_fid_pass;
  Value<Int_t> g_stage0, g_vh_class, g_final_state_true, g_fid_final_state, g_fid_njets;
  Value<Float_t> g_h_pt, g_h_eta, g_h_phi, g_h_mass, g_h_y, g_fid_m4l, g_fid_mz1, g_fid_mz2, g_fid_pt4l, g_fid_y4l,
      g_fid_jet1_pt;
  explicit Reader(TTree* t)
      : reader(t), file_key(reader, "file_key"), entry(reader, "entry"), weight(reader, "weight"), met(reader, "met"),
        met_phi(reader, "met_phi"), npv(reader, "npv"), nlep(reader, "nlep"), njet(reader, "njet"), ncand(reader, "ncand"),
        nlep_all(reader, "nlep_all"), njet_all(reader, "njet_all"), ncand_overflow(reader, "ncand_overflow"),
        trig_bits(reader, "trig_bits"), trigger(reader, "trigger"), l_pdg(reader, "l_pdg"), l_index(reader, "l_index"),
        l_lost_hits(reader, "l_lost_hits"), l_gen_flav(reader, "l_gen_flav"),
        l_trig_bits(reader, "l_trig_bits"), l_flags(reader, "l_flags"), l_overlap_mask(reader, "l_overlap_mask"), l_pt(reader, "l_pt"), l_pt_raw(reader, "l_pt_raw"),
        l_eta(reader, "l_eta"), l_eta_sc(reader, "l_eta_sc"), l_phi(reader, "l_phi"), l_normal(reader, "l_normal"),
        l_rel_err(reader, "l_rel_err"), l_iso(reader, "l_iso"), l_iso_chg(reader, "l_iso_chg"), l_sip(reader, "l_sip"),
        l_mva_noiso(reader, "l_mva_noiso"), l_gen_pt(reader, "l_gen_pt"), l_trig_pt(reader, "l_trig_pt"), l_fsr_pt(reader, "l_fsr_pt"),
        l_fsr_eta(reader, "l_fsr_eta"), l_fsr_phi(reader, "l_fsr_phi"), j_pt(reader, "j_pt"), j_eta(reader, "j_eta"),
        j_phi(reader, "j_phi"), j_mass(reader, "j_mass"), j_deepjet(reader, "j_deepjet"), j_flavour(reader, "j_flavour"),
        c_l1(reader, "c_l1"), c_l2(reader, "c_l2"), c_l3(reader, "c_l3"), c_l4(reader, "c_l4"), c_ss(reader, "c_ss"),
        c_mz1(reader, "c_mz1"), c_mz2(reader, "c_mz2"), c_m4l(reader, "c_m4l"), n_scale(reader, "n_scale"),
        n_pdf(reader, "n_pdf"), n_ps(reader, "n_ps"), w_scale(reader, "w_scale"), w_pdf(reader, "w_pdf"), w_ps(reader, "w_ps"),
        g_present(reader, "g_present"), g_fid_pass(reader, "g_fid_pass"), g_stage0(reader, "g_stage0"),
        g_vh_class(reader, "g_vh_class"), g_final_state_true(reader, "g_final_state_true"),
        g_fid_final_state(reader, "g_fid_final_state"), g_fid_njets(reader, "g_fid_njets"), g_h_pt(reader, "g_h_pt"),
        g_h_eta(reader, "g_h_eta"), g_h_phi(reader, "g_h_phi"), g_h_mass(reader, "g_h_mass"), g_h_y(reader, "g_h_y"),
        g_fid_m4l(reader, "g_fid_m4l"), g_fid_mz1(reader, "g_fid_mz1"), g_fid_mz2(reader, "g_fid_mz2"),
        g_fid_pt4l(reader, "g_fid_pt4l"), g_fid_y4l(reader, "g_fid_y4l"), g_fid_jet1_pt(reader, "g_fid_jet1_pt") {}
};

struct Lepton {
  int k = 0, pdg = 0, gen_flav = -1, lost_hits = 0;
  unsigned flags = 0;
  double pt = 0, pt_reco = 0, pt_raw = 0, eta = 0, eta_sc = 0, phi = 0, normal = 0, rel_err = 0, iso = 0, iso_all = 0, iso_chg = 0, sip = 0,
         gen_pt = -1, sigma = 0;
  unsigned overlap_mask = 0;
  double fsr_pt = -1, fsr_eta = 0, fsr_phi = 0, fsr_sigma = 0;
  bool loose = false, tight_id = false, selected = false, cleaned = false, ghost = false;
  bool muon() const { return std::abs(pdg) == 13; }
  int charge() const { return pdg > 0 ? -1 : 1; }
  double mass() const { return muon() ? h4l::kMuonMass : h4l::kElectronMass; }
};

P4 lepton_p4(const Lepton& l, double pt) { return P4(pt, l.eta, l.phi, l.mass()); }
P4 photon_p4(const Lepton& l) { return l.fsr_pt > 0 ? P4(l.fsr_pt, l.fsr_eta, l.fsr_phi, 0.0) : P4(); }
P4 dressed(const Lepton& l, double pt) { return lepton_p4(l, pt) + photon_p4(l); }

// Candidate kinematics recomputed from the leptons (pT may be varied).
struct Kinematics {
  double mz1 = 0, mz2 = 0, m4l = 0, mza = -1, mzb = -1;
  bool z1_closer = true, pass_cuts = false;
  P4 total;
};

struct Candidate {
  int index = -1;
  std::array<int, 4> leg{};
  bool ss = false;
};

// ---------------------------------------------------------------- output rows
struct CandidateRow {
  ULong64_t file_key = 0;
  Long64_t entry = 0;
  Float_t weight = 1, met = 0, met_phi = 0;
  Int_t npv = 0, cr_type = -1, final_state = -1, z1_flavour = 0, z2_flavour = 0, ncand_region = 0;
  UInt_t trig_bits = 0;
  Float_t m4l = 0, m4l_err = 0, m4l_refit = 0, m4l_refit_err = 0, mz1 = 0, mz2 = 0, mz1_refit = 0, pt4l = 0, eta4l = 0,
          phi4l = 0, y4l = 0;
  Int_t refit_status = -1;
  Int_t l_pdg[4], l_gen_flav[4], l_lost_hits[4], l_pass[4], l_index[4], l_trig_bits[4], l_flags[4];
  Float_t l_pt[4], l_pt_raw[4], l_eta[4], l_eta_sc[4], l_phi[4], l_sigma[4], l_iso[4], l_sip[4], l_mva[4], l_gen_pt[4],
      l_trig_pt[4], l_fsr_pt[4], l_fsr_eta[4], l_fsr_phi[4];
  Float_t p_sig_ggH = 0, p_bkg_qqZZ = 0, c_bkg_kin = 0, d_bkg_kin = -1, p_JJVBF = 0, p_JJQCD = 0, p_HadWH = 0, p_HadZH = 0,
          p_HadWH_mavjj = 0, p_HadWH_mavjj_true = 0, p_HadZH_mavjj = 0, p_HadZH_mavjj_true = 0, p_JQCD = 0, p_JVBF = 0,
          pAux_JVBF = 0, c_2jet = 0, c_1jet = 0, c_WH = 0, c_ZH = 0, d_2jet = -1, d_1jet = -1, d_wh_mela = -1, d_zh_mela = -1;
  Int_t mela_njets = 0, mela_charge_flipped = 0, njets = 0, nbjets = 0, n_extra_leptons = 0, n_extra_os_sf_pairs = 0,
        region = 0;
  Float_t jet_pt[4], jet_eta[4], jet_phi[4], jet_mass[4], jet_btag[4];
  Int_t jet_flavour[4];
  Float_t mjj = -1, detajj = -1;
  // MC: responses of m4l to the lepton scale (d m4l / d ln pT of every lepton of the flavour) and
  // smearing (d m4l / d r for an absolute change of the smear r of every lepton of the flavour); theory weights.
  Float_t dm4l_scale_mu = 0, dm4l_scale_e = 0, dm4l_smear_mu = 0, dm4l_smear_e = 0;
  Int_t n_scale = 0, n_pdf = 0, n_ps = 0;
  Float_t w_scale[16], w_pdf[128], w_ps[8];
  Bool_t g_present = false, g_fid_pass = false;
  Int_t g_stage0 = 0, g_vh_class = 0, g_final_state_true = -1, g_fid_final_state = -1, g_fid_njets = -1;
  Float_t g_h_pt = -1, g_h_eta = 0, g_h_phi = 0, g_h_mass = -1, g_h_y = 0, g_fid_m4l = -1, g_fid_mz1 = -1, g_fid_mz2 = -1,
          g_fid_pt4l = -1, g_fid_y4l = 0, g_fid_jet1_pt = -1;

  void book(TTree& t) {
    t.Branch("file_key", &file_key, "file_key/l");
    t.Branch("entry", &entry, "entry/L");
    t.Branch("trig_bits", &trig_bits, "trig_bits/i");
    for (auto [n, a] : std::initializer_list<std::pair<const char*, Float_t*>>{
             {"weight", &weight}, {"met", &met}, {"met_phi", &met_phi}, {"m4l", &m4l}, {"m4l_err", &m4l_err},
             {"m4l_refit", &m4l_refit}, {"m4l_refit_err", &m4l_refit_err}, {"mz1", &mz1}, {"mz2", &mz2},
             {"mz1_refit", &mz1_refit}, {"pt4l", &pt4l}, {"eta4l", &eta4l}, {"phi4l", &phi4l}, {"y4l", &y4l},
             {"p_sig_ggH", &p_sig_ggH}, {"p_bkg_qqZZ", &p_bkg_qqZZ}, {"c_bkg_kin", &c_bkg_kin}, {"d_bkg_kin", &d_bkg_kin},
             {"p_JJVBF", &p_JJVBF}, {"p_JJQCD", &p_JJQCD}, {"p_HadWH", &p_HadWH}, {"p_HadZH", &p_HadZH},
             {"p_HadWH_mavjj", &p_HadWH_mavjj}, {"p_HadWH_mavjj_true", &p_HadWH_mavjj_true},
             {"p_HadZH_mavjj", &p_HadZH_mavjj}, {"p_HadZH_mavjj_true", &p_HadZH_mavjj_true}, {"p_JQCD", &p_JQCD},
             {"p_JVBF", &p_JVBF}, {"pAux_JVBF", &pAux_JVBF}, {"c_2jet", &c_2jet}, {"c_1jet", &c_1jet}, {"c_WH", &c_WH},
             {"c_ZH", &c_ZH}, {"d_2jet", &d_2jet}, {"d_1jet", &d_1jet}, {"d_wh_mela", &d_wh_mela}, {"d_zh_mela", &d_zh_mela},
             {"mjj", &mjj}, {"detajj", &detajj}, {"dm4l_scale_mu", &dm4l_scale_mu}, {"dm4l_scale_e", &dm4l_scale_e},
             {"dm4l_smear_mu", &dm4l_smear_mu}, {"dm4l_smear_e", &dm4l_smear_e}, {"g_h_pt", &g_h_pt}, {"g_h_eta", &g_h_eta},
             {"g_h_phi", &g_h_phi}, {"g_h_mass", &g_h_mass}, {"g_h_y", &g_h_y}, {"g_fid_m4l", &g_fid_m4l},
             {"g_fid_mz1", &g_fid_mz1}, {"g_fid_mz2", &g_fid_mz2}, {"g_fid_pt4l", &g_fid_pt4l}, {"g_fid_y4l", &g_fid_y4l},
             {"g_fid_jet1_pt", &g_fid_jet1_pt}})
      t.Branch(n, a, (std::string(n) + "/F").c_str());
    for (auto [n, a] : std::initializer_list<std::pair<const char*, Int_t*>>{
             {"npv", &npv}, {"cr_type", &cr_type}, {"final_state", &final_state}, {"z1_flavour", &z1_flavour},
             {"refit_status", &refit_status},
             {"z2_flavour", &z2_flavour}, {"ncand_region", &ncand_region}, {"mela_njets", &mela_njets},
             {"mela_charge_flipped", &mela_charge_flipped}, {"region", &region}, {"njets", &njets},
             {"nbjets", &nbjets}, {"n_extra_leptons", &n_extra_leptons}, {"n_extra_os_sf_pairs", &n_extra_os_sf_pairs},
             {"n_scale", &n_scale}, {"n_pdf", &n_pdf}, {"n_ps", &n_ps}, {"g_stage0", &g_stage0}, {"g_vh_class", &g_vh_class},
             {"g_final_state_true", &g_final_state_true}, {"g_fid_final_state", &g_fid_final_state},
             {"g_fid_njets", &g_fid_njets}})
      t.Branch(n, a, (std::string(n) + "/I").c_str());
    t.Branch("g_present", &g_present, "g_present/O");
    t.Branch("g_fid_pass", &g_fid_pass, "g_fid_pass/O");
    for (auto [n, a] : std::initializer_list<std::pair<const char*, Int_t*>>{
             {"l_pdg", l_pdg}, {"l_gen_flav", l_gen_flav}, {"l_lost_hits", l_lost_hits}, {"l_pass", l_pass}, {"l_index", l_index},
             {"l_trig_bits", l_trig_bits}, {"l_flags", l_flags}})
      t.Branch(n, a, (std::string(n) + "[4]/I").c_str());
    for (auto [n, a] : std::initializer_list<std::pair<const char*, Float_t*>>{
             {"l_pt", l_pt}, {"l_pt_raw", l_pt_raw}, {"l_eta", l_eta}, {"l_eta_sc", l_eta_sc}, {"l_phi", l_phi},
             {"l_sigma", l_sigma}, {"l_iso", l_iso}, {"l_sip", l_sip}, {"l_mva", l_mva}, {"l_gen_pt", l_gen_pt},
             {"l_trig_pt", l_trig_pt}, {"l_fsr_pt", l_fsr_pt}, {"l_fsr_eta", l_fsr_eta}, {"l_fsr_phi", l_fsr_phi}})
      t.Branch(n, a, (std::string(n) + "[4]/F").c_str());
    for (auto [n, a] : std::initializer_list<std::pair<const char*, Float_t*>>{
             {"jet_pt", jet_pt}, {"jet_eta", jet_eta}, {"jet_phi", jet_phi}, {"jet_mass", jet_mass}, {"jet_btag", jet_btag}})
      t.Branch(n, a, (std::string(n) + "[4]/F").c_str());
    t.Branch("jet_flavour", jet_flavour, "jet_flavour[4]/I");
    t.Branch("w_scale", w_scale, "w_scale[n_scale]/F");
    t.Branch("w_pdf", w_pdf, "w_pdf[n_pdf]/F");
    t.Branch("w_ps", w_ps, "w_ps[n_ps]/F");
  }
};

struct ZlRow {
  ULong64_t file_key = 0;
  Long64_t entry = 0;
  Float_t weight = 1, met = 0, mz1 = 0, probe_pt = 0, probe_eta = 0, probe_eta_sc = 0, probe_phi = 0, probe_iso = 0,
          probe_sip = 0, probe_mva = 0, probe_os_mass = 0, tag1_pt = 0, tag2_pt = 0, m3l = 0, probe_dr_tag = 0;
  Int_t npv = 0, z1_flavour = 0, probe_pdg = 0, probe_gen_flav = -1, probe_lost_hits = 0, probe_id = 0, probe_pass = 0,
        n_loose = 0, n_loose_sip = 0, probe_flags = 0;
  UInt_t trig_bits = 0;
  void book(TTree& t) {
    t.Branch("file_key", &file_key, "file_key/l");
    t.Branch("entry", &entry, "entry/L");
    t.Branch("trig_bits", &trig_bits, "trig_bits/i");
    for (auto [n, a] : std::initializer_list<std::pair<const char*, Float_t*>>{
             {"weight", &weight}, {"met", &met}, {"mz1", &mz1}, {"probe_pt", &probe_pt}, {"probe_eta", &probe_eta},
             {"probe_eta_sc", &probe_eta_sc}, {"probe_phi", &probe_phi}, {"probe_iso", &probe_iso}, {"probe_sip", &probe_sip},
             {"probe_mva", &probe_mva}, {"probe_os_mass", &probe_os_mass}, {"tag1_pt", &tag1_pt}, {"tag2_pt", &tag2_pt},
             {"m3l", &m3l}, {"probe_dr_tag", &probe_dr_tag}})
      t.Branch(n, a, (std::string(n) + "/F").c_str());
    for (auto [n, a] : std::initializer_list<std::pair<const char*, Int_t*>>{
             {"npv", &npv}, {"z1_flavour", &z1_flavour}, {"probe_pdg", &probe_pdg}, {"probe_gen_flav", &probe_gen_flav},
             {"probe_lost_hits", &probe_lost_hits}, {"probe_id", &probe_id}, {"probe_pass", &probe_pass}, {"n_loose", &n_loose},
             {"n_loose_sip", &n_loose_sip}, {"probe_flags", &probe_flags}})
      t.Branch(n, a, (std::string(n) + "/I").c_str());
  }
};

// ---------------------------------------------------------------- physics helpers
Kinematics kinematics(const std::vector<Lepton>& leptons, const std::array<int, 4>& leg, const std::array<double, 4>& pt,
                      bool ss, const Config& cfg, double z2_low) {
  Kinematics k;
  const Lepton* l[4] = {&leptons[leg[0]], &leptons[leg[1]], &leptons[leg[2]], &leptons[leg[3]]};
  P4 d[4];
  for (int i = 0; i < 4; ++i) d[i] = dressed(*l[i], pt[i]);
  k.mz1 = (d[0] + d[1]).M();
  k.mz2 = (d[2] + d[3]).M();
  k.total = d[0] + d[1] + d[2] + d[3];
  k.m4l = k.total.M();
  k.z1_closer = ss || std::fabs(k.mz1 - kZ) < std::fabs(k.mz2 - kZ);
  bool pass = k.mz1 > cfg.z1_low && k.mz1 < cfg.z1_high && k.mz2 > z2_low && k.mz2 < cfg.z2_high && k.m4l > cfg.min_m4l;
  std::array<double, 4> pts = pt;
  std::sort(pts.begin(), pts.end(), std::greater<double>());
  pass = pass && pts[0] > cfg.lead_pt && pts[1] > cfg.sublead_pt;
  for (int i = 0; i < 4 && pass; ++i)
    for (int j = i + 1; j < 4 && pass; ++j) {
      if (h4l::delta_r(l[i]->eta, l[i]->phi, l[j]->eta, l[j]->phi) <= cfg.min_dr) pass = false;
      if (l[i]->charge() != l[j]->charge() && (lepton_p4(*l[i], pt[i]) + lepton_p4(*l[j], pt[j])).M() <= cfg.min_os_mass)
        pass = false;
    }
  // Smart cut (4e, 4mu with an opposite-sign Z2): the alternative pairing.
  if (!ss && std::abs(l[0]->pdg) == std::abs(l[2]->pdg)) {
    const int partner_a = l[2]->charge() != l[0]->charge() ? 2 : 3;
    const int partner_b = partner_a == 2 ? 3 : 2;
    const double m1 = (d[0] + d[partner_a]).M(), m2 = (d[1] + d[partner_b]).M();
    const bool first = std::fabs(m1 - kZ) < std::fabs(m2 - kZ);
    k.mza = first ? m1 : m2;
    k.mzb = first ? m2 : m1;
    if (std::fabs(k.mza - kZ) < std::fabs(k.mz1 - kZ) && k.mzb < z2_low) pass = false;
  }
  k.pass_cuts = pass;
  return k;
}

// Per-event m4l uncertainty: quadrature of the m4l changes for each lepton (and photon) moved by its
// error; cov01 is the covariance of the pT of legs 0 and 1 (non-zero after the Z1 refit).
double mass_error(const std::vector<Lepton>& leptons, const std::array<int, 4>& leg, const std::array<double, 4>& pt,
                  const std::array<double, 4>& sigma, bool with_photons, double cov01 = 0.0) {
  P4 d[4];
  for (int i = 0; i < 4; ++i) d[i] = dressed(leptons[leg[i]], pt[i]);
  const double m0 = (d[0] + d[1] + d[2] + d[3]).M();
  double sum = 0, jacobian[2] = {0, 0};
  for (int i = 0; i < 4; ++i) {
    P4 moved[4] = {d[0], d[1], d[2], d[3]};
    moved[i] = dressed(leptons[leg[i]], pt[i] + sigma[i]);
    const double dm = (moved[0] + moved[1] + moved[2] + moved[3]).M() - m0;
    sum += dm * dm;
    if (i < 2 && sigma[i] > 0) jacobian[i] = dm / sigma[i];
    const Lepton& l = leptons[leg[i]];
    if (with_photons && l.fsr_pt > 0) {
      moved[i] = lepton_p4(l, pt[i]) + P4(l.fsr_pt + l.fsr_sigma, l.fsr_eta, l.fsr_phi, 0.0);
      const double dg = (moved[0] + moved[1] + moved[2] + moved[3]).M() - m0;
      sum += dg * dg;
    }
  }
  sum += 2.0 * jacobian[0] * jacobian[1] * cov01;
  return std::sqrt(std::max(sum, 0.0));
}

// Z1 kinematic refit (AN 5.4): maximize Gauss(pT1) Gauss(pT2) f(m12), photons fixed.
struct Refit {
  bool ok = false;
  int status = -1;  // 10 x MIGRAD status + covariance status
  double pt1 = 0, pt2 = 0, err1 = 0, err2 = 0, cov12 = 0;
};
Refit refit_z1(const Lepton& a, const Lepton& b, const Lineshape& shape) {
  Refit r;
  const double s1 = a.sigma, s2 = b.sigma;
  if (!(s1 > 0) || !(s2 > 0)) return r;
  auto nll = [&](const double* x) {
    const double p1 = a.pt * (1 + x[0]), p2 = b.pt * (1 + x[1]);
    if (p1 <= 0 || p2 <= 0) return 1e30;
    const double m = (dressed(a, p1) + dressed(b, p2)).M();
    const double g1 = (p1 - a.pt) / s1, g2 = (p2 - b.pt) / s2;
    return 0.5 * (g1 * g1 + g2 * g2) - shape.log_at(m);
  };
  std::unique_ptr<ROOT::Math::Minimizer> minimizer(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
  ROOT::Math::Functor functor(nll, 2);
  minimizer->SetFunction(functor);
  minimizer->SetErrorDef(0.5);
  minimizer->SetStrategy(1);
  minimizer->SetPrintLevel(-1);
  minimizer->SetTolerance(0.01);  // EDM goal 1e-5: pT shifts to about 0.01 sigma
  minimizer->SetMaxFunctionCalls(20000);
  const double w1 = s1 / a.pt, w2 = s2 / b.pt;
  minimizer->SetLimitedVariable(0, "x1", 0.0, 0.1 * w1, -std::min(0.9, 8 * w1), 8 * w1);
  minimizer->SetLimitedVariable(1, "x2", 0.0, 0.1 * w2, -std::min(0.9, 8 * w2), 8 * w2);
  minimizer->Minimize();
  const int migrad_status = minimizer->Status();
  const double* xm = minimizer->X();
  const double x[2] = {xm[0], xm[1]};
  // Covariance: the inverse of the Hessian of the NLL by central differences with steps of 0.3 sigma
  // of each lepton (the line shape varies on the scale of Gamma_Z, far above the interpolation grid).
  const double h[2] = {0.3 * w1, 0.3 * w2};
  auto f = [&](double d0, double d1) {
    const double y[2] = {x[0] + d0, x[1] + d1};
    return nll(y);
  };
  const double f0 = f(0, 0);
  const double h00 = (f(h[0], 0) - 2 * f0 + f(-h[0], 0)) / (h[0] * h[0]);
  const double h11 = (f(0, h[1]) - 2 * f0 + f(0, -h[1])) / (h[1] * h[1]);
  const double h01 = (f(h[0], h[1]) - f(h[0], -h[1]) - f(-h[0], h[1]) + f(-h[0], -h[1])) / (4 * h[0] * h[1]);
  const double det = h00 * h11 - h01 * h01;
  const bool positive = h00 > 0 && h11 > 0 && det > 0 && std::isfinite(det);
  r.status = 10 * migrad_status + (positive ? 3 : 0);
  r.ok = (migrad_status == 0 || migrad_status == 1) && positive && std::isfinite(x[0]) && std::isfinite(x[1]);
  r.pt1 = a.pt * (1 + x[0]);
  r.pt2 = b.pt * (1 + x[1]);
  if (positive) {
    r.err1 = a.pt * std::sqrt(h11 / det);
    r.err2 = b.pt * std::sqrt(h00 / det);
    r.cov12 = -a.pt * b.pt * h01 / det;
  }
  return r;
}

int flavour_code(const Lepton& l) { return std::abs(l.pdg); }
int final_state_of(int z1, int z2) {
  if (z1 == 13 && z2 == 13) return 0;  // 4mu
  if (z1 == 11 && z2 == 11) return 1;  // 4e
  return 2;                            // 2e2mu (either Z1 flavour)
}
}  // namespace

int main(int argc, char** argv) {
  std::string task_path, out_json, out_root, selftest_path;
  double selftest_sqrts = 13.0, selftest_mh = 125.0;
  for (int index = 1; index < argc; index += 2) {
    const std::string key = argv[index];
    if (index + 1 >= argc) {
      std::cerr << "missing value for " << key << "\n";
      return 64;
    }
    if (key == "--task") task_path = argv[index + 1];
    else if (key == "--out-json") out_json = argv[index + 1];
    else if (key == "--out-root") out_root = argv[index + 1];
    else if (key == "--mela-selftest") selftest_path = argv[index + 1];
    else if (key == "--sqrts") selftest_sqrts = std::stod(argv[index + 1]);
    else if (key == "--mh") selftest_mh = std::stod(argv[index + 1]);
    else {
      std::cerr << "unknown argument " << key << "\n";
      return 64;
    }
  }
  gROOT->SetBatch(true);
  SetErrorHandler(record_errors);
  if (!selftest_path.empty()) {
    try {
      if (!h4lmela::init(selftest_sqrts, selftest_mh, 1)) throw std::runtime_error("MELA initialization failed");
      if (g_read_error) throw std::runtime_error("ROOT error during the MELA initialization: " + g_read_message);
      const json values = mela_selftest();
      h4l::write_json_file(selftest_path, {{"sqrts_tev", selftest_sqrts}, {"mh", selftest_mh}, {"cases", values}});
    } catch (const std::exception& error) {
      std::cerr << "ERROR: " << error.what() << "\n";
      return 1;
    }
    return 0;
  }
  if (task_path.empty() || out_json.empty() || out_root.empty()) {
    std::cerr << "usage: h4l_select --task TASK.json --out-json OUT.json --out-root OUT.root\n"
                 "       h4l_select --mela-selftest OUT.json --sqrts TEV --mh GEV\n";
    return 64;
  }
  try {
    const json task = h4l::read_json(task_path);
    const Config cfg = read_config(task.at("selection_config"));
    const bool is_mc = task.at("config").at("kind").get<std::string>() == "mc";
    const std::string role = is_mc ? "mc" : "data";
    const auto muon_payload = h4l::FactorizedPayload::from_json(task.at("calibration_payload").at("muon"));
    const auto electron_payload = h4l::FactorizedPayload::from_json(task.at("calibration_payload").at("electron"));
    // The event records store the pT calibrated with their own payload; this task recalibrates from the raw pT.
    const bool same_payload_as_reco = task.at("calibration").at("same_as_event_records").get<bool>();
    const auto lambda_mu = read_lambda(task.at("lambda_payload").at("flavours").at("muon"), role);
    const auto lambda_el = read_lambda(task.at("lambda_payload").at("flavours").at("electron"), role);
    Lineshape shape;
    {
      const json& z1 = task.at("refit_inputs").at("z1_lineshape");
      const auto& edges = z1.at("edges");
      if (edges.size() != z1.at("density").size() + 1 || edges.size() < 3)
        throw std::runtime_error("the Z1 line shape needs one edge more than densities");
      shape.lo = edges.at(0).get<double>();
      shape.width = edges.at(1).get<double>() - shape.lo;
      for (std::size_t i = 1; i < edges.size(); ++i)
        if (std::fabs(edges.at(i).get<double>() - edges.at(i - 1).get<double>() - shape.width) > 1e-9 * shape.width)
          throw std::runtime_error("the Z1 line shape needs uniform bins");
      for (const auto& v : z1.at("density")) shape.log_density.push_back(std::log(std::max(v.get<double>(), 1e-12)));
    }
    const double fsr_a = task.at("refit_inputs").at("fsr_photon_resolution").at("a").get<double>();
    const double fsr_b = task.at("refit_inputs").at("fsr_photon_resolution").at("b").get<double>();
    // MELA: verify the libraries, dictionary and data files against the plan; copy the files MELA opens for
    // update into the (empty) working directory; initialize; repeat the self-test.
    const json& mela_cfg = task.at("mela");
    json mela_files = json::array();
    std::map<std::string, std::string> mela_digest;
    if (mela_cfg.at("files").empty()) throw std::runtime_error("the task lists no MELA files");
    for (const auto& item : mela_cfg.at("files")) {
      const std::string path = item.at("path").get<std::string>();
      const std::string digest = h4l::sha256_file(path);
      if (digest != item.at("sha256").get<std::string>()) throw std::runtime_error("MELA file changed: " + path);
      mela_files.push_back({{"path", path}, {"sha256", digest}});
      mela_digest[path] = digest;
    }
    json working_copies = json::array();
    for (const auto& item : mela_cfg.at("working_copies")) {
      const std::string source = item.at("source").get<std::string>(), name = item.at("name").get<std::string>();
      if (!mela_digest.count(source)) throw std::runtime_error("working copy of an unlisted MELA file: " + source);
      const std::filesystem::path local(name);
      if (local.is_absolute() || name.find("..") != std::string::npos) throw std::runtime_error("bad working-copy name " + name);
      // Never write through an existing name (a symbolic link would lead back to the installation).
      std::error_code status_error;
      if (std::filesystem::exists(std::filesystem::symlink_status(local, status_error)))
        throw std::runtime_error("the working directory already holds " + name + " (a fresh directory is required)");
      if (local.has_parent_path()) std::filesystem::create_directories(local.parent_path());
      std::filesystem::copy_file(source, local, std::filesystem::copy_options::none);
      if (h4l::sha256_file(local.string()) != mela_digest.at(source)) throw std::runtime_error("bad working copy of " + source);
      working_copies.push_back({{"name", name}, {"source", source}});
    }
    g_read_error = false;
    if (!h4lmela::init(mela_cfg.at("sqrts_tev").get<double>(), mela_cfg.at("mh").get<double>(), 1))
      throw std::runtime_error("MELA initialization failed");
    // ROOT errors during the initialization (e.g. an XRootD open that fell back to the FUSE mount) are
    // recorded; the self-test decides whether MELA works.
    const json mela_init_errors = g_read_error ? json(g_read_message) : json();
    g_read_error = false;
    const json selftest = mela_selftest();
    const std::string mismatch = selftest_mismatch(selftest, mela_cfg.at("selftest").at("cases"));
    if (!mismatch.empty()) throw std::runtime_error("MELA self-test failed: " + mismatch);
    std::set<ULong64_t> planned_keys;
    for (const auto& key : task.at("source_file_keys")) planned_keys.insert(h4l::parse_file_key(key.get<std::string>()));

    std::unique_ptr<TFile> output(TFile::Open(out_root.c_str(), "RECREATE"));
    if (!output || output->IsZombie()) throw std::runtime_error("cannot create " + out_root);
    output->cd();
    std::vector<std::unique_ptr<CandidateRow>> region_rows;
    std::vector<TTree*> region_trees;
    for (const auto& [name, z2_low] : cfg.regions) {
      (void)z2_low;
      region_rows.push_back(std::make_unique<CandidateRow>());
      region_trees.push_back(new TTree(name.c_str(), ("selected candidates, " + name).c_str()));
      region_rows.back()->book(*region_trees.back());
    }
    CandidateRow cr_row;
    TTree& cr_tree = *new TTree("CR", "Z+X control regions: cr_type 0 = 2P2F, 1 = 3P1F, 2 = SS");
    cr_row.book(cr_tree);
    ZlRow zl_row;
    TTree& zl_tree = *new TTree("ZL", "Z + 1 loose lepton (fake rates)");
    zl_row.book(zl_tree);
    TTree& inputs_tree = *new TTree("Inputs", "per-input counters");
    std::string input_path;
    Long64_t input_rows = 0;
    inputs_tree.Branch("path", &input_path);
    inputs_tree.Branch("rows", &input_rows, "rows/L");

    std::map<std::string, long long> counts;
    for (const char* key : {"rows", "trigger_pv", "refit_ok", "refit_failed", "mela_failed", "mela_decay_failed", "zl",
                            "lepton_truncated", "jet_truncated", "candidate_truncated", "lepton_unmeasured"})
      counts[key] = 0;
    for (const auto& [name, z2_low] : cfg.regions) {
      (void)z2_low;
      const std::string tag = name == "SR" ? "" : "_" + name;
      counts["selected_" + name] = 0;
      for (const char* cr : {"cr_2p2f", "cr_3p1f", "cr_ss"}) counts[cr + tag] = 0;
    }
    json inputs_report = json::array(), open_errors = json::array();
    for (const auto& input : task.at("inputs")) {
      input_path = input.at("root").get<std::string>();
      const long long expected = input.at("rows").get<long long>();
      std::string open_error;
      auto file = h4l::open_input(input_path, 5, open_error);
      if (!file) throw InputError("cannot open " + input_path + ": " + open_error);
      // Errors of failed open attempts (retried) are recorded, then cleared.
      if (g_read_error) open_errors.push_back({{"input", input_path}, {"message", g_read_message}});
      g_read_error = false;
      auto* tree = dynamic_cast<TTree*>(file->Get("Events4l"));
      if (!tree) throw InputError("no Events4l in " + input_path);
      if (tree->GetEntries() != expected) throw InputError(input_path + ": Events4l entries differ from the plan");
      input_rows = tree->GetEntries();
      {
      Reader in(tree);
      for (Long64_t e = 0; e < input_rows; ++e) {
        if (in.reader.SetEntry(e) != TTreeReader::kEntryValid)
          throw InputError("read error at entry " + std::to_string(e) + " of " + input_path);
        ++counts["rows"];
        if (!planned_keys.count(*in.file_key)) throw InputError("event of a file key outside the plan in " + input_path);
        if (*in.nlep_all > *in.nlep) ++counts["lepton_truncated"];
        if (*in.njet_all > *in.njet) ++counts["jet_truncated"];
        if (*in.ncand_overflow > 0) ++counts["candidate_truncated"];
        if (!*in.trigger || *in.npv < cfg.min_good_pv) continue;
        ++counts["trigger_pv"];
        // Leptons.
        const int n = *in.nlep;
        std::vector<Lepton> leptons(n);
        for (int k = 0; k < n; ++k) {
          Lepton& l = leptons[k];
          l.k = k;
          l.pdg = in.l_pdg[k];
          l.flags = in.l_flags[k];
          l.pt = in.l_pt[k];
          l.pt_raw = in.l_pt_raw[k];
          l.eta = in.l_eta[k];
          l.eta_sc = in.l_eta_sc[k];
          l.phi = in.l_phi[k];
          l.normal = in.l_normal[k];
          l.rel_err = in.l_rel_err[k];
          l.iso_all = in.l_iso[k];
          l.iso_chg = in.l_iso_chg[k];
          l.overlap_mask = in.l_overlap_mask[k];
          l.sip = in.l_sip[k];
          l.gen_flav = in.l_gen_flav[k];
          l.gen_pt = in.l_gen_pt[k];
          l.lost_hits = in.l_lost_hits[k];
          l.fsr_pt = in.l_fsr_pt[k];
          l.fsr_eta = in.l_fsr_eta[k];
          l.fsr_phi = in.l_fsr_phi[k];
          if (l.fsr_pt > 0) l.fsr_sigma = l.fsr_pt * std::sqrt(fsr_a * fsr_a / l.fsr_pt + fsr_b * fsr_b);
          // A lepton without finite kinematics, a finite positive momentum error, isolation and SIP (one
          // muon with a NaN error in 564k leptons of TTBar_0013) has no usable measurement: it is not loose.
          const bool measured = std::isfinite(l.pt) && std::isfinite(l.pt_raw) && std::isfinite(l.eta) && std::isfinite(l.eta_sc) &&
                                std::isfinite(l.phi) && std::isfinite(l.rel_err) && l.rel_err > 0 && std::isfinite(l.iso_all) &&
                                std::isfinite(l.iso_chg) && std::isfinite(l.sip) &&
                                (!(l.fsr_pt > 0) || (std::isfinite(l.fsr_eta) && std::isfinite(l.fsr_phi)));
          if (!measured) ++counts["lepton_unmeasured"];
          const double abs_eta = std::fabs(l.muon() ? l.eta : l.eta_sc);
          // The calibration of this task, applied to the raw pT (data pT exp(-u); MC pT (1 + r N) with the stored
          // per-lepton deviate); with the payload of the event records it must reproduce the stored pT.
          l.pt_reco = l.pt;
          if (measured) {
            const auto& payload = l.muon() ? muon_payload : electron_payload;
            l.pt = is_mc ? payload.smeared_mc_pt(l.pt_raw, abs_eta, l.normal) : payload.corrected_data_pt(l.pt_raw, abs_eta);
            if (same_payload_as_reco && std::fabs(l.pt - l.pt_reco) > 1e-4 * l.pt_reco)
              throw std::runtime_error("the calibration of the event records is not reproduced in " + input_path);
          }
          l.sigma = measured ? lambda_of(l.muon() ? lambda_mu : lambda_el, abs_eta, l.rel_err) * l.rel_err * l.pt : 0.0;
          l.loose = measured && (l.muon() ? (l.pt > cfg.mu_pt && std::fabs(l.eta) < cfg.mu_eta)
                                          : (l.pt > cfg.el_pt && std::fabs(l.eta) < cfg.el_eta));
          l.tight_id = l.muon() ? (!cfg.muon_id || h4l::an_tight_muon(l.flags, l.pt, cfg.muon_high_pt))
                                : (cfg.electron_id_bit < 0 || (l.flags >> cfg.electron_id_bit & 1U));
        }
        // Ghost cleaning beyond dR < 0.02 (substitute for the AN segment arbitration): a tracker-only muon
        // (neither global nor PF) within dR < 0.05 of a same-charge PF muon is a duplicate.
        for (auto& a : leptons) {
          if (!a.muon() || (a.flags & 1U) || (a.flags >> 3 & 1U)) continue;
          for (const auto& b : leptons)
            if (&a != &b && b.muon() && (b.flags >> 3 & 1U) && b.pdg == a.pdg &&
                h4l::delta_r(a.eta, a.phi, b.eta, b.phi) < cfg.ghost_dr)
              a.ghost = true;
          if (a.ghost) a.loose = false;
        }
        // AN FSR-subtracted isolation (3.3): the photons of loose muons passing SIP are removed from the
        // neutral isolation of every loose lepton passing SIP (flavour-specific veto); other leptons keep
        // their plain isolation.
        for (auto& l : leptons) {
          l.iso = l.iso_all;
          if (!l.loose || l.sip >= cfg.max_sip) continue;
          double photons = 0;
          for (const auto& o : leptons) {
            if (!o.muon() || !o.loose || o.sip >= cfg.max_sip || !(o.fsr_pt > 0)) continue;
            const double dr = h4l::delta_r(l.eta, l.phi, o.fsr_eta, o.fsr_phi);
            if (h4l::fsr_in_isolation(l.muon(), std::fabs(l.eta_sc), dr)) photons += o.fsr_pt;
          }
          l.iso = l.iso_chg + std::max(0.0, (l.iso_all - l.iso_chg) - photons / l.pt_raw);
        }
        for (auto& l : leptons) l.selected = l.loose && l.tight_id && l.sip < cfg.max_sip && l.iso < cfg.max_iso;
        // Cross cleaning: an electron within dR < 0.05 of any selected muon is removed.
        for (auto& l : leptons) {
          if (l.muon()) continue;
          for (int m = 0; m < n; ++m)
            if ((l.overlap_mask >> m & 1U) && leptons[m].muon() && leptons[m].selected) {
              l.cleaned = true;
              l.loose = l.selected = false;
            }
        }
        // Candidates.
        std::vector<Candidate> candidates;
        for (int c = 0; c < *in.ncand; ++c)
          candidates.push_back({c, {in.c_l1[c], in.c_l2[c], in.c_l3[c], in.c_l4[c]}, static_cast<bool>(in.c_ss[c])});
        auto nominal_pt = [&](const Candidate& c) {
          return std::array<double, 4>{leptons[c.leg[0]].pt, leptons[c.leg[1]].pt, leptons[c.leg[2]].pt, leptons[c.leg[3]].pt};
        };
        // Cross-check of the stored masses (the pT of the event records, the reco floors).
        for (const auto& c : candidates) {
          const std::array<double, 4> reco_pt = {leptons[c.leg[0]].pt_reco, leptons[c.leg[1]].pt_reco, leptons[c.leg[2]].pt_reco,
                                                 leptons[c.leg[3]].pt_reco};
          const Kinematics k = kinematics(leptons, c.leg, reco_pt, c.ss, cfg, cfg.regions.front().second);
          if (std::fabs(k.m4l - in.c_m4l[c.index]) > 2e-3 * k.m4l || std::fabs(k.mz1 - in.c_mz1[c.index]) > 2e-3 * k.mz1 ||
              std::fabs(k.mz2 - in.c_mz2[c.index]) > 2e-3 * k.mz2)
            throw std::runtime_error("recomputed candidate masses differ from h4l_reco in " + input_path);
        }
        // MELA decay-only discriminant per candidate (cached).
        std::map<int, double> dbkg_cache;
        // MELA input; a same-sign Z2 (SS region) gets its second lepton's charge flipped: the production
        // matrix elements treat the 4l system as the Higgs and do not depend on the charges, D_bkg^kin is then
        // only an approximation.
        auto mela_input = [&](const Candidate& c, const std::vector<std::pair<P4, int>>& jets) {
          h4lmela::MelaInput mi;
          for (int i = 0; i < 4; ++i) {
            const Lepton& l = leptons[c.leg[i]];
            const P4 p = dressed(l, l.pt);
            const int pdg = (c.ss && i == 3) ? -l.pdg : l.pdg;
            mi.leptons[i] = {p.Pt(), p.Eta(), p.Phi(), std::max(p.M(), 0.0), pdg};
          }
          for (const auto& [p, index] : jets) {
            (void)index;
            mi.jets.push_back({p.Pt(), p.Eta(), p.Phi(), p.M()});
          }
          return mi;
        };
        auto dbkg = [&](const Candidate& c) {
          const auto found = dbkg_cache.find(c.index);
          if (found != dbkg_cache.end()) return found->second;
          h4lmela::MelaOutput mo;
          const bool ok = h4lmela::compute(mela_input(c, {}), mo, false);
          if (!ok) ++counts["mela_decay_failed"];
          const double value = ok ? mo.D_bkg_kin : -1.0;
          dbkg_cache[c.index] = value;
          return value;
        };
        // Cleaned jets (AN 3.4): away from the selected leptons, from the leptons of the candidate (the
        // failing leptons of a control-region candidate pass in its signal-region analog) and from the
        // selected FSR photons (those of loose muons passing SIP).
        auto clean_jets = [&](const std::array<int, 4>& legs) {
          std::vector<std::pair<P4, int>> out;
          for (int j = 0; j < *in.njet; ++j) {
            if (in.j_pt[j] <= cfg.jet_pt || std::fabs(in.j_eta[j]) >= cfg.jet_eta) continue;
            bool clean = true;
            for (const auto& l : leptons) {
              const bool in_candidate = l.k == legs[0] || l.k == legs[1] || l.k == legs[2] || l.k == legs[3];
              if ((l.selected || in_candidate) && h4l::delta_r(l.eta, l.phi, in.j_eta[j], in.j_phi[j]) < cfg.jet_clean_dr)
                clean = false;
              const bool photon_selected = l.muon() && l.loose && l.sip < cfg.max_sip;
              if ((photon_selected || in_candidate) && l.fsr_pt > 0 &&
                  h4l::delta_r(l.fsr_eta, l.fsr_phi, in.j_eta[j], in.j_phi[j]) < cfg.jet_clean_dr)
                clean = false;
            }
            if (clean) out.push_back({P4(in.j_pt[j], in.j_eta[j], in.j_phi[j], in.j_mass[j]), j});
          }
          std::sort(out.begin(), out.end(), [](const auto& a, const auto& b) { return a.first.Pt() > b.first.Pt(); });
          return out;
        };

        // Fill the observables of one chosen candidate.  Every field starts from its default (no value of
        // an earlier row survives); the observables do not depend on the region, so a candidate chosen in
        // several regions is computed once per event (cache) and only region, cr_type and ncand_region differ.
        struct Filled {
          CandidateRow row;
          bool refit_ok = false, mela_ok = false;
        };
        std::map<int, Filled> filled_cache;
        auto fill = [&](CandidateRow& row, const Candidate& c, int cr_type, int ncand_region, double z2_low, int region) {
          const auto cached = filled_cache.find(c.index);
          if (cached != filled_cache.end()) {
            row = cached->second.row;
            row.region = region;
            row.cr_type = cr_type;
            row.ncand_region = ncand_region;
            ++counts[cached->second.refit_ok ? "refit_ok" : "refit_failed"];
            if (!cached->second.mela_ok) ++counts["mela_failed"];
            return;
          }
          row = CandidateRow{};
          row.region = region;
          bool refit_ok = false, mela_ok = false;
          const auto pt = nominal_pt(c);
          const Kinematics k = kinematics(leptons, c.leg, pt, c.ss, cfg, z2_low);
          row.file_key = *in.file_key;
          row.entry = *in.entry;
          row.weight = *in.weight;
          row.met = *in.met;
          row.met_phi = *in.met_phi;
          row.npv = *in.npv;
          row.trig_bits = *in.trig_bits;
          row.cr_type = cr_type;
          row.ncand_region = ncand_region;
          row.z1_flavour = flavour_code(leptons[c.leg[0]]);
          row.z2_flavour = flavour_code(leptons[c.leg[2]]);
          row.final_state = final_state_of(row.z1_flavour, row.z2_flavour);
          row.m4l = k.m4l;
          row.mz1 = k.mz1;
          row.mz2 = k.mz2;
          row.pt4l = k.total.Pt();
          row.eta4l = k.total.Eta();
          row.phi4l = k.total.Phi();
          row.y4l = k.total.Rapidity();
          std::array<double, 4> sigma{};
          for (int i = 0; i < 4; ++i) {
            const Lepton& l = leptons[c.leg[i]];
            sigma[i] = l.sigma;
            row.l_pdg[i] = l.pdg;
            row.l_gen_flav[i] = l.gen_flav;
            row.l_lost_hits[i] = l.lost_hits;
            row.l_pass[i] = (l.tight_id && l.iso < cfg.max_iso) ? 1 : 0;
            row.l_index[i] = l.k;
            row.l_flags[i] = static_cast<Int_t>(l.flags);
            row.l_pt[i] = l.pt;
            row.l_pt_raw[i] = l.pt_raw;
            row.l_eta[i] = l.eta;
            row.l_eta_sc[i] = l.eta_sc;
            row.l_phi[i] = l.phi;
            row.l_sigma[i] = l.sigma;
            row.l_iso[i] = l.iso;
            row.l_sip[i] = l.sip;
            row.l_mva[i] = in.l_mva_noiso[l.k];
            row.l_gen_pt[i] = l.gen_pt;
            row.l_trig_bits[i] = in.l_trig_bits[l.k];
            row.l_trig_pt[i] = in.l_trig_pt[l.k];
            row.l_fsr_pt[i] = l.fsr_pt;
            row.l_fsr_eta[i] = l.fsr_eta;
            row.l_fsr_phi[i] = l.fsr_phi;
          }
          row.m4l_err = mass_error(leptons, c.leg, pt, sigma, true);
          // Z1 refit.
          const Refit r = refit_z1(leptons[c.leg[0]], leptons[c.leg[1]], shape);
          row.refit_status = r.status;
          if (r.ok) {
            std::array<double, 4> refit_pt = pt;
            refit_pt[0] = r.pt1;
            refit_pt[1] = r.pt2;
            const Kinematics kr = kinematics(leptons, c.leg, refit_pt, c.ss, cfg, z2_low);
            row.m4l_refit = kr.m4l;
            row.mz1_refit = kr.mz1;
            std::array<double, 4> refit_sigma = sigma;
            refit_sigma[0] = r.err1;
            refit_sigma[1] = r.err2;
            row.m4l_refit_err = mass_error(leptons, c.leg, refit_pt, refit_sigma, true, r.cov12);
            ++counts["refit_ok"];
            refit_ok = true;
          } else {
            row.m4l_refit = k.m4l;
            row.mz1_refit = k.mz1;
            row.m4l_refit_err = row.m4l_err;
            ++counts["refit_failed"];
          }
          // Additional selected leptons.
          std::vector<int> extra;
          for (const auto& l : leptons)
            if (l.selected && l.k != c.leg[0] && l.k != c.leg[1] && l.k != c.leg[2] && l.k != c.leg[3]) extra.push_back(l.k);
          row.n_extra_leptons = static_cast<int>(extra.size());
          row.n_extra_os_sf_pairs = 0;
          for (std::size_t a = 0; a < extra.size(); ++a)
            for (std::size_t b = a + 1; b < extra.size(); ++b)
              if (leptons[extra[a]].pdg == -leptons[extra[b]].pdg) ++row.n_extra_os_sf_pairs;
          // Jets.
          const auto jets = clean_jets(c.leg);
          row.njets = static_cast<int>(jets.size());
          row.nbjets = 0;
          for (int i = 0; i < 4; ++i) {
            row.jet_pt[i] = row.jet_eta[i] = row.jet_phi[i] = row.jet_mass[i] = row.jet_btag[i] = -99;
            row.jet_flavour[i] = -1;
          }
          for (std::size_t j = 0; j < jets.size(); ++j) {
            const int index = jets[j].second;
            if (std::fabs(in.j_eta[index]) < cfg.btag_eta && in.j_deepjet[index] > cfg.btag_threshold) ++row.nbjets;
            if (j < 4) {
              row.jet_pt[j] = jets[j].first.Pt();
              row.jet_eta[j] = jets[j].first.Eta();
              row.jet_phi[j] = jets[j].first.Phi();
              row.jet_mass[j] = jets[j].first.M();
              row.jet_btag[j] = in.j_deepjet[index];
              row.jet_flavour[j] = in.j_flavour[index];
            }
          }
          row.mjj = jets.size() >= 2 ? (jets[0].first + jets[1].first).M() : -1;
          row.detajj = jets.size() >= 2 ? std::fabs(jets[0].first.Eta() - jets[1].first.Eta()) : -1;
          // MELA with production (the fields start from their defaults).
          row.mela_charge_flipped = c.ss ? 1 : 0;
          {
            h4lmela::MelaOutput mo;
            if (h4lmela::compute(mela_input(c, jets), mo, true)) {
              row.p_sig_ggH = mo.p_sig_ggH;
              row.p_bkg_qqZZ = mo.p_bkg_qqZZ;
              row.c_bkg_kin = mo.c_bkg_kin;
              row.d_bkg_kin = mo.D_bkg_kin;
              row.p_JJVBF = mo.p_JJVBF_sig;
              row.p_JJQCD = mo.p_JJQCD_sig;
              row.p_HadWH = mo.p_HadWH_sig;
              row.p_HadZH = mo.p_HadZH_sig;
              row.p_HadWH_mavjj = mo.p_HadWH_mavjj;
              row.p_HadWH_mavjj_true = mo.p_HadWH_mavjj_true;
              row.p_HadZH_mavjj = mo.p_HadZH_mavjj;
              row.p_HadZH_mavjj_true = mo.p_HadZH_mavjj_true;
              row.p_JQCD = mo.p_JQCD_sig;
              row.p_JVBF = mo.p_JVBF_sig;
              row.pAux_JVBF = mo.pAux_JVBF_sig;
              row.c_2jet = mo.c_2jet;
              row.c_1jet = mo.c_1jet;
              row.c_WH = mo.c_WH;
              row.c_ZH = mo.c_ZH;
              row.d_2jet = mo.D_2jet;
              row.d_1jet = mo.D_1jet;
              row.d_wh_mela = mo.D_WH;
              row.d_zh_mela = mo.D_ZH;
              row.mela_njets = mo.nJets;
              mela_ok = true;
            } else {
              ++counts["mela_failed"];
            }
          }
          // MC: responses of m4l to the lepton scale and smearing per flavour.
          row.dm4l_scale_mu = row.dm4l_scale_e = row.dm4l_smear_mu = row.dm4l_smear_e = 0;
          if (is_mc) {
            for (int flavour : {13, 11}) {
              const double epsilon = 1e-3, delta_r = 1e-3;
              std::array<double, 4> scaled = pt, smeared = pt;
              for (int i = 0; i < 4; ++i) {
                const Lepton& l = leptons[c.leg[i]];
                if (std::abs(l.pdg) != flavour) continue;
                scaled[i] = pt[i] * (1 + epsilon);
                const double abs_eta = std::fabs(l.muon() ? l.eta : l.eta_sc);
                const double r = (l.muon() ? muon_payload : electron_payload).smear_at(l.pt_raw, abs_eta);
                smeared[i] = l.pt_raw * (1 + (r + delta_r) * l.normal);
              }
              const double ds = (kinematics(leptons, c.leg, scaled, c.ss, cfg, z2_low).m4l - k.m4l) / epsilon;
              const double dk = (kinematics(leptons, c.leg, smeared, c.ss, cfg, z2_low).m4l - k.m4l) / delta_r;
              (flavour == 13 ? row.dm4l_scale_mu : row.dm4l_scale_e) = static_cast<Float_t>(ds);
              (flavour == 13 ? row.dm4l_smear_mu : row.dm4l_smear_e) = static_cast<Float_t>(dk);
            }
          }
          row.n_scale = std::min<int>(*in.n_scale, 16);
          row.n_pdf = std::min<int>(*in.n_pdf, 128);
          row.n_ps = std::min<int>(*in.n_ps, 8);
          for (int i = 0; i < row.n_scale; ++i) row.w_scale[i] = in.w_scale[i];
          for (int i = 0; i < row.n_pdf; ++i) row.w_pdf[i] = in.w_pdf[i];
          for (int i = 0; i < row.n_ps; ++i) row.w_ps[i] = in.w_ps[i];
          row.g_present = *in.g_present;
          row.g_fid_pass = *in.g_fid_pass;
          row.g_stage0 = *in.g_stage0;
          row.g_vh_class = *in.g_vh_class;
          row.g_final_state_true = *in.g_final_state_true;
          row.g_fid_final_state = *in.g_fid_final_state;
          row.g_fid_njets = *in.g_fid_njets;
          row.g_h_pt = *in.g_h_pt;
          row.g_h_eta = *in.g_h_eta;
          row.g_h_phi = *in.g_h_phi;
          row.g_h_mass = *in.g_h_mass;
          row.g_h_y = *in.g_h_y;
          row.g_fid_m4l = *in.g_fid_m4l;
          row.g_fid_mz1 = *in.g_fid_mz1;
          row.g_fid_mz2 = *in.g_fid_mz2;
          row.g_fid_pt4l = *in.g_fid_pt4l;
          row.g_fid_y4l = *in.g_fid_y4l;
          row.g_fid_jet1_pt = *in.g_fid_jet1_pt;
          filled_cache[c.index] = {row, refit_ok, mela_ok};
        };

        // Best candidate among those passing: per set of four leptons the Z1 closest to m_Z, then max D_bkg^kin.
        auto choose = [&](const std::vector<std::pair<Candidate, Kinematics>>& pass) -> int {
          std::map<std::array<int, 4>, int> best_of_set;
          for (int i = 0; i < static_cast<int>(pass.size()); ++i) {
            std::array<int, 4> set = pass[i].first.leg;
            std::sort(set.begin(), set.end());
            const auto found = best_of_set.find(set);
            if (found == best_of_set.end() ||
                std::fabs(pass[i].second.mz1 - kZ) < std::fabs(pass[found->second].second.mz1 - kZ))
              best_of_set[set] = i;
          }
          int best = -1;
          double best_d = -2;
          for (const auto& [set, i] : best_of_set) {
            (void)set;
            const double d = dbkg(pass[i].first);
            if (d > best_d) {
              best_d = d;
              best = i;
            }
          }
          return best;
        };

        std::vector<bool> in_region(cfg.regions.size(), false);
        for (std::size_t r = 0; r < cfg.regions.size(); ++r) {
          const double z2_low = cfg.regions[r].second;
          std::vector<std::pair<Candidate, Kinematics>> pass;
          for (const auto& c : candidates) {
            if (c.ss) continue;
            bool all = true;
            for (int i : c.leg) all = all && leptons[i].selected;
            if (!all) continue;
            const Kinematics k = kinematics(leptons, c.leg, nominal_pt(c), c.ss, cfg, z2_low);
            if (!k.z1_closer || !k.pass_cuts) continue;
            pass.push_back({c, k});
          }
          if (pass.empty()) continue;
          const int best = choose(pass);
          fill(*region_rows[r], pass[best].first, -1, static_cast<int>(pass.size()), z2_low, static_cast<int>(r));
          output->cd();
          if (region_trees[r]->Fill() <= 0) throw std::runtime_error("region Fill failed");
          ++counts["selected_" + cfg.regions[r].first];
          in_region[r] = true;
        }
        // Z+X control regions of every signal region (events of that signal region excluded); the row's
        // region field is the index of the signal region whose selection it mirrors.
        for (std::size_t r = 0; r < cfg.regions.size(); ++r) {
          if (in_region[r]) continue;
          const double z2_low = cfg.regions[r].second;
          const std::string tag = r == 0 ? "" : "_" + cfg.regions[r].first;
          std::vector<std::pair<Candidate, Kinematics>> os[2], ss;
          for (const auto& c : candidates) {
            const Lepton &a = leptons[c.leg[0]], &b = leptons[c.leg[1]], &x = leptons[c.leg[2]], &y = leptons[c.leg[3]];
            if (!a.selected || !b.selected) continue;
            if (!x.loose || !y.loose || x.sip >= cfg.max_sip || y.sip >= cfg.max_sip) continue;
            const Kinematics k = kinematics(leptons, c.leg, nominal_pt(c), c.ss, cfg, z2_low);
            if (!k.pass_cuts) continue;
            if (c.ss) {
              if (k.m4l > cfg.ss_min_m4l) ss.push_back({c, k});  // the AN 7.2.2 cut m4l > 100 GeV is applied downstream
              continue;
            }
            const int fails = static_cast<int>(!(x.tight_id && x.iso < cfg.max_iso)) + static_cast<int>(!(y.tight_id && y.iso < cfg.max_iso));
            if (fails == 2) os[0].push_back({c, k});
            if (fails == 1) os[1].push_back({c, k});
          }
          for (int type = 0; type < 2; ++type) {
            if (os[type].empty()) continue;
            const int best = choose(os[type]);
            fill(cr_row, os[type][best].first, type, static_cast<int>(os[type].size()), z2_low, static_cast<int>(r));
            output->cd();
            if (cr_tree.Fill() <= 0) throw std::runtime_error("CR Fill failed");
            ++counts[(type == 0 ? "cr_2p2f" : "cr_3p1f") + tag];
          }
          if (!ss.empty()) {
            int best = 0;
            for (int i = 1; i < static_cast<int>(ss.size()); ++i) {
              const double di = std::fabs(ss[i].second.mz1 - kZ), db = std::fabs(ss[best].second.mz1 - kZ);
              const double si = leptons[ss[i].first.leg[2]].pt + leptons[ss[i].first.leg[3]].pt;
              const double sb = leptons[ss[best].first.leg[2]].pt + leptons[ss[best].first.leg[3]].pt;
              if (di < db - 1e-9 || (std::fabs(di - db) <= 1e-9 && si > sb)) best = i;
            }
            fill(cr_row, ss[best].first, 2, static_cast<int>(ss.size()), z2_low, static_cast<int>(r));
            output->cd();
            if (cr_tree.Fill() <= 0) throw std::runtime_error("CR Fill failed");
            ++counts["cr_ss" + tag];
          }
        }
        // Z + 1 loose lepton.
        {
          int za = -1, zb = -1;
          double best = 1e9;
          for (int a = 0; a < n; ++a)
            for (int b = a + 1; b < n; ++b) {
              const Lepton &x = leptons[a], &y = leptons[b];
              if (!x.selected || !y.selected || x.pdg != -y.pdg) continue;
              const double lead = std::max(x.pt, y.pt), sub = std::min(x.pt, y.pt);
              if (lead <= cfg.zl_lead || sub <= cfg.zl_sublead) continue;
              const double m = (dressed(x, x.pt) + dressed(y, y.pt)).M();
              if (m > cfg.zl_z1_low && m < cfg.zl_z1_high && std::fabs(m - kZ) < best) {
                best = std::fabs(m - kZ);
                za = a;
                zb = b;
              }
            }
          if (za >= 0) {
            std::vector<int> others;
            for (int k = 0; k < n; ++k)
              if (k != za && k != zb && leptons[k].loose) others.push_back(k);
            if (others.size() == 1) {
              const Lepton& p = leptons[others[0]];
              const Lepton& os_tag = leptons[za].charge() != p.charge() ? leptons[za] : leptons[zb];
              const double m_os = (lepton_p4(p, p.pt) + lepton_p4(os_tag, os_tag.pt)).M();
              if (m_os > cfg.zl_probe_os_mass) {
                zl_row.file_key = *in.file_key;
                zl_row.entry = *in.entry;
                zl_row.weight = *in.weight;
                zl_row.met = *in.met;
                zl_row.npv = *in.npv;
                zl_row.trig_bits = *in.trig_bits;
                zl_row.mz1 = (dressed(leptons[za], leptons[za].pt) + dressed(leptons[zb], leptons[zb].pt)).M();
                zl_row.z1_flavour = flavour_code(leptons[za]);
                zl_row.tag1_pt = std::max(leptons[za].pt, leptons[zb].pt);
                zl_row.tag2_pt = std::min(leptons[za].pt, leptons[zb].pt);
                zl_row.probe_pt = p.pt;
                zl_row.probe_eta = p.eta;
                zl_row.probe_eta_sc = p.eta_sc;
                zl_row.probe_phi = p.phi;
                zl_row.probe_iso = p.iso;
                zl_row.probe_sip = p.sip;
                zl_row.probe_mva = in.l_mva_noiso[p.k];
                zl_row.probe_os_mass = m_os;
                zl_row.probe_pdg = p.pdg;
                zl_row.probe_gen_flav = p.gen_flav;
                zl_row.probe_lost_hits = p.lost_hits;
                zl_row.probe_id = p.tight_id ? 1 : 0;
                zl_row.probe_flags = static_cast<Int_t>(p.flags);
                zl_row.probe_pass = p.selected ? 1 : 0;
                zl_row.n_loose = 0;
                zl_row.n_loose_sip = 0;
                for (const auto& l : leptons) {
                  zl_row.n_loose += l.loose ? 1 : 0;
                  zl_row.n_loose_sip += (l.loose && l.sip < cfg.max_sip) ? 1 : 0;
                }
                zl_row.m3l = (dressed(leptons[za], leptons[za].pt) + dressed(leptons[zb], leptons[zb].pt) + lepton_p4(p, p.pt)).M();
                zl_row.probe_dr_tag = std::min(h4l::delta_r(p.eta, p.phi, leptons[za].eta, leptons[za].phi),
                                               h4l::delta_r(p.eta, p.phi, leptons[zb].eta, leptons[zb].phi));
                output->cd();
                if (zl_tree.Fill() <= 0) throw std::runtime_error("ZL Fill failed");
                ++counts["zl"];
              }
            }
          }
        }
        check_read(input_path + " entry " + std::to_string(e));
      }
      }  // the reader is destroyed before its file
      check_read(input_path);
      if (inputs_tree.Fill() <= 0) throw std::runtime_error("Inputs Fill failed");
      inputs_report.push_back({{"root", input_path}, {"reco_task_id", input.at("reco_task_id")},
                               {"root_sha256", input.at("root_sha256")}, {"rows", input_rows}});
      std::cout << "[select] " << input_path << ": " << input_rows << " rows\n";
      file.reset();
      check_read(input_path + " (close)");
    }
    output->cd();
    std::map<std::string, long long> trees;
    for (std::size_t r = 0; r < region_trees.size(); ++r) {
      if (region_trees[r]->Write("", TObject::kOverwrite) <= 0) throw std::runtime_error("cannot write " + cfg.regions[r].first);
      trees[cfg.regions[r].first] = region_trees[r]->GetEntries();
    }
    for (TTree* t : {&cr_tree, &zl_tree, &inputs_tree}) {
      if (t->Write("", TObject::kOverwrite) <= 0) throw std::runtime_error(std::string("cannot write ") + t->GetName());
      trees[t->GetName()] = t->GetEntries();
    }
    output->Close();
    output.reset();
    if (g_read_error) throw std::runtime_error("ROOT error while writing: " + g_read_message);
    h4l::validate_root_output_full(out_root, trees);
    if (g_read_error) throw std::runtime_error("ROOT error while validating the output: " + g_read_message);
    char host[256] = {};
    if (gethostname(host, sizeof(host) - 1) != 0) std::snprintf(host, sizeof(host), "unknown");
    json report = {{"schema", "h4l_v3_select_report/2"},
                   {"task_id", task.at("task_id")},
                   {"sample", task.at("config").at("sample")},
                   {"kind", task.at("config").at("kind")},
                   {"role", task.at("config").at("role")},
                   {"selection_version", task.at("selection_config").at("version")},
                   {"selection_config_fnv1a64", config_hash(task.at("selection_config"))},
                   {"calibration_payload_fnv1a64", config_hash(task.at("calibration_payload"))},
                   {"calibration", task.at("calibration")},
                   {"lambda_payload_fnv1a64", config_hash(task.at("lambda_payload"))},
                   {"refit_inputs_fnv1a64", config_hash(task.at("refit_inputs"))},
                   {"mela_fnv1a64", config_hash(task.at("mela"))},
                   {"frozen_program", task.value("frozen_program", json())},
                   {"root_version", gROOT->GetVersion()},
                   {"host", std::string(host)},
                   {"inputs_root_sha256_origin", "planned"},
                   {"mela_files", mela_files},
                   {"mela_working_copies", working_copies},
                   {"mela_selftest", "passed"},
                   {"mela_init_root_error", mela_init_errors},
                   {"input_open_root_errors", open_errors},
                   {"inputs", inputs_report},
                   {"source_file_keys", task.at("source_file_keys")},
                   {"counts", counts},
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
