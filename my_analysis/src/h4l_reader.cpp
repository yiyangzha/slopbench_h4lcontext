// h4l_reader: the ROOT (C++) pass over the large input files of the evaluation run.
//
//   h4l_reader OUT.root LIST.txt
//
// LIST.txt: one input per line, tab separated: path, kind (data|mc), file key (unsigned 64-bit), wanted outputs (a
// comma-separated subset of calib, tnp, rec_data, rec_sig, rec_zz, rec_bkg).  For every event: the analysis trigger OR
// and a good vertex; AN loose leptons (muons: |eta| < 2.4, |dxy| < 0.5, |dz| < 1, global or tracker with a matched
// station, ghost cleaning; electrons: |eta| < 2.5, |dxy| < 0.5, |dz| < 1) with raw pT above the record floors 3 / 5 GeV,
// the AN FSR photon of every loose muon (the lowest dR/ET^2), for MC a frozen N(0, 1) deviate per lepton (file key,
// entry, collection, index), the prompt flag (genPartFlav == 1) and the event genWeight.  Outputs (all quantities raw;
// the calibration, the FSR-subtracted isolation of the selection and everything downstream are applied in Python):
//   Calib: the Z -> ll control pair of the lepton calibration (the main analysis's control-pair legs: muon PF, electron
//          mvaFall17V2noIso WP90, SIP < 4, FSR-subtracted isolation < 0.35 with the FSR photons of both legs inside
//          0.01 < dR < 0.3; opposite-sign same flavour, raw mass without FSR in 55-125 GeV, closest to m_Z), with the
//          legs' relative momentum errors (lambda calibration) and their charge x signed-eta x phi slice codes (the
//          residual-width diagnostics of the resolution systematic);
//   TnP:   tag-and-probe pairs (tag: tight ID, single-lepton path fired and matched to its HLT object; probe: any loose
//          lepton of the same flavour and opposite sign at dR > 0.2, mass 55-125 GeV, per tag the probe closest to m_Z);
//   Rec:   loose-lepton event records.  data (rec_data): >= 3 loose leptons.  MC: rec_sig, a loose four-lepton
//          candidate with 85 < m4l < 170 GeV (the window and the m_H morphing); rec_zz, that or a loose candidate with
//          m4l > 65 GeV and a lepton with SIP < 4 that can fail the tight ID or isolation (the 3P1F / 2P2F rows of the
//          prompt subtraction) or a Z + 1 lepton candidate; rec_bkg (DY, ttbar), a loose candidate with m4l > 65 GeV or a
//          Z + 1 lepton candidate.  Loose candidate (the floors of the main analysis's event records): Z1 an opposite-sign
//          same-flavour pair with 30 < m(ll gamma) < 130 GeV, Z2 another same-flavour pair of either charge with
//          3 < m(ll gamma) < 130 GeV, SIP < 8, dR > 0.02, pT 15 / 7 GeV, m(OS) > 2 GeV.  Z + 1 lepton candidate: 3 or 4
//          leptons above 80 % of the analysis thresholds with at most 3 above 125 %, MET < 25 GeV, an opposite-sign
//          same-flavour pair of tight-ID SIP < 4 leptons with pT above 16 / 8 GeV and 36 < m(ll gamma) < 132 GeV;
//   Files: path, entries, generated entries before the production preselection and the sum of genEventSumw (MC).
// Keep the thresholds in sync with h4l_eval/config.py.  Any read error aborts with a non-zero exit code.

#include <TError.h>
#include <TFile.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <Math/ProbFuncMathCore.h>
#include <Math/QuantFuncMathCore.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kMZ = 91.1876, kMuonMass = 0.1056584, kElectronMass = 0.000511;
const std::vector<std::string> kAnalysisOr = {
    "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL",
    "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
    "HLT_TripleMu_12_10_5", "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
    "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
    "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", "HLT_Mu8_DiEle12_CaloIdL_TrackIdL", "HLT_DiMu9_Ele9_CaloIdL_TrackIdL",
    "HLT_Ele25_eta2p1_WPTight_Gsf", "HLT_Ele27_WPTight_Gsf", "HLT_Ele27_eta2p1_WPLoose_Gsf", "HLT_IsoMu20", "HLT_IsoTkMu20",
    "HLT_IsoMu22", "HLT_IsoTkMu22"};
// Tag paths: name, HLT object id, filter bit mask, object pT threshold, max |eta_SC| (electrons).
struct TagPath {
  std::string name;
  int id, bits;
  double pt, eta_max;
};
const std::vector<TagPath> kTagPaths = {{"HLT_IsoMu24", 13, 2, 24.0, 9.0}, {"HLT_IsoTkMu24", 13, 8, 24.0, 9.0},
                                        {"HLT_Ele27_WPTight_Gsf", 11, 2, 27.0, 2.5}, {"HLT_Ele25_eta2p1_WPTight_Gsf", 11, 2, 25.0, 2.1}};

bool g_error = false;
std::string g_message;
void handler(int level, Bool_t abort, const char* location, const char* message) {
  if (level >= kError) {
    g_error = true;
    g_message = std::string(location ? location : "?") + ": " + (message ? message : "");
  }
  DefaultErrorHandler(level, abort, location, message);
}

std::uint64_t splitmix64(std::uint64_t x) {
  std::uint64_t z = x + 0x9E3779B97F4A7C15ULL;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
double object_normal(std::uint64_t key, std::uint64_t entry, std::uint64_t collection, std::uint64_t index) {
  std::uint64_t h = splitmix64(key ^ splitmix64(entry * 4ULL + collection));
  h = splitmix64(h ^ splitmix64(index + 0x5bd1e995ULL));
  const double u = (static_cast<double>(h >> 11) + 0.5) * (1.0 / 9007199254740992.0);
  return ROOT::Math::normal_quantile(u, 1.0);
}

double delta_r(double e1, double p1, double e2, double p2) {
  double dp = std::fmod(p1 - p2 + M_PI, 2 * M_PI);
  if (dp < 0) dp += 2 * M_PI;
  dp -= M_PI;
  return std::hypot(e1 - e2, dp);
}

struct V4 {
  double px = 0, py = 0, pz = 0, e = 0;
  V4() = default;
  V4(double pt, double eta, double phi, double m) {
    px = pt * std::cos(phi);
    py = pt * std::sin(phi);
    pz = pt * std::sinh(eta);
    e = std::sqrt(px * px + py * py + pz * pz + m * m);
  }
  V4 operator+(const V4& o) const {
    V4 r;
    r.px = px + o.px;
    r.py = py + o.py;
    r.pz = pz + o.pz;
    r.e = e + o.e;
    return r;
  }
  double m() const { return std::sqrt(std::max(0.0, e * e - px * px - py * py - pz * pz)); }
};

struct Lepton {
  int flavour = 0, charge = 0, index = 0, lost_hits = 0;
  double pt = 0, eta = 0, eta_sc = 0, phi = 0, rel_err = 0, iso_all = 0, iso_chg = 0, iso_fsr = 0, sip = 0;
  bool pf = false, high_pt_id = false, wpl = false, wp90 = false, tag_id = false, prompt = false;
  double fsr_pt = -1, fsr_eta = 0, fsr_phi = 0, fsr_dret2 = -1, g = 0;
  double mass() const { return flavour == 13 ? kMuonMass : kElectronMass; }
  V4 bare() const { return V4(pt, eta, phi, mass()); }
  V4 dressed() const { return fsr_pt > 0 ? bare() + V4(fsr_pt, fsr_eta, fsr_phi, 0.0) : bare(); }
  double abs_eta_var() const { return flavour == 13 ? std::fabs(eta) : std::fabs(eta_sc); }
  double threshold() const { return flavour == 13 ? 5.0 : 7.0; }
  bool tight_id(double pt_value) const { return flavour == 13 ? (pf || (pt_value > 200.0 && high_pt_id)) : wpl; }
  bool calib_id() const { return flavour == 13 ? pf : wp90; }
  bool tight_iso() const { return tight_id(pt) && sip < 4.0 && iso_fsr < 0.35; }
};

// FSR-subtracted isolation of a control-pair leg with the FSR photons of both legs (main analysis calib_extract).
double pair_isolation(const Lepton& leg, const Lepton& a, const Lepton& b) {
  double photons = 0;
  for (const Lepton* owner : {&a, &b}) {
    if (!(owner->fsr_pt > 0)) continue;
    const double dr = delta_r(leg.eta, leg.phi, owner->fsr_eta, owner->fsr_phi);
    if (dr > 0.01 && dr < 0.3) photons += owner->fsr_pt;
  }
  return leg.iso_chg + std::max(0.0, (leg.iso_all - leg.iso_chg) - photons / leg.pt);
}

// Charge x signed-eta x phi slice code of a control-pair leg: (charge > 0) * 144 + eta slice * 12 + phi slice, 12 slices
// of the signed lepton eta in [-2.4, 2.4] (muons) or [-2.5, 2.5] (electrons) and of phi in [-pi, pi].
int slice_code(const Lepton& l) {
  const double eta_max = l.flavour == 13 ? 2.4 : 2.5;
  const int ie = std::clamp(static_cast<int>(std::floor((l.eta + eta_max) / (2 * eta_max / 12))), 0, 11);
  const int ip = std::clamp(static_cast<int>(std::floor((l.phi + M_PI) / (2 * M_PI / 12))), 0, 11);
  return (l.charge > 0 ? 144 : 0) + ie * 12 + ip;
}

template <typename T>
using Arr = TTreeReaderArray<T>;
template <typename T>
using Val = TTreeReaderValue<T>;

struct Out {
  TTree calib{"Calib", "Z->ll control pairs"}, tnp{"TnP", "tag-and-probe pairs"}, rec{"Rec", "event records"},
      files{"Files", "inputs"};
  // Calib
  Char_t c_fl = 0;
  Float_t c_m = 0, c_pt1 = 0, c_eta1 = 0, c_pt2 = 0, c_eta2 = 0, c_g1 = 0, c_g2 = 0, c_w = 1, c_rel1 = 0, c_rel2 = 0;
  UShort_t c_code1 = 0, c_code2 = 0;
  // TnP
  Char_t t_fl = 0;
  Bool_t t_pass = false, t_prompt = false;
  Float_t t_m = 0, t_tpt = 0, t_teta = 0, t_ppt = 0, t_peta = 0, t_tg = 0, t_pg = 0, t_w = 1;
  // Rec
  static constexpr int kMax = 32;
  Int_t r_n = 0, r_file = 0;
  Float_t r_met = 0, r_w = 1;
  Long64_t r_entry = 0;
  Char_t r_fl[kMax], r_ch[kMax], r_lost[kMax];
  Bool_t r_pf[kMax], r_hpt[kMax], r_wpl[kMax];
  Float_t r_pt[kMax], r_eta[kMax], r_etasc[kMax], r_phi[kMax], r_relerr[kMax], r_isoall[kMax], r_isochg[kMax], r_sip[kMax],
      r_fpt[kMax], r_feta[kMax], r_fphi[kMax], r_g[kMax];
  // Files
  std::string f_path;
  Long64_t f_entries = 0, f_npre = 0;
  Double_t f_sumw = 0;

  Out() {
    calib.Branch("flavour", &c_fl, "flavour/B");
    for (auto [n, a] : std::initializer_list<std::pair<const char*, Float_t*>>{
             {"mass", &c_m}, {"pt1", &c_pt1}, {"eta1", &c_eta1}, {"pt2", &c_pt2}, {"eta2", &c_eta2}, {"g1", &c_g1}, {"g2", &c_g2},
             {"w", &c_w}, {"rel1", &c_rel1}, {"rel2", &c_rel2}})
      calib.Branch(n, a, (std::string(n) + "/F").c_str());
    calib.Branch("code1", &c_code1, "code1/s");
    calib.Branch("code2", &c_code2, "code2/s");
    tnp.Branch("flavour", &t_fl, "flavour/B");
    tnp.Branch("pass", &t_pass, "pass/O");
    tnp.Branch("prompt", &t_prompt, "prompt/O");
    for (auto [n, a] : std::initializer_list<std::pair<const char*, Float_t*>>{{"mass", &t_m}, {"tag_pt", &t_tpt}, {"tag_eta", &t_teta},
                                                                              {"probe_pt", &t_ppt}, {"probe_eta", &t_peta},
                                                                              {"tag_g", &t_tg}, {"probe_g", &t_pg}, {"w", &t_w}})
      tnp.Branch(n, a, (std::string(n) + "/F").c_str());
    rec.Branch("n_lep", &r_n, "n_lep/I");
    rec.Branch("e_file", &r_file, "e_file/I");
    rec.Branch("e_met", &r_met, "e_met/F");
    rec.Branch("e_w", &r_w, "e_w/F");
    rec.Branch("e_entry", &r_entry, "e_entry/L");
    rec.Branch("l_flavour", r_fl, "l_flavour[n_lep]/B");
    rec.Branch("l_charge", r_ch, "l_charge[n_lep]/B");
    rec.Branch("l_lost_hits", r_lost, "l_lost_hits[n_lep]/B");
    for (auto [n, a] : std::initializer_list<std::pair<const char*, Bool_t*>>{{"l_pf", r_pf}, {"l_high_pt_id", r_hpt}, {"l_wpl", r_wpl}})
      rec.Branch(n, a, (std::string(n) + "[n_lep]/O").c_str());
    for (auto [n, a] : std::initializer_list<std::pair<const char*, Float_t*>>{
             {"l_pt_raw", r_pt}, {"l_eta", r_eta}, {"l_eta_sc", r_etasc}, {"l_phi", r_phi}, {"l_rel_err", r_relerr},
             {"l_iso_all", r_isoall}, {"l_iso_chg", r_isochg}, {"l_sip", r_sip}, {"l_fsr_pt", r_fpt}, {"l_fsr_eta", r_feta},
             {"l_fsr_phi", r_fphi}, {"l_g", r_g}})
      rec.Branch(n, a, (std::string(n) + "[n_lep]/F").c_str());
    files.Branch("path", &f_path);
    files.Branch("entries", &f_entries, "entries/L");
    files.Branch("n_preselection", &f_npre, "n_preselection/L");
    files.Branch("sumw", &f_sumw, "sumw/D");
  }
};

long long generated_entries(TFile& f) {
  if (auto* named = dynamic_cast<TNamed*>(f.Get("PFnanoFinalMCProvenance"))) {
    const std::string title = named->GetTitle();
    for (const char* key : {"\"entries_before_selection\"", "\"n_preselection\"", "\"genEventCount\""}) {
      auto pos = title.find(key);
      if (pos == std::string::npos) continue;
      pos = title.find(':', pos);
      return std::stoll(title.substr(pos + 1));
    }
  }
  if (auto* runs = dynamic_cast<TTree*>(f.Get("Runs"))) {
    if (runs->GetBranch("genEventCount")) {
      Long64_t count = 0;
      long long total = 0;
      runs->SetBranchAddress("genEventCount", &count);
      for (Long64_t i = 0; i < runs->GetEntries(); ++i) {
        if (runs->GetEntry(i) <= 0) throw std::runtime_error("cannot read Runs");
        total += count;
      }
      runs->ResetBranchAddresses();
      return total;
    }
  }
  throw std::runtime_error("no generated-entry bookkeeping");
}

// Sum of genEventSumw over the Runs entries (0 when the file carries none: the count-based normalization is used then).
double generated_sumw(TFile& f) {
  auto* runs = dynamic_cast<TTree*>(f.Get("Runs"));
  if (!runs || !runs->GetBranch("genEventSumw")) return 0.0;
  Double_t sumw = 0;
  double total = 0;
  runs->SetBranchAddress("genEventSumw", &sumw);
  for (Long64_t i = 0; i < runs->GetEntries(); ++i) {
    if (runs->GetEntry(i) <= 0) throw std::runtime_error("cannot read Runs");
    total += sumw;
  }
  runs->ResetBranchAddresses();
  return total;
}

// A loose four-lepton candidate (the floors of the main analysis's event records, raw pT): Z1 an opposite-sign same-flavour
// pair with 30 < m(ll gamma) < 130, Z2 another same-flavour pair of either charge with 3 < m(ll gamma) < 130, all four
// leptons with SIP < 8, dR > 0.02 between them, leading / subleading pT above 15 / 7 GeV, every opposite-sign pair above
// 2 GeV without FSR, m4l in [lo, hi].  need_failing: one of the Z2 leptons has SIP < 4 and can fail the tight ID or the
// isolation (tight ID failed or the plain isolation above 0.35, an upper bound of the FSR-subtracted one).
bool loose_candidate(const std::vector<Lepton>& l, double lo, double hi, bool need_failing) {
  const int n = static_cast<int>(l.size());
  auto can_fail = [&](const Lepton& x) { return x.sip < 4.0 && (!x.tight_id(x.pt) || x.iso_all >= 0.35); };
  for (int a = 0; a < n; ++a)
    for (int b = a + 1; b < n; ++b) {
      if (l[a].flavour != l[b].flavour || l[a].charge == l[b].charge || l[a].sip >= 8 || l[b].sip >= 8) continue;
      const V4 z1 = l[a].dressed() + l[b].dressed();
      const double m1 = z1.m();
      if (m1 <= 30 || m1 >= 130) continue;
      for (int c = 0; c < n; ++c)
        for (int d = c + 1; d < n; ++d) {
          if (c == a || c == b || d == a || d == b) continue;
          if (l[c].flavour != l[d].flavour || l[c].sip >= 8 || l[d].sip >= 8) continue;
          if (need_failing && !can_fail(l[c]) && !can_fail(l[d])) continue;
          const V4 z2 = l[c].dressed() + l[d].dressed();
          const double m2 = z2.m();
          if (m2 <= 3 || m2 >= 130) continue;
          const int idx[4] = {a, b, c, d};
          bool pass = true;
          double pts[4];
          for (int u = 0; u < 4 && pass; ++u) {
            pts[u] = l[idx[u]].pt;
            for (int v = u + 1; v < 4 && pass; ++v) {
              const Lepton &x = l[idx[u]], &y = l[idx[v]];
              if (delta_r(x.eta, x.phi, y.eta, y.phi) <= 0.02) pass = false;
              if (x.charge != y.charge && (x.bare() + y.bare()).m() <= 2.0) pass = false;
            }
          }
          if (!pass) continue;
          std::sort(pts, pts + 4, [](double p, double q) { return p > q; });
          if (pts[0] <= 15.0 || pts[1] <= 7.0) continue;
          const double m4 = (z1 + z2).m();
          if (m4 > lo && m4 < hi) return true;
        }
    }
  return false;
}

// A Z + 1 lepton candidate for the fake-rate rows (margins around the selection, which is applied downstream on the
// calibrated leptons): 3 or 4 leptons above 80 % of the analysis thresholds, at most 3 above 125 % of them (an electron
// within dR < 0.05 of a muon not counted there), MET < 25 GeV, an opposite-sign same-flavour pair of tight-ID leptons with
// SIP < 4, pT above 16 / 8 GeV and 36 < m(ll gamma) < 132 GeV.
bool zl_candidate(const std::vector<Lepton>& l, double met) {
  if (!(met < 25.0)) return false;
  const int n = static_cast<int>(l.size());
  int above = 0, sure = 0;
  for (int k = 0; k < n; ++k) {
    if (l[k].pt > 0.8 * l[k].threshold()) ++above;
    if (l[k].pt > 1.25 * l[k].threshold()) {
      bool overlap = false;
      if (l[k].flavour == 11)
        for (int j = 0; j < n; ++j)
          if (l[j].flavour == 13 && delta_r(l[k].eta, l[k].phi, l[j].eta, l[j].phi) < 0.05) overlap = true;
      if (!overlap) ++sure;
    }
  }
  if (above < 3 || above > 4 || sure > 3) return false;
  for (int a = 0; a < n; ++a)
    for (int b = a + 1; b < n; ++b) {
      const Lepton &x = l[a], &y = l[b];
      if (x.flavour != y.flavour || x.charge == y.charge) continue;
      if (!x.tight_id(x.pt) || !y.tight_id(y.pt) || x.sip >= 4.0 || y.sip >= 4.0) continue;
      if (std::max(x.pt, y.pt) <= 16.0 || std::min(x.pt, y.pt) <= 8.0) continue;
      const double m = (x.dressed() + y.dressed()).m();
      if (m > 36.0 && m < 132.0) return true;
    }
  return false;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: h4l_reader OUT.root LIST.txt\n";
    return 64;
  }
  gROOT->SetBatch(true);
  SetErrorHandler(handler);
  try {
    std::ifstream list(argv[2]);
    std::vector<std::array<std::string, 4>> inputs;
    std::string line;
    while (std::getline(list, line)) {
      if (line.empty()) continue;
      std::istringstream ss(line);
      std::array<std::string, 4> f;
      for (auto& x : f) std::getline(ss, x, '\t');
      inputs.push_back(f);
    }
    std::unique_ptr<TFile> out_file(TFile::Open(argv[1], "RECREATE", "", 404));  // LZ4: fast to read back
    if (!out_file || out_file->IsZombie()) throw std::runtime_error(std::string("cannot create ") + argv[1]);
    out_file->cd();
    auto* out = new Out();
    for (std::size_t fi = 0; fi < inputs.size(); ++fi) {
      const auto& [path, kind, key_text, wants] = inputs[fi];
      const bool is_mc = kind == "mc";
      const std::uint64_t key = std::stoull(key_text);
      const bool w_calib = wants.find("calib") != std::string::npos, w_tnp = wants.find("tnp") != std::string::npos,
                 w_rdata = wants.find("rec_data") != std::string::npos, w_rsig = wants.find("rec_sig") != std::string::npos,
                 w_rzz = wants.find("rec_zz") != std::string::npos, w_rbkg = wants.find("rec_bkg") != std::string::npos;
      g_error = false;
      std::unique_ptr<TFile> in(TFile::Open(path.c_str(), "READ"));
      if (!in || in->IsZombie()) throw std::runtime_error("cannot open " + path);
      auto* tree = dynamic_cast<TTree*>(in->Get("Events"));
      if (!tree) throw std::runtime_error("no Events in " + path);
      out->f_path = path;
      out->f_entries = tree->GetEntries();
      out->f_npre = is_mc ? generated_entries(*in) : 0;
      out->f_sumw = is_mc ? generated_sumw(*in) : 0.0;
      TTreeReader r(tree);
      Arr<Float_t> mu_pt(r, "Muon_pt"), mu_eta(r, "Muon_eta"), mu_phi(r, "Muon_phi"), mu_dxy(r, "Muon_dxy"), mu_dz(r, "Muon_dz"),
          mu_sip(r, "Muon_sip3d"), mu_iso(r, "Muon_pfRelIso03_all"), mu_isochg(r, "Muon_pfRelIso03_chg"), mu_pterr(r, "Muon_ptErr"),
          mu_iso04(r, "Muon_pfRelIso04_all");
      Arr<Int_t> mu_charge(r, "Muon_charge"), mu_nst(r, "Muon_nStations");
      Arr<Bool_t> mu_global(r, "Muon_isGlobal"), mu_tracker(r, "Muon_isTracker"), mu_pf(r, "Muon_isPFcand"), mu_tight(r, "Muon_tightId");
      Arr<UChar_t> mu_hpt(r, "Muon_highPtId");
      Arr<Float_t> el_pt(r, "Electron_pt"), el_eta(r, "Electron_eta"), el_desc(r, "Electron_deltaEtaSC"), el_phi(r, "Electron_phi"),
          el_dxy(r, "Electron_dxy"), el_dz(r, "Electron_dz"), el_sip(r, "Electron_sip3d"), el_iso(r, "Electron_pfRelIso03_all"),
          el_isochg(r, "Electron_pfRelIso03_chg"), el_eerr(r, "Electron_energyErr");
      Arr<Int_t> el_charge(r, "Electron_charge"), el_cb(r, "Electron_cutBased");
      Arr<Bool_t> el_wpl(r, "Electron_mvaFall17V2noIso_WPL"), el_wp90(r, "Electron_mvaFall17V2noIso_WP90");
      Arr<UChar_t> el_lost(r, "Electron_lostHits");
      Val<Float_t> met(r, "MET_pt");
      Val<Int_t> npv(r, "PV_npvsGood");
      std::unique_ptr<Val<Float_t>> gen_weight;
      if (is_mc && tree->GetBranch("genWeight")) gen_weight = std::make_unique<Val<Float_t>>(r, "genWeight");
      const bool has_fsr = tree->GetBranch("FsrPhoton_pt") != nullptr;
      std::unique_ptr<Arr<Float_t>> f_pt, f_eta, f_phi, f_iso, f_dret2;
      std::unique_ptr<Arr<Int_t>> f_mu;
      if (has_fsr) {
        f_pt = std::make_unique<Arr<Float_t>>(r, "FsrPhoton_pt");
        f_eta = std::make_unique<Arr<Float_t>>(r, "FsrPhoton_eta");
        f_phi = std::make_unique<Arr<Float_t>>(r, "FsrPhoton_phi");
        f_iso = std::make_unique<Arr<Float_t>>(r, "FsrPhoton_relIso03");
        f_dret2 = std::make_unique<Arr<Float_t>>(r, "FsrPhoton_dROverEt2");
        f_mu = std::make_unique<Arr<Int_t>>(r, "FsrPhoton_muonIdx");
      }
      const bool has_to = w_tnp && tree->GetBranch("TrigObj_pt") != nullptr;
      std::unique_ptr<Arr<Float_t>> to_pt, to_eta, to_phi;
      std::unique_ptr<Arr<Int_t>> to_id, to_bits;
      if (has_to) {
        to_pt = std::make_unique<Arr<Float_t>>(r, "TrigObj_pt");
        to_eta = std::make_unique<Arr<Float_t>>(r, "TrigObj_eta");
        to_phi = std::make_unique<Arr<Float_t>>(r, "TrigObj_phi");
        to_id = std::make_unique<Arr<Int_t>>(r, "TrigObj_id");
        to_bits = std::make_unique<Arr<Int_t>>(r, "TrigObj_filterBits");
      }
      std::unique_ptr<Arr<UChar_t>> mu_flav, el_flav;
      if (is_mc && tree->GetBranch("Muon_genPartFlav")) mu_flav = std::make_unique<Arr<UChar_t>>(r, "Muon_genPartFlav");
      if (is_mc && tree->GetBranch("Electron_genPartFlav")) el_flav = std::make_unique<Arr<UChar_t>>(r, "Electron_genPartFlav");
      std::vector<std::unique_ptr<Val<Bool_t>>> paths;
      for (const auto& p : kAnalysisOr)
        if (tree->GetBranch(p.c_str())) paths.push_back(std::make_unique<Val<Bool_t>>(r, p.c_str()));
      if (paths.empty()) throw std::runtime_error("no analysis trigger path in " + path);
      std::vector<std::unique_ptr<Val<Bool_t>>> tag_fired(kTagPaths.size());
      if (w_tnp)
        for (std::size_t k = 0; k < kTagPaths.size(); ++k)
          if (tree->GetBranch(kTagPaths[k].name.c_str())) tag_fired[k] = std::make_unique<Val<Bool_t>>(r, kTagPaths[k].name.c_str());

      Long64_t entry = -1;
      while (r.Next()) {
        ++entry;
        if (g_error) throw std::runtime_error("ROOT error while reading " + path + ": " + g_message);
        bool trig = false;
        for (auto& p : paths) trig = trig || **p;
        if (!trig || *npv < 1) continue;
        const double weight = gen_weight ? static_cast<double>(**gen_weight) : 1.0;
        // Loose muons.
        std::vector<Lepton> muons;
        for (std::size_t i = 0; i < mu_pt.GetSize(); ++i) {
          if (!(mu_pt[i] > 3.0) || std::fabs(mu_eta[i]) >= 2.4 || std::fabs(mu_dxy[i]) >= 0.5 || std::fabs(mu_dz[i]) >= 1.0) continue;
          if (!(mu_global[i] || (mu_tracker[i] && mu_nst[i] > 0))) continue;
          if (!(mu_pterr[i] > 0) || !std::isfinite(mu_sip[i]) || !std::isfinite(mu_iso[i]) || !std::isfinite(mu_isochg[i])) continue;
          Lepton l;
          l.flavour = 13;
          l.charge = mu_charge[i];
          l.index = static_cast<int>(i);
          l.pt = mu_pt[i];
          l.eta = l.eta_sc = mu_eta[i];
          l.phi = mu_phi[i];
          l.rel_err = mu_pterr[i] / mu_pt[i];
          l.iso_all = mu_iso[i];
          l.iso_chg = mu_isochg[i];
          l.sip = mu_sip[i];
          l.pf = mu_pf[i];
          l.high_pt_id = mu_hpt[i] >= 1;
          l.tag_id = mu_tight[i] && mu_iso04[i] < 0.15;
          l.prompt = mu_flav && (*mu_flav)[i] == 1;
          l.g = is_mc ? object_normal(key, static_cast<std::uint64_t>(entry), 0, i) : 0.0;
          if (has_fsr) {
            for (std::size_t k = 0; k < f_pt->GetSize(); ++k) {
              if ((*f_mu)[k] != static_cast<int>(i)) continue;
              if (!((*f_pt)[k] > 2.0) || std::fabs((*f_eta)[k]) >= 2.4 || (*f_iso)[k] >= 1.8 || (*f_dret2)[k] >= 0.012) continue;
              if (delta_r(l.eta, l.phi, (*f_eta)[k], (*f_phi)[k]) >= 0.5) continue;
              if (l.fsr_dret2 < 0 || (*f_dret2)[k] < l.fsr_dret2) {
                l.fsr_pt = (*f_pt)[k];
                l.fsr_eta = (*f_eta)[k];
                l.fsr_phi = (*f_phi)[k];
                l.fsr_dret2 = (*f_dret2)[k];
              }
            }
          }
          muons.push_back(l);
        }
        // Ghost cleaning: dR < 0.02 keeps the PF (else the higher-pT) muon; a tracker-only muon within dR < 0.05 of a
        // same-charge PF muon is a ghost.
        std::vector<char> ghost(muons.size(), 0);
        for (std::size_t a = 0; a < muons.size(); ++a)
          for (std::size_t b = a + 1; b < muons.size(); ++b) {
            const double dr = delta_r(muons[a].eta, muons[a].phi, muons[b].eta, muons[b].phi);
            if (dr < 0.02 && !ghost[a] && !ghost[b]) {
              const bool keep_a = muons[a].pf != muons[b].pf ? muons[a].pf : muons[a].pt >= muons[b].pt;
              ghost[keep_a ? b : a] = 1;
            }
            if (dr < 0.05 && muons[a].charge == muons[b].charge) {
              const bool trk_a = mu_tracker[muons[a].index] && !mu_global[muons[a].index];
              const bool trk_b = mu_tracker[muons[b].index] && !mu_global[muons[b].index];
              if (trk_a && !muons[a].pf && muons[b].pf) ghost[a] = 1;
              if (trk_b && !muons[b].pf && muons[a].pf) ghost[b] = 1;
            }
          }
        std::vector<Lepton> lep;
        for (std::size_t a = 0; a < muons.size(); ++a)
          if (!ghost[a]) lep.push_back(muons[a]);
        for (std::size_t i = 0; i < el_pt.GetSize(); ++i) {
          if (!(el_pt[i] > 5.0) || std::fabs(el_eta[i]) >= 2.5 || std::fabs(el_dxy[i]) >= 0.5 || std::fabs(el_dz[i]) >= 1.0) continue;
          if (!(el_eerr[i] > 0) || !std::isfinite(el_sip[i]) || !std::isfinite(el_iso[i]) || !std::isfinite(el_isochg[i])) continue;
          Lepton l;
          l.flavour = 11;
          l.charge = el_charge[i];
          l.index = static_cast<int>(i);
          l.lost_hits = el_lost[i];
          l.pt = el_pt[i];
          l.eta = el_eta[i];
          l.eta_sc = el_eta[i] + el_desc[i];
          l.phi = el_phi[i];
          l.rel_err = el_eerr[i] / (el_pt[i] * std::cosh(el_eta[i]));
          l.iso_all = el_iso[i];
          l.iso_chg = el_isochg[i];
          l.sip = el_sip[i];
          l.wpl = el_wpl[i];
          l.wp90 = el_wp90[i];
          l.tag_id = el_cb[i] >= 4;
          l.prompt = el_flav && (*el_flav)[i] == 1;
          l.g = is_mc ? object_normal(key, static_cast<std::uint64_t>(entry), 1, i) : 0.0;
          lep.push_back(l);
        }
        if (lep.size() < 2) continue;
        // FSR-subtracted isolation of the tag-and-probe definitions (photons of the loose muons with SIP < 4).
        for (auto& l : lep) {
          double photons = 0;
          for (const auto& o : lep) {
            if (o.flavour != 13 || !(o.fsr_pt > 0) || o.sip >= 4.0) continue;
            const double dr = delta_r(l.eta, l.phi, o.fsr_eta, o.fsr_phi);
            const bool veto_ok = l.flavour == 13 ? dr > 0.01 : (std::fabs(l.eta_sc) < 1.479 || dr > 0.08);
            if (dr < 0.3 && veto_ok) photons += o.fsr_pt;
          }
          l.iso_fsr = l.iso_chg + std::max(0.0, (l.iso_all - l.iso_chg) - photons / l.pt);
        }
        const int n = static_cast<int>(lep.size());
        // Calibration pair.
        if (w_calib) {
          int ba = -1, bb = -1;
          double bm = 0, best = 1e9;
          for (int a = 0; a < n; ++a)
            for (int b = a + 1; b < n; ++b) {
              if (lep[a].flavour != lep[b].flavour || lep[a].charge == lep[b].charge) continue;
              if (!lep[a].calib_id() || !lep[b].calib_id() || lep[a].sip >= 4.0 || lep[b].sip >= 4.0) continue;
              if (pair_isolation(lep[a], lep[a], lep[b]) >= 0.35 || pair_isolation(lep[b], lep[a], lep[b]) >= 0.35) continue;
              const double m = (lep[a].bare() + lep[b].bare()).m();
              if (m <= 55 || m >= 125) continue;
              if (std::fabs(m - kMZ) < best) {
                best = std::fabs(m - kMZ);
                ba = a;
                bb = b;
                bm = m;
              }
            }
          if (ba >= 0) {
            out->c_fl = static_cast<Char_t>(lep[ba].flavour);
            out->c_m = bm;
            out->c_pt1 = lep[ba].pt;
            out->c_eta1 = lep[ba].abs_eta_var();
            out->c_pt2 = lep[bb].pt;
            out->c_eta2 = lep[bb].abs_eta_var();
            out->c_g1 = lep[ba].g;
            out->c_g2 = lep[bb].g;
            out->c_w = weight;
            out->c_rel1 = lep[ba].rel_err;
            out->c_rel2 = lep[bb].rel_err;
            out->c_code1 = static_cast<UShort_t>(slice_code(lep[ba]));
            out->c_code2 = static_cast<UShort_t>(slice_code(lep[bb]));
            out->calib.Fill();
          }
        }
        // Tag and probe.
        if (w_tnp && has_to) {
          for (int a = 0; a < n; ++a) {
            const Lepton& t = lep[a];
            if (!t.tag_id) continue;
            if (t.flavour == 13 && !(t.pt > 26.0 && std::fabs(t.eta) < 2.4)) continue;
            if (t.flavour == 11) {
              const double s = std::fabs(t.eta_sc);
              if (!(t.pt > 30.0 && s < 2.1 && !(s > 1.4442 && s < 1.566))) continue;
            }
            bool matched = false;
            for (std::size_t k = 0; k < kTagPaths.size() && !matched; ++k) {
              const auto& tp = kTagPaths[k];
              if (tp.id != t.flavour || !tag_fired[k] || !**tag_fired[k]) continue;
              if (t.flavour == 11 && std::fabs(t.eta_sc) >= tp.eta_max) continue;
              for (std::size_t o = 0; o < to_pt->GetSize(); ++o) {
                if ((*to_id)[o] != tp.id || ((*to_bits)[o] & tp.bits) == 0 || (*to_pt)[o] < tp.pt) continue;
                if (delta_r(t.eta, t.phi, (*to_eta)[o], (*to_phi)[o]) < 0.1) {
                  matched = true;
                  break;
                }
              }
            }
            if (!matched) continue;
            int bp = -1;
            double bm = 0, best = 1e9;
            for (int b = 0; b < n; ++b) {
              if (b == a || lep[b].flavour != t.flavour || lep[b].charge == t.charge) continue;
              if (delta_r(t.eta, t.phi, lep[b].eta, lep[b].phi) <= 0.2) continue;
              const double m = (t.bare() + lep[b].bare()).m();
              if (m <= 55 || m >= 125) continue;
              if (std::fabs(m - kMZ) < best) {
                best = std::fabs(m - kMZ);
                bp = b;
                bm = m;
              }
            }
            if (bp < 0) continue;
            const Lepton& p = lep[bp];
            out->t_fl = static_cast<Char_t>(t.flavour);
            out->t_m = bm;
            out->t_tpt = t.pt;
            out->t_teta = t.abs_eta_var();
            out->t_ppt = p.pt;
            out->t_peta = p.abs_eta_var();
            out->t_pass = p.tight_iso();
            out->t_prompt = t.prompt && p.prompt;
            out->t_tg = t.g;
            out->t_pg = p.g;
            out->t_w = weight;
            out->tnp.Fill();
          }
        }
        // Event records.
        bool keep = false;
        if (w_rdata && n >= 3) keep = true;
        if (!keep && w_rsig && n >= 4) keep = loose_candidate(lep, 85.0, 170.0, false);
        if (!keep && w_rzz) {
          if (n >= 4) keep = loose_candidate(lep, 85.0, 170.0, false) || loose_candidate(lep, 65.0, 1e9, true);
          if (!keep && n >= 3) keep = zl_candidate(lep, *met);
        }
        if (!keep && w_rbkg) {
          if (n >= 4) keep = loose_candidate(lep, 65.0, 1e9, false);
          if (!keep && n >= 3) keep = zl_candidate(lep, *met);
        }
        if (keep) {
          const int m = std::min(n, static_cast<int>(Out::kMax));
          out->r_n = m;
          out->r_file = static_cast<Int_t>(fi);
          out->r_met = *met;
          out->r_w = weight;
          out->r_entry = entry;
          for (int k = 0; k < m; ++k) {
            const Lepton& l = lep[k];
            out->r_fl[k] = static_cast<Char_t>(l.flavour);
            out->r_ch[k] = static_cast<Char_t>(l.charge);
            out->r_lost[k] = static_cast<Char_t>(std::min(l.lost_hits, 127));
            out->r_pf[k] = l.pf;
            out->r_hpt[k] = l.high_pt_id;
            out->r_wpl[k] = l.wpl;
            out->r_pt[k] = l.pt;
            out->r_eta[k] = l.eta;
            out->r_etasc[k] = l.eta_sc;
            out->r_phi[k] = l.phi;
            out->r_relerr[k] = l.rel_err;
            out->r_isoall[k] = l.iso_all;
            out->r_isochg[k] = l.iso_chg;
            out->r_sip[k] = l.sip;
            out->r_fpt[k] = l.fsr_pt;
            out->r_feta[k] = l.fsr_eta;
            out->r_fphi[k] = l.fsr_phi;
            out->r_g[k] = l.g;
          }
          out->rec.Fill();
        }
      }
      if (r.GetEntryStatus() != TTreeReader::kEntryNotFound && r.GetEntryStatus() != TTreeReader::kEntryValid &&
          r.GetEntryStatus() != TTreeReader::kEntryBeyondEnd)
        throw std::runtime_error("read error in " + path);
      if (entry + 1 != out->f_entries) throw std::runtime_error("read " + std::to_string(entry + 1) + " of " + std::to_string(out->f_entries) + " entries of " + path);
      out->files.Fill();
    }
    out_file->cd();
    for (TTree* t : {&out->calib, &out->tnp, &out->rec, &out->files}) {
      if (t->Write() <= 0 && t->GetEntries() > 0) throw std::runtime_error("cannot write tree");
    }
    out_file->Close();
    if (g_error) throw std::runtime_error("ROOT error: " + g_message);
  } catch (const std::exception& e) {
    std::cerr << "h4l_reader: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
