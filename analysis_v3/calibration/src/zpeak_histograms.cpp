// zpeak_histograms: Z -> ll mass histograms and template pairs of the
// lepton calibration (factorized model, h4l/calibration.h).
//
// Usage: zpeak_histograms --job JOB.json --out OUT.root
//
// The job names CalibPairs files of one sample (calib_extract outputs, as
// {path, calib_pairs, root_sha256} from the extraction scan), the
// calibration configuration, the role ("data" or "mc") and the optional
// transformations.  For every pair, each leg's pT is transformed in this
// order:
//   inject  (closure tests and response passes): pT *= (1 + s(pT, |eta|)) *
//           (1 + r(pT, |eta|) * N_inject), s and r the documented test
//           functions, N_inject the per-lepton deviate of the inject stream;
//   payload: role data: pT *= exp(-u(pT, |eta|)); role mc:
//           pT *= (1 + sqrt(max(v(pT, |eta|), 0)) * N), N the per-lepton
//           deviate of stream 0 (factorized payload of h4l/calibration.h).
// |eta| is abs(eta) for muons and abs(eta_SC) for electrons.  The tight
// control-pair selection (FSR-subtracted isolation) is then applied to the
// transformed legs and the pair mass is rescaled by
// sqrt(pT1' pT2' / (pT1 pT2)) (exact for massless legs).  Every pair enters
// two category families of its transformed legs:
//   family A: the unordered pair (i <= j) of fine |eta| bins, only if both
//             legs have pT >= family_a_min_pt;
//   family B: the unordered pair (k <= l) of (pT bin, region) bins.
// Per category <F>_<f>_<i>_<j> (F = A or B, f = mm or ee):
//   m<F>     the mass;
//   md<F>    ("delta_smear") the mass smeared once by the common relative
//            delta of the template fit, m (1 + delta eps_d), eps_d the pair
//            deviate of stream 2;
//   w<F>     the leg composition of the pairs with 80 < m < 100 GeV: for A
//            the pT bins of the leg in bin i (entries 0 .. n_pt-1) and of the
//            leg in bin j (n_pt .. 2 n_pt-1; for i = j both legs in the first
//            half); for B the fine |eta| bins of the leg in bin k and in bin l;
//   lpt<F>, ppt<F>   both legs' pT and the pair pT.
// Per flavour, for the legs of pairs with 80 < m < 100 GeV:
//   leg_<f>_{w,wpt,w2,ws,wr2}   per (pT bin, region) bin: sums of weights, of
//                               w pT, of w^2 and of the injected functions;
//   legeta_<f>_{w,weta,ws,wr2}  per fine |eta| bin.
// Trees: TemplatePairs ("write_template_pairs", the MC role): flavour, the
// family A and B categories (-1 if none), mass, weight and the pair deviate
// of stream 3; FullPairs ("write_full_pairs"): the transformed legs of every
// selected pair (pT, signed eta, phi, charge, ...) for the per-event
// mass-uncertainty calibration and diagnostics.
// Weights: MC genWeight * weight_scale; data 1.  "half" (0 or 1) keeps only
// the rows whose original file key has that parity.  Every input must be
// distinct and hold exactly its scanned number of rows.
// Publication: OUT.root.lock is created exclusively first (exit code 17 if it
// or OUT.root exists); the output is written to an attempt-qualified
// temporary file, closed, reopened and checked with bounded retries (the
// XRootD view of a fresh FUSE write can lag) and renamed onto OUT.root.  The
// lock file is kept as the record of the run that produced the output.

#include "h4l/calibration.h"
#include "h4l/hash.h"
#include "h4l/io.h"
#include "h4l/root_io.h"

#include <TError.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TTree.h>

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

namespace {
using h4l::json;

bool g_read_error = false;
std::string g_read_message;
void record_errors(int level, Bool_t abort, const char* location, const char* message) {
  if (level >= kError) {
    g_read_error = true;
    g_read_message = std::string(location ? location : "?") + ": " + (message ? message : "");
  }
  DefaultErrorHandler(level, abort, location, message);
}

// Closure-test functions of (pT, |eta|): a + b (|eta| - eta0) + c (pT - pt0) / pt0 + d (|eta| - eta0)^2.
struct TestFunction {
  double a = 0, b = 0, c = 0, d = 0, pt0 = 45, eta0 = 1.2;
  double operator()(double pt, double abs_eta) const {
    const double x = abs_eta - eta0;
    return a + b * x + c * (pt - pt0) / pt0 + d * x * x;
  }
  static TestFunction from_json(const json& node) {
    TestFunction f;
    f.a = node.value("a", 0.0);
    f.b = node.value("b_eta", 0.0);
    f.c = node.value("c_pt", 0.0);
    f.d = node.value("d_eta2", 0.0);
    f.pt0 = node.value("pt0", 45.0);
    f.eta0 = node.value("eta0", 1.2);
    return f;
  }
};

struct Flavour {
  std::string name;  // "muon" or "electron"
  std::string tag;   // "mm" or "ee"
  h4l::Collection collection = h4l::Collection::Muon;
  h4l::CalibModel model;
  double min_pt = 0, max_abs_eta = 0, max_iso = 0, max_sip = 0;
  std::string id;
  bool inject = false;
  TestFunction inject_scale, inject_smear;
  h4l::FactorizedPayload payload;

  bool pass(double pt, double abs_eta_lepton, float iso_fsr, float sip, unsigned flags) const {
    if (pt < min_pt || abs_eta_lepton >= max_abs_eta || iso_fsr >= max_iso || sip >= max_sip) return false;
    return collection == h4l::Collection::Muon ? h4l::muon_id(id, flags) : h4l::electron_id(id, flags);
  }
};

struct LegIn {
  Short_t index = 0;
  Float_t pt = 0, eta = 0, eta_sc = 0, phi = 0, pt_err = 0, iso_fsr = 0, sip = 0, r9 = -9, gen_pt = -1;
  UShort_t flags = 0;
  Char_t charge = 0;
  UChar_t gen_flav = 0;
};

std::unique_ptr<TH1D> make_h1(const std::string& name, int bins, double low, double high) {
  auto histogram = std::make_unique<TH1D>(name.c_str(), name.c_str(), bins, low, high);
  histogram->SetDirectory(nullptr);
  histogram->Sumw2();
  return histogram;
}

// Deviate of one pair (the two original lepton indices) in one stream.
double pair_normal(ULong64_t key, Long64_t entry, h4l::Collection collection, Short_t i1, Short_t i2, std::uint64_t stream) {
  const std::uint64_t index = 1000 + 64 * static_cast<std::uint64_t>(i1) + static_cast<std::uint64_t>(i2);
  return h4l::object_normal(key, static_cast<std::uint64_t>(entry), collection, index, stream);
}

constexpr std::uint64_t kDataDeltaStream = 2, kTemplateStream = 3;
}  // namespace

int main(int argc, char** argv) {
  std::string job_path, out_path;
  for (int index = 1; index < argc; index += 2) {
    const std::string key = argv[index];
    if (index + 1 >= argc) {
      std::cerr << "missing value for " << key << "\n";
      return 64;
    }
    if (key == "--job") job_path = argv[index + 1];
    else if (key == "--out") out_path = argv[index + 1];
    else {
      std::cerr << "unknown argument " << key << "\n";
      return 64;
    }
  }
  if (job_path.empty() || out_path.empty()) {
    std::cerr << "usage: zpeak_histograms --job JOB.json --out OUT.root\n";
    return 64;
  }
  gROOT->SetBatch(true);
  SetErrorHandler(record_errors);
  try {
    const json job = h4l::read_json(job_path);
    const json& config = job.at("calibration_config");
    const std::string role = job.at("role").get<std::string>();
    if (role != "data" && role != "mc") throw std::runtime_error("role must be data or mc");
    const bool weighted = job.at("weighted").get<bool>();
    const double weight_scale = job.at("weight_scale").get<double>();
    const int half = job.value("half", -1);
    const std::uint64_t inject_stream = job.value("inject_stream", 1ULL);
    const bool delta_smear = job.value("delta_smear", false), write_template = job.value("write_template_pairs", false),
               write_full = job.value("write_full_pairs", false);
    if (write_template && role != "mc") throw std::runtime_error("template pairs are written for the MC role");
    const double delta = config.at("template_fit").at("common_delta").get<double>();
    const double min_pt_a = config.at("categories").at("family_a_min_pt").get<double>();
    const json& histogram_config = config.at("histogram");
    const double mass_low = histogram_config.at("mass_range").at(0).get<double>();
    const double mass_high = histogram_config.at("mass_range").at(1).get<double>();
    const int mass_bins = histogram_config.at("mass_bins").get<int>();
    const json& leg_pt_bins = histogram_config.at("leg_pt_bins");
    const json& pair_pt_bins = histogram_config.at("pair_pt_bins");
    const json payload = job.value("payload", json());
    const json inject = job.value("inject", json());

    std::map<int, Flavour> flavours;
    for (const auto& [name, pdg, tag] : {std::tuple{"muon", 13, "mm"}, std::tuple{"electron", 11, "ee"}}) {
      Flavour f;
      f.name = name;
      f.tag = tag;
      f.collection = pdg == 13 ? h4l::Collection::Muon : h4l::Collection::Electron;
      f.model = h4l::CalibModel::from_json(config.at("model").at(name));
      const json& selection = config.at("selection").at(name);
      f.min_pt = selection.at("min_pt").get<double>();
      f.max_abs_eta = selection.at("max_abs_eta").get<double>();
      f.max_iso = selection.at("max_iso_fsr").get<double>();
      f.max_sip = selection.at("max_sip").get<double>();
      f.id = selection.at("id").get<std::string>();
      if (inject.is_object() && inject.contains(name)) {
        f.inject = true;
        f.inject_scale = TestFunction::from_json(inject.at(name).at("scale"));
        f.inject_smear = TestFunction::from_json(inject.at(name).at("smear"));
      }
      if (payload.is_object() && payload.contains(name)) {
        f.payload = h4l::FactorizedPayload::from_json(payload.at(name));
        if (!(f.payload.model == f.model)) throw std::runtime_error("payload model differs from the configuration model");
      }
      flavours[pdg] = f;
    }

    // Reserve the output: exactly one run may produce it.
    const h4l::fs::path target(out_path);
    if (!target.is_absolute()) throw std::runtime_error("--out must be an absolute path");
    h4l::ensure_directory(target.parent_path());
    if (h4l::fs::exists(target)) {
      std::cerr << "ERROR: output exists: " << out_path << "\n";
      return 17;
    }
    {
      const std::string lock = out_path + ".lock";
      const int fd = ::open(lock.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0644);
      if (fd < 0) {
        std::cerr << "ERROR: output reserved by another run: " << lock << "\n";
        return 17;
      }
      const std::string note = "started " + h4l::utc_now() + " pid " + std::to_string(::getpid()) + "\n";
      const bool written = ::write(fd, note.data(), note.size()) == static_cast<ssize_t>(note.size());
      if (::close(fd) != 0 || !written) throw std::runtime_error("cannot write " + lock);
    }
    const std::string temporary = h4l::attempt_path(target, "hist").string() + ".tmp";
    std::unique_ptr<TFile> output(TFile::Open(temporary.c_str(), "CREATE"));
    if (!output || output->IsZombie()) throw std::runtime_error("cannot create " + temporary);

    // Trees live in the output file; histograms are held in memory.
    output->cd();
    Char_t t_flavour = 0;
    Short_t t_ai = -1, t_aj = -1, t_bi = -1, t_bj = -1;
    Float_t t_m = 0, t_w = 0, t_eps = 0;
    TTree* template_tree = nullptr;
    if (write_template) {
      template_tree = new TTree("TemplatePairs", "MC pairs of the template fits");
      template_tree->Branch("flavour", &t_flavour, "flavour/B");
      template_tree->Branch("a_i", &t_ai, "a_i/S");
      template_tree->Branch("a_j", &t_aj, "a_j/S");
      template_tree->Branch("b_i", &t_bi, "b_i/S");
      template_tree->Branch("b_j", &t_bj, "b_j/S");
      template_tree->Branch("m", &t_m, "m/F");
      template_tree->Branch("w", &t_w, "w/F");
      template_tree->Branch("eps", &t_eps, "eps/F");
    }
    Float_t f_pt[2], f_eta[2], f_phi[2], f_aeta[2], f_rel_err[2], f_r9[2], f_gen_pt[2], f_pair_pt = 0;
    Char_t f_charge[2];
    UChar_t f_gen_flav[2];
    TTree* full_tree = nullptr;
    if (write_full) {
      full_tree = new TTree("FullPairs", "selected transformed pairs");
      full_tree->Branch("flavour", &t_flavour, "flavour/B");
      full_tree->Branch("a_i", &t_ai, "a_i/S");
      full_tree->Branch("a_j", &t_aj, "a_j/S");
      full_tree->Branch("b_i", &t_bi, "b_i/S");
      full_tree->Branch("b_j", &t_bj, "b_j/S");
      full_tree->Branch("m", &t_m, "m/F");
      full_tree->Branch("w", &t_w, "w/F");
      full_tree->Branch("pair_pt", &f_pair_pt, "pair_pt/F");
      full_tree->Branch("pt", f_pt, "pt[2]/F");
      full_tree->Branch("eta", f_eta, "eta[2]/F");
      full_tree->Branch("phi", f_phi, "phi[2]/F");
      full_tree->Branch("abs_eta_bin", f_aeta, "abs_eta_bin[2]/F");
      full_tree->Branch("rel_err", f_rel_err, "rel_err[2]/F");
      full_tree->Branch("r9", f_r9, "r9[2]/F");
      full_tree->Branch("gen_pt", f_gen_pt, "gen_pt[2]/F");
      full_tree->Branch("charge", f_charge, "charge[2]/B");
      full_tree->Branch("gen_flav", f_gen_flav, "gen_flav[2]/b");
    }

    std::map<std::string, std::unique_ptr<TH1D>> hists;
    auto histogram = [&](const std::string& name, int bins, double low, double high) -> TH1D& {
      auto& h = hists[name];
      if (!h) h = make_h1(name, bins, low, high);
      return *h;
    };
    std::map<std::string, std::unique_ptr<TH1D>> stats;
    for (const auto& [pdg, f] : flavours) {
      const int nb = f.model.n_pt() * f.model.n_regions(), ne = f.model.n_eta();
      for (const char* stat : {"w", "wpt", "w2", "ws", "wr2"}) {
        const std::string name = std::string("leg_") + f.tag + "_" + stat;
        stats[name] = make_h1(name, nb, -0.5, nb - 0.5);
      }
      for (const char* stat : {"w", "weta", "ws", "wr2"}) {
        const std::string name = std::string("legeta_") + f.tag + "_" + stat;
        stats[name] = make_h1(name, ne, -0.5, ne - 0.5);
      }
    }
    auto counts = make_h1("counts", 8, -0.5, 7.5);  // 0 rows, 1 kept by half, 2/3 mm/ee selected (w), 4/5 mm/ee selected rows

    long long rows_total = 0, template_rows = 0, full_rows = 0;
    std::set<std::string> seen;
    for (const auto& input : job.at("inputs")) {
      const std::string path = input.at("path").get<std::string>();
      const long long expected_rows = input.at("calib_pairs").get<long long>();
      if (!seen.insert(path).second) throw std::runtime_error("duplicate input " + path);
      std::string open_error;
      auto file = h4l::open_input(path, 5, open_error);
      if (!file) throw std::runtime_error("cannot open " + path + ": " + open_error);
      g_read_error = false;
      auto* tree = dynamic_cast<TTree*>(file->Get("CalibPairs"));
      if (!tree) throw std::runtime_error("no CalibPairs in " + path);
      if (tree->GetEntries() != expected_rows)
        throw std::runtime_error(path + " has " + std::to_string(tree->GetEntries()) + " rows, the scan " +
                                 std::to_string(expected_rows));
      ULong64_t key = 0;
      Long64_t entry = 0;
      Float_t weight = 1, pair_mass = 0, pair_pt = 0;
      Char_t flavour = 0;
      LegIn legs[2];
      tree->SetBranchStatus("*", false);
      auto address = [&](const std::string& name, void* where) {
        tree->SetBranchStatus(name.c_str(), true);
        if (tree->SetBranchAddress(name.c_str(), where) < 0) throw std::runtime_error("cannot address " + name);
      };
      address("file_key", &key);
      address("entry", &entry);
      address("weight", &weight);
      address("flavour", &flavour);
      address("mass", &pair_mass);
      address("pair_pt", &pair_pt);
      for (int l = 0; l < 2; ++l) {
        const std::string p = l == 0 ? "l1_" : "l2_";
        address(p + "index", &legs[l].index);
        address(p + "pt", &legs[l].pt);
        address(p + "eta", &legs[l].eta);
        address(p + "eta_sc", &legs[l].eta_sc);
        address(p + "phi", &legs[l].phi);
        address(p + "pt_err", &legs[l].pt_err);
        address(p + "iso_fsr", &legs[l].iso_fsr);
        address(p + "sip", &legs[l].sip);
        address(p + "r9", &legs[l].r9);
        address(p + "gen_pt", &legs[l].gen_pt);
        address(p + "flags", &legs[l].flags);
        address(p + "charge", &legs[l].charge);
        address(p + "gen_flav", &legs[l].gen_flav);
      }
      output->cd();
      for (Long64_t row = 0; row < expected_rows; ++row) {
        if (tree->GetEntry(row) <= 0) throw std::runtime_error("cannot read row " + std::to_string(row) + " of " + path);
        if (g_read_error) throw std::runtime_error("ROOT error while reading " + path + ": " + g_read_message);
        ++rows_total;
        counts->Fill(0);
        if (half >= 0 && static_cast<int>(key % 2) != half) continue;
        counts->Fill(1);
        const auto found = flavours.find(flavour);
        if (found == flavours.end()) throw std::runtime_error("unknown flavour in " + path);
        const Flavour& f = found->second;
        double pt[2], abs_eta[2], injected_s[2] = {0, 0}, injected_r2[2] = {0, 0};
        bool pass = true;
        for (int l = 0; l < 2; ++l) {
          const LegIn& leg = legs[l];
          abs_eta[l] = std::fabs(f.model.eta_sc ? leg.eta_sc : leg.eta);
          pt[l] = leg.pt;
          if (f.inject) {
            const double normal = h4l::object_normal(key, static_cast<std::uint64_t>(entry), f.collection,
                                                     static_cast<std::uint64_t>(leg.index), inject_stream);
            injected_s[l] = f.inject_scale(pt[l], abs_eta[l]);
            const double r = std::max(0.0, f.inject_smear(pt[l], abs_eta[l]));
            injected_r2[l] = r * r;
            pt[l] *= (1.0 + injected_s[l]) * (1.0 + r * normal);
          }
          if (!f.payload.empty()) {
            if (role == "data") {
              pt[l] = f.payload.corrected_data_pt(pt[l], abs_eta[l]);
            } else {
              const double normal = h4l::object_normal(key, static_cast<std::uint64_t>(entry), f.collection,
                                                       static_cast<std::uint64_t>(leg.index), 0);
              pt[l] = f.payload.smeared_mc_pt(pt[l], abs_eta[l], normal);
            }
          }
          pass = pass && f.pass(pt[l], std::fabs(leg.eta), leg.iso_fsr, leg.sip, leg.flags);
        }
        if (!pass) continue;
        const int e1 = f.model.eta_bin(abs_eta[0]), e2 = f.model.eta_bin(abs_eta[1]);
        const int q1 = f.model.pt_region_bin(pt[0], abs_eta[0]), q2 = f.model.pt_region_bin(pt[1], abs_eta[1]);
        if (e1 < 0 || e2 < 0 || q1 < 0 || q2 < 0) throw std::runtime_error("selected leg outside the calibration bins");
        const double m = pair_mass * std::sqrt(pt[0] * pt[1] / (static_cast<double>(legs[0].pt) * legs[1].pt));
        const double w = (weighted ? weight : 1.0) * weight_scale;
        const bool muon = f.collection == h4l::Collection::Muon;
        counts->Fill(muon ? 2 : 3, w);
        counts->Fill(muon ? 4 : 5);
        const bool in_a = pt[0] >= min_pt_a && pt[1] >= min_pt_a;
        const bool peak = m > 80.0 && m < 100.0;
        const double m_delta = delta_smear ? m * (1.0 + delta * pair_normal(key, entry, f.collection, legs[0].index,
                                                                            legs[1].index, kDataDeltaStream))
                                           : 0.0;
        const int n_pt = f.model.n_pt(), n_eta = f.model.n_eta();
        for (const char family : {'A', 'B'}) {
          if (family == 'A' && !in_a) continue;
          const int b1 = family == 'A' ? e1 : q1, b2 = family == 'A' ? e2 : q2;
          const int i = std::min(b1, b2), j = std::max(b1, b2);
          const int first = b1 <= b2 ? 0 : 1;  // the leg in bin i
          const std::string suffix = std::string(1, family) + "_" + f.tag + "_" + std::to_string(i) + "_" + std::to_string(j);
          histogram("m" + suffix, mass_bins, mass_low, mass_high).Fill(m, w);
          if (delta_smear) histogram("md" + suffix, mass_bins, mass_low, mass_high).Fill(m_delta, w);
          TH1D& lpt = histogram("lpt" + suffix, leg_pt_bins.at(0).get<int>(), leg_pt_bins.at(1).get<double>(),
                                leg_pt_bins.at(2).get<double>());
          lpt.Fill(pt[0], w);
          lpt.Fill(pt[1], w);
          histogram("ppt" + suffix, pair_pt_bins.at(0).get<int>(), pair_pt_bins.at(1).get<double>(),
                    pair_pt_bins.at(2).get<double>()).Fill(pair_pt, w);
          if (peak) {
            // Composition: A records the pT bins, B the fine |eta| bins of its legs.
            const int size = family == 'A' ? n_pt : n_eta;
            TH1D& comp = histogram("w" + suffix, 2 * size, -0.5, 2 * size - 0.5);
            for (int l = 0; l < 2; ++l) {
              const int set = (i == j || l == first) ? 0 : 1;
              const int component = family == 'A' ? f.model.pt_bin(pt[l]) : (l == 0 ? e1 : e2);
              comp.Fill(set * size + component, w);
            }
          }
        }
        if (peak) {
          for (int l = 0; l < 2; ++l) {
            const int q = l == 0 ? q1 : q2, e = l == 0 ? e1 : e2;
            stats["leg_" + f.tag + "_w"]->Fill(q, w);
            stats["leg_" + f.tag + "_wpt"]->Fill(q, w * pt[l]);
            stats["leg_" + f.tag + "_w2"]->Fill(q, w * w);
            stats["leg_" + f.tag + "_ws"]->Fill(q, w * injected_s[l]);
            stats["leg_" + f.tag + "_wr2"]->Fill(q, w * injected_r2[l]);
            stats["legeta_" + f.tag + "_w"]->Fill(e, w);
            stats["legeta_" + f.tag + "_weta"]->Fill(e, w * abs_eta[l]);
            stats["legeta_" + f.tag + "_ws"]->Fill(e, w * injected_s[l]);
            stats["legeta_" + f.tag + "_wr2"]->Fill(e, w * injected_r2[l]);
          }
        }
        t_flavour = flavour;
        t_ai = in_a ? static_cast<Short_t>(std::min(e1, e2)) : -1;
        t_aj = in_a ? static_cast<Short_t>(std::max(e1, e2)) : -1;
        t_bi = static_cast<Short_t>(std::min(q1, q2));
        t_bj = static_cast<Short_t>(std::max(q1, q2));
        t_m = static_cast<Float_t>(m);
        t_w = static_cast<Float_t>(w);
        if (template_tree) {
          t_eps = static_cast<Float_t>(pair_normal(key, entry, f.collection, legs[0].index, legs[1].index, kTemplateStream));
          if (template_tree->Fill() <= 0) throw std::runtime_error("TemplatePairs Fill failed");
          ++template_rows;
        }
        if (full_tree) {
          f_pair_pt = pair_pt;
          for (int l = 0; l < 2; ++l) {
            f_pt[l] = static_cast<Float_t>(pt[l]);
            f_eta[l] = legs[l].eta;
            f_phi[l] = legs[l].phi;
            f_aeta[l] = static_cast<Float_t>(abs_eta[l]);
            f_rel_err[l] = legs[l].pt_err / legs[l].pt;
            f_r9[l] = legs[l].r9;
            f_gen_pt[l] = legs[l].gen_pt;
            f_charge[l] = legs[l].charge;
            f_gen_flav[l] = legs[l].gen_flav;
          }
          if (full_tree->Fill() <= 0) throw std::runtime_error("FullPairs Fill failed");
          ++full_rows;
        }
      }
      if (g_read_error) throw std::runtime_error("ROOT error while reading " + path + ": " + g_read_message);
    }

    // Write, close, reopen with retries, rename.
    std::vector<std::string> names;
    output->cd();
    auto write = [&](TObject& object) {
      if (object.Write() <= 0) throw std::runtime_error(std::string("cannot write ") + object.GetName());
      names.push_back(object.GetName());
    };
    for (auto& [name, h] : hists) write(*h);
    for (auto& [name, h] : stats) write(*h);
    write(*counts);
    TNamed job_text("job", h4l::dump_json(job).c_str());
    write(job_text);
    std::map<std::string, long long> trees;
    if (template_tree) {
      if (template_tree->Write() <= 0) throw std::runtime_error("cannot write TemplatePairs");
      trees["TemplatePairs"] = template_rows;
    }
    if (full_tree) {
      if (full_tree->Write() <= 0) throw std::runtime_error("cannot write FullPairs");
      trees["FullPairs"] = full_rows;
    }
    output->Close();
    output.reset();
    if (g_read_error) throw std::runtime_error("ROOT error while writing: " + g_read_message);
    std::string last_error;
    bool checked = false;
    for (int attempt = 1; attempt <= 6 && !checked; ++attempt) {
      try {
        h4l::validate_root_output(temporary, trees, names);
        std::unique_ptr<TFile> check(TFile::Open(temporary.c_str(), "READ"));
        auto* reread = check ? dynamic_cast<TH1D*>(check->Get("counts")) : nullptr;
        if (!reread || reread->GetBinContent(1) != static_cast<double>(rows_total))
          throw std::runtime_error("re-read counts differ in " + temporary);
        checked = true;
      } catch (const std::exception& error) {
        last_error = error.what();
        std::cerr << "WARN: check attempt " << attempt << " failed: " << last_error << "\n";
        std::this_thread::sleep_for(std::chrono::seconds(10 * attempt));
      }
    }
    if (!checked) throw std::runtime_error("the written output does not verify: " + last_error);
    if (h4l::fs::exists(target)) throw std::runtime_error("output appeared meanwhile: " + out_path + "; kept " + temporary);
    h4l::fs::rename(temporary, target);
    std::cout << "[histograms] " << out_path << ": " << rows_total << " rows, " << hists.size() << " histograms, "
              << template_rows << " template pairs, " << full_rows << " full pairs\n";
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
