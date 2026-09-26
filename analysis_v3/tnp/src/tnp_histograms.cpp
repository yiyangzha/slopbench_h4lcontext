// tnp_histograms: pass and fail mass histograms of the Z tag-and-probe.
//
// Usage: tnp_histograms --job JOB.json --out OUT.root
//
// The job names TnPPairs files of one sample (tnp_extract outputs, as
// {path, tnp_pairs, root_sha256} from the extraction scan), the
// tag-and-probe configuration (tnp_ul16_v3.json), the calibration payload
// (h4l/calibration.h), the bit of every trigger path in the skim trigger word,
// the role ("data" or "mc") and the weight scale.  For every pair:
//   1. both legs are calibrated: data pT *= exp(-u), MC pT *= (1 + r N) with
//      the per-lepton deviate of stream 0 (as in every stage); the pair mass
//      is rescaled by sqrt(pT1' pT2' / (pT1 pT2));
//   2. each leg that is a tag candidate of the extraction is tried as the tag,
//      with the final path-consistent tag selection on the calibrated pT;
//   3. the other leg is a probe if it passes the AN loose selection and lies
//      at dR > min_dr_tag_probe from the tag (the pairs are opposite-sign);
//   4. the selection chain of the analysis, each step conditional on the
//      previous one: id | loose, sip | id, iso | id and sip, and the direct
//      full = id, sip, iso | loose.  The isolation is the AN FSR-subtracted
//      one (AN 3.3; h4l::fsr_in_isolation): the FSR photons of the pair's
//      muon legs that pass the loose selection and SIP are removed from the
//      neutral isolation of a leg that passes the loose selection and SIP;
//      the tight muon is the AN one (h4l::an_tight_muon);
//   5. per tag and per step, only the probe whose pair mass is closest to
//      m_Z is kept (the step's own denominator applied before the choice);
//      the rows of one event are consecutive in TnPPairs;
//   6. the chosen probe fills pass_<f>_<step>_<b> or fail_<f>_<step>_<b> for
//      its bin (calibrated pT, |eta| or |eta_SC|; lower edges, the last bin
//      unbounded) and for <b> = "all".
// MC probes also fill truth_<f>_<step>_<b> (entries 1-2: fail/pass weights of
// probes matched to a prompt generator lepton, gen_flav 1, with
// 80 < m < 100 GeV) as a diagnostic of the fit method only.  Jobs with
// fill_templates (the DY MC) fill the signal templates of the fits,
// tpass_/tfail_<f>_<step>_<b>, from the pairs whose tag and probe are both
// matched to prompt generator leptons, in the wider template mass range
// (histogram.template_range, so that the Gaussian smearing of the fit has no
// edge effect), with the same per-(tag, step) choice of the probe closest to
// m_Z among those pairs; <b> is the bin, "pt<i>" (the pT row of the bin, all
// |eta|) and "all" (fallbacks of sparse bins).  Per bin,
// pstat_<f>_<b> holds the sums of w, w pT, w |eta| of the loose probes chosen
// for the full step in 80 < m < 100 GeV (the starting widths of the fits).
// Closure mode (MC jobs with a "closure" block, the efficiency closure of
// run_tnp.py --closure): only the pairs whose original file key has the given
// parity are used, and a probe that passes the identification loses it with
// probability 1 - keep_id (per flavour; per lepton from the uniform of stream 4
// of its object seed, so the same lepton is thinned in every pair): an
// efficiency loss of the id step, and so of the full selection.
// The output is written to an attempt-qualified temporary file, checked and
// renamed; an existing output is never overwritten (exit code 17).

#include "h4l/calibration.h"
#include "h4l/hash.h"
#include "h4l/io.h"
#include "h4l/kinematics.h"
#include "h4l/root_io.h"

#include <TError.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TTree.h>

#include <fcntl.h>
#include <unistd.h>

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
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

constexpr int kSteps = 4;
const char* const kStepNames[kSteps] = {"id", "sip", "iso", "full"};

struct Leg {
  Short_t index = 0;
  Float_t pt = 0, eta = 0, eta_sc = 0, phi = 0, iso03 = 0, iso03_chg = 0, iso04 = 0, sip = 0, trig_pt = -1, fsr_pt = -1,
          fsr_eta = 0, fsr_phi = 0;
  UShort_t flags = 0;
  Int_t trig_bits = 0;
  UChar_t cut_based = 0, gen_flav = 0;
};

struct TagPath {
  unsigned bit = 0;
  int filter = 0;
  double trig_pt_min = 0, max_abs_eta_sc = 99;
};

struct Flavour {
  std::string name, tag;
  h4l::Collection collection = h4l::Collection::Muon;
  h4l::FactorizedPayload payload;
  bool eta_sc = false;
  std::vector<TagPath> tag_paths;
  double tag_min_pt = 0, tag_max_eta = 0, tag_max_iso04 = 99, gap_low = 9, gap_high = 9;
  int tag_min_cut_based = 0;
  double probe_min_pt = 0, probe_max_eta = 0, min_dr = 0.2;
  std::string probe_id, id;
  double max_iso = 0, max_sip = 0, muon_high_pt = 200;
  std::vector<double> pt_lower, eta_lower;

  bool muon() const { return collection == h4l::Collection::Muon; }
  bool tag_ok(const Leg& leg, double pt, UInt_t trig) const {
    const double abs_eta = std::fabs(leg.eta), abs_eta_sc = std::fabs(leg.eta_sc);
    bool fired = false;
    for (const auto& path : tag_paths)
      if ((trig >> path.bit & 1U) && (leg.trig_bits & path.filter) && leg.trig_pt >= path.trig_pt_min &&
          abs_eta_sc < path.max_abs_eta_sc)
        fired = true;
    if (!fired || pt <= tag_min_pt) return false;
    if (muon()) return abs_eta < tag_max_eta && (leg.flags >> 6 & 1U) && leg.iso04 < tag_max_iso04;
    return abs_eta_sc < tag_max_eta && !(abs_eta_sc >= gap_low && abs_eta_sc < gap_high) && leg.cut_based >= tag_min_cut_based;
  }
  bool loose(const Leg& leg, double pt) const {
    if (pt <= probe_min_pt || std::fabs(leg.eta) >= probe_max_eta) return false;
    return muon() ? h4l::muon_id(probe_id, leg.flags) : h4l::electron_id(probe_id, leg.flags);
  }
  bool pass_id(const Leg& leg, double pt) const {
    if (muon()) {
      if (id == "an_tight_hzz") return h4l::an_tight_muon(leg.flags, pt, muon_high_pt);
      return h4l::muon_id(id, leg.flags);
    }
    return h4l::electron_id(id, leg.flags);
  }
  int bin(double pt, double abs_eta) const {
    const int ip = h4l::lower_edge_bin(pt_lower, pt), ie = h4l::lower_edge_bin(eta_lower, abs_eta);
    return ip < 0 || ie < 0 ? -1 : ip * static_cast<int>(eta_lower.size()) + ie;
  }
  int pt_row(double pt) const { return h4l::lower_edge_bin(pt_lower, pt); }
};

std::vector<double> vector_of(const json& node) {
  std::vector<double> out;
  for (const auto& value : node) out.push_back(value.get<double>());
  return out;
}

// One tag-probe combination of an event.
struct Combination {
  int pdg = 0, tag_index = 0, bin = -1, pt_row = -1;
  double m = 0, w = 0, pt = 0, abs_eta = 0;
  bool in[kSteps] = {}, pass[kSteps] = {}, truth = false, truth_pair = false, in_window = false;
};
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
    std::cerr << "usage: tnp_histograms --job JOB.json --out OUT.root\n";
    return 64;
  }
  gROOT->SetBatch(true);
  TH1::AddDirectory(false);
  SetErrorHandler(record_errors);
  try {
    const json job = h4l::read_json(job_path);
    const json& cfg = job.at("tnp_config");
    const std::string role = job.at("role").get<std::string>();
    if (role != "data" && role != "mc") throw std::runtime_error("role must be data or mc");
    const bool is_mc = role == "mc";
    const double weight_scale = job.at("weight_scale").get<double>();
    const json& bits = job.at("trigger_bits");
    const double mass_low = cfg.at("histogram").at("mass_range").at(0).get<double>();
    const double mass_high = cfg.at("histogram").at("mass_range").at(1).get<double>();
    const int mass_bins = cfg.at("histogram").at("mass_bins").get<int>();
    const double z_mass = cfg.at("fit").at("z_mass").get<double>();
    const double template_low = cfg.at("histogram").at("template_range").at(0).get<double>();
    const double template_high = cfg.at("histogram").at("template_range").at(1).get<double>();
    const int template_bins = cfg.at("histogram").at("template_bins").get<int>();
    const bool fill_templates = job.value("fill_templates", false);
    if (fill_templates && !is_mc) throw std::runtime_error("templates are filled from MC only");
    int closure_parity = -1;
    std::map<int, double> keep_id;  // pdg -> probability to keep a passed identification
    if (job.contains("closure")) {
      if (!is_mc) throw std::runtime_error("the closure mode runs on MC only");
      const json& closure = job.at("closure");
      closure_parity = closure.at("parity").get<int>();
      if (closure_parity != 0 && closure_parity != 1) throw std::runtime_error("the closure parity must be 0 or 1");
      keep_id[13] = closure.at("keep_id").at("muon").get<double>();
      keep_id[11] = closure.at("keep_id").at("electron").get<double>();
      for (const auto& [pdg, k] : keep_id) {
        if (!(k > 0.0 && k <= 1.0)) throw std::runtime_error("keep_id must lie in (0, 1]");
        if (fill_templates && k < 1.0) throw std::runtime_error("templates come from the unthinned role only");
      }
    }
    constexpr std::uint64_t kThinningStream = 4;
    if (template_low > mass_low || template_high < mass_high) throw std::runtime_error("the template range must contain the mass range");

    std::map<int, Flavour> flavours;
    for (const auto& [name, pdg, tag] : {std::tuple{"muon", 13, "mm"}, std::tuple{"electron", 11, "ee"}}) {
      Flavour f;
      f.name = name;
      f.tag = tag;
      f.collection = pdg == 13 ? h4l::Collection::Muon : h4l::Collection::Electron;
      f.payload = h4l::FactorizedPayload::from_json(job.at("payload").at(name));
      f.eta_sc = f.payload.model.eta_sc;
      const json& t = cfg.at("tag").at(name);
      for (const auto& path : t.at("paths")) {
        TagPath p;
        p.bit = bits.at(path.at(0).get<std::string>()).get<unsigned>();
        p.filter = path.at(1).get<int>();
        p.trig_pt_min = path.at(2).get<double>();
        if (path.size() > 3) p.max_abs_eta_sc = path.at(3).get<double>();
        f.tag_paths.push_back(p);
      }
      f.tag_min_pt = t.at("min_pt").get<double>();
      if (f.muon()) {
        f.tag_max_eta = t.at("max_abs_eta").get<double>();
        f.tag_max_iso04 = t.at("max_iso04").get<double>();
      } else {
        f.tag_max_eta = t.at("max_abs_eta_sc").get<double>();
        f.gap_low = t.at("gap").at(0).get<double>();
        f.gap_high = t.at("gap").at(1).get<double>();
        f.tag_min_cut_based = t.at("min_cut_based").get<int>();
      }
      const json& probe = cfg.at("probe").at(name);
      f.probe_min_pt = probe.at("min_pt").get<double>();
      f.probe_max_eta = probe.at("max_abs_eta").get<double>();
      f.probe_id = probe.at("id").get<std::string>();
      f.min_dr = cfg.at("probe").at("min_dr_tag_probe").get<double>();
      const json& steps = cfg.at("steps").at(name);
      f.id = steps.at("id").get<std::string>();
      f.max_iso = steps.at("max_iso_fsr").get<double>();
      f.max_sip = steps.at("max_sip").get<double>();
      if (f.muon()) f.muon_high_pt = steps.at("high_pt_alternative_above").get<double>();
      f.pt_lower = vector_of(cfg.at("bins").at(name).at("pt_edges"));
      f.eta_lower = vector_of(cfg.at("bins").at(name).at("eta_edges"));
      flavours[pdg] = f;
    }

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

    std::map<std::string, std::unique_ptr<TH1D>> hists;
    auto histogram = [&](const std::string& name, int n, double low, double high) -> TH1D& {
      auto& h = hists[name];
      if (!h) {
        h = std::make_unique<TH1D>(name.c_str(), name.c_str(), n, low, high);
        h->Sumw2();
      }
      return *h;
    };
    long long rows_total = 0, combinations_total = 0, chosen_total = 0;
    // Fill the chosen probe of every (tag, step) of one event.
    std::vector<Combination> event;
    auto flush = [&]() {
      for (int s = 0; s < kSteps; ++s) {
        // Signal templates: among the prompt-prompt pairs in the template range, the probe closest to m_Z
        // per tag.
        if (fill_templates) {
          std::map<std::pair<int, int>, int> best_truth;
          for (int i = 0; i < static_cast<int>(event.size()); ++i) {
            const Combination& c = event[i];
            if (!c.in[s] || !c.truth_pair) continue;
            const auto key = std::make_pair(c.pdg, c.tag_index);
            const auto found = best_truth.find(key);
            if (found == best_truth.end() || std::fabs(c.m - z_mass) < std::fabs(event[found->second].m - z_mass))
              best_truth[key] = i;
          }
          for (const auto& [key, i] : best_truth) {
            (void)key;
            const Combination& c = event[i];
            const Flavour& f = flavours.at(c.pdg);
            for (const std::string& label : {std::to_string(c.bin), "pt" + std::to_string(c.pt_row), std::string("all")}) {
              const std::string name = (c.pass[s] ? "tpass_" : "tfail_") + f.tag + "_" + kStepNames[s] + "_" + label;
              histogram(name, template_bins, template_low, template_high).Fill(c.m, c.w);
            }
          }
        }
        std::map<std::pair<int, int>, int> best;  // (pdg, tag index) -> combination
        for (int i = 0; i < static_cast<int>(event.size()); ++i) {
          const Combination& c = event[i];
          if (!c.in[s] || !c.in_window) continue;
          const auto key = std::make_pair(c.pdg, c.tag_index);
          const auto found = best.find(key);
          if (found == best.end() || std::fabs(c.m - z_mass) < std::fabs(event[found->second].m - z_mass)) best[key] = i;
        }
        for (const auto& [key, i] : best) {
          (void)key;
          const Combination& c = event[i];
          const Flavour& f = flavours.at(c.pdg);
          const bool peak = c.m > 80.0 && c.m < 100.0;
          for (const std::string& label : {std::to_string(c.bin), std::string("all")}) {
            const std::string suffix = f.tag + "_" + kStepNames[s] + "_" + label;
            histogram((c.pass[s] ? "pass_" : "fail_") + suffix, mass_bins, mass_low, mass_high).Fill(c.m, c.w);
            if (is_mc && peak && c.truth) histogram("truth_" + suffix, 2, 0.5, 2.5).Fill(c.pass[s] ? 2 : 1, c.w);
            if (s == 3 && peak) {
              TH1D& stat = histogram("pstat_" + f.tag + "_" + label, 3, 0.5, 3.5);
              stat.Fill(1, c.w);
              stat.Fill(2, c.w * c.pt);
              stat.Fill(3, c.w * c.abs_eta);
            }
          }
          ++chosen_total;
        }
      }
      event.clear();
    };

    std::set<std::string> seen;
    for (const auto& input : job.at("inputs")) {
      const std::string path = input.at("path").get<std::string>();
      const long long expected = input.at("tnp_pairs").get<long long>();
      if (!seen.insert(path).second) throw std::runtime_error("duplicate input " + path);
      std::string open_error;
      auto file = h4l::open_input(path, 5, open_error);
      if (!file) throw std::runtime_error("cannot open " + path + ": " + open_error);
      g_read_error = false;
      auto* tree = dynamic_cast<TTree*>(file->Get("TnPPairs"));
      if (!tree) throw std::runtime_error("no TnPPairs in " + path);
      if (tree->GetEntries() != expected)
        throw std::runtime_error(path + " has " + std::to_string(tree->GetEntries()) + " rows, the scan " + std::to_string(expected));
      ULong64_t key = 0, current_key = 0;
      Long64_t entry = 0, current_entry = -1;
      Float_t weight = 1, mass = 0;
      Char_t flavour = 0;
      UInt_t trig = 0;
      UChar_t candidates = 0;
      Leg legs[2];
      tree->SetBranchStatus("*", false);
      auto address = [&](const std::string& name, void* where) {
        tree->SetBranchStatus(name.c_str(), true);
        if (tree->SetBranchAddress(name.c_str(), where) < 0) throw std::runtime_error("cannot address " + name);
      };
      address("file_key", &key);
      address("entry", &entry);
      address("weight", &weight);
      address("flavour", &flavour);
      address("mass", &mass);
      address("trig", &trig);
      address("tag_candidates", &candidates);
      for (int l = 0; l < 2; ++l) {
        const std::string p = l == 0 ? "l1_" : "l2_";
        address(p + "index", &legs[l].index);
        for (auto [field, where] : std::initializer_list<std::pair<const char*, Float_t*>>{
                 {"pt", &legs[l].pt}, {"eta", &legs[l].eta}, {"eta_sc", &legs[l].eta_sc}, {"phi", &legs[l].phi},
                 {"iso03", &legs[l].iso03}, {"iso03_chg", &legs[l].iso03_chg}, {"iso04", &legs[l].iso04},
                 {"sip", &legs[l].sip}, {"trig_pt", &legs[l].trig_pt}, {"fsr_pt", &legs[l].fsr_pt},
                 {"fsr_eta", &legs[l].fsr_eta}, {"fsr_phi", &legs[l].fsr_phi}})
          address(p + field, where);
        address(p + "flags", &legs[l].flags);
        address(p + "trig_bits", &legs[l].trig_bits);
        address(p + "cut_based", &legs[l].cut_based);
        address(p + "gen_flav", &legs[l].gen_flav);
      }
      std::set<std::pair<ULong64_t, Long64_t>> events_seen;
      for (Long64_t row = 0; row < expected; ++row) {
        if (tree->GetEntry(row) <= 0) throw std::runtime_error("cannot read row " + std::to_string(row) + " of " + path);
        if (g_read_error) throw std::runtime_error("ROOT error while reading " + path + ": " + g_read_message);
        ++rows_total;
        if (key != current_key || entry != current_entry) {
          flush();
          if (!events_seen.insert({key, entry}).second)
            throw std::runtime_error("the pairs of one event are not consecutive in " + path);
          current_key = key;
          current_entry = entry;
        }
        if (closure_parity >= 0 && static_cast<int>(key % 2) != closure_parity) continue;
        const auto found = flavours.find(flavour);
        if (found == flavours.end()) throw std::runtime_error("unknown flavour in " + path);
        const Flavour& f = found->second;
        double pt[2], abs_eta_bin[2];
        for (int l = 0; l < 2; ++l) {
          abs_eta_bin[l] = std::fabs(f.eta_sc ? legs[l].eta_sc : legs[l].eta);
          if (is_mc) {
            const double normal = h4l::object_normal(key, static_cast<std::uint64_t>(entry), f.collection,
                                                     static_cast<std::uint64_t>(legs[l].index), 0);
            pt[l] = f.payload.smeared_mc_pt(legs[l].pt, abs_eta_bin[l], normal);
          } else {
            pt[l] = f.payload.corrected_data_pt(legs[l].pt, abs_eta_bin[l]);
          }
        }
        const double m = mass * std::sqrt(pt[0] * pt[1] / (static_cast<double>(legs[0].pt) * legs[1].pt));
        const bool in_window = m > mass_low && m < mass_high;
        // Outside the fit histograms a pair only feeds the templates (prompt-prompt pairs of the DY MC).
        const bool truth_pair_event = fill_templates && legs[0].gen_flav == 1 && legs[1].gen_flav == 1;
        if (!in_window && !(truth_pair_event && m > template_low && m < template_high)) continue;
        const double w = (is_mc ? weight : 1.0) * weight_scale;
        // AN FSR-subtracted isolation of each leg (defined for loose legs passing SIP).
        bool loose_sip[2];
        for (int l = 0; l < 2; ++l) loose_sip[l] = f.loose(legs[l], pt[l]) && legs[l].sip < f.max_sip;
        double iso[2];
        for (int l = 0; l < 2; ++l) {
          iso[l] = legs[l].iso03;
          if (!loose_sip[l]) continue;
          double photons = 0;
          for (int o = 0; o < 2; ++o) {
            if (!f.muon() || !loose_sip[o] || !(legs[o].fsr_pt > 0)) continue;
            const double dr = h4l::delta_r(legs[l].eta, legs[l].phi, legs[o].fsr_eta, legs[o].fsr_phi);
            if (h4l::fsr_in_isolation(f.muon(), std::fabs(legs[l].eta_sc), dr)) photons += legs[o].fsr_pt;
          }
          iso[l] = legs[l].iso03_chg + std::max(0.0, static_cast<double>(legs[l].iso03 - legs[l].iso03_chg) - photons / legs[l].pt);
        }
        for (int t = 0; t < 2; ++t) {
          if (!(candidates >> t & 1U)) continue;
          if (!f.tag_ok(legs[t], pt[t], trig)) continue;
          const int p = 1 - t;
          const Leg& probe = legs[p];
          if (!f.loose(probe, pt[p])) continue;
          if (h4l::delta_r(legs[t].eta, legs[t].phi, probe.eta, probe.phi) <= f.min_dr) continue;
          const int b = f.bin(pt[p], abs_eta_bin[p]);
          if (b < 0) continue;
          bool id = f.pass_id(probe, pt[p]);
          if (id && closure_parity >= 0) {
            const double k = keep_id.at(flavour);
            if (k < 1.0 && h4l::object_uniform(key, static_cast<std::uint64_t>(entry), f.collection,
                                               static_cast<std::uint64_t>(probe.index), kThinningStream) >= k)
              id = false;
          }
          const bool sip = probe.sip < f.max_sip, iso_ok = iso[p] < f.max_iso;
          Combination c;
          c.pdg = flavour;
          c.tag_index = legs[t].index;
          c.bin = b;
          c.m = m;
          c.w = w;
          c.pt = pt[p];
          c.abs_eta = abs_eta_bin[p];
          c.truth = is_mc && probe.gen_flav == 1;
          c.truth_pair = truth_pair_event && m > template_low && m < template_high;
          c.in_window = in_window;
          c.pt_row = f.pt_row(pt[p]);
          // The chain id | loose, sip | id, iso | id and sip; full | loose.
          c.in[0] = true;
          c.pass[0] = id;
          c.in[1] = id;
          c.pass[1] = sip;
          c.in[2] = id && sip;
          c.pass[2] = iso_ok;
          c.in[3] = true;
          c.pass[3] = id && sip && iso_ok;
          event.push_back(c);
          if (in_window) ++combinations_total;
        }
      }
      flush();
      if (g_read_error) throw std::runtime_error("ROOT error while reading " + path + ": " + g_read_message);
    }
    auto counts = std::make_unique<TH1D>("counts", "counts", 3, -0.5, 2.5);
    counts->SetBinContent(1, static_cast<double>(rows_total));
    counts->SetBinContent(2, static_cast<double>(combinations_total));
    counts->SetBinContent(3, static_cast<double>(chosen_total));
    const std::string temporary = h4l::attempt_path(target, "tnp").string() + ".tmp";
    std::vector<std::string> names;
    {
      std::unique_ptr<TFile> output(TFile::Open(temporary.c_str(), "CREATE"));
      if (!output || output->IsZombie()) throw std::runtime_error("cannot create " + temporary);
      output->cd();
      for (auto& [name, h] : hists) {
        if (h->Write() <= 0) throw std::runtime_error("cannot write " + name);
        names.push_back(name);
      }
      if (counts->Write() <= 0) throw std::runtime_error("cannot write counts");
      names.push_back("counts");
      TNamed job_text("job", h4l::dump_json(job).c_str());
      if (job_text.Write() <= 0) throw std::runtime_error("cannot write the job");
      names.push_back("job");
      output->Close();
    }
    if (g_read_error) throw std::runtime_error("ROOT error while writing: " + g_read_message);
    h4l::validate_root_output(temporary, {}, names);
    if (h4l::fs::exists(target)) throw std::runtime_error("output appeared meanwhile: " + out_path);
    h4l::fs::rename(temporary, target);
    std::cout << "[tnp] " << out_path << ": " << rows_total << " pairs, " << combinations_total
              << " tag-probe combinations, " << chosen_total << " chosen (tag, step) probes, " << hists.size()
              << " histograms\n";
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
