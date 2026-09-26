// lambda_histograms: Z -> ll mass histograms in bins of the predicted
// per-event mass uncertainty, for the lambda calibration of AN-16-442 5.3.
//
// Usage: lambda_histograms --job JOB.json --out OUT.root
//
// The job names FullPairs files (zpeak_histograms outputs of the final
// calibration payload: data corrected, MC smeared), the lambda configuration
// (lambda_ul16_v2.json, embedded as "lambda_config") and the weighting.
// Every leg has a relative momentum error d = dpT/pT (Muon_ptErr / pT;
// Electron_energyErr / E for electrons) and a lambda region: the first
// configured region whose |eta| (|eta_SC| for electrons) and d ranges contain
// it (ranges half-open, a null upper edge unbounded).  A pair with legs in
// regions a <= b and predicted relative mass uncertainty
//   e = 0.5 sqrt(d1^2 + d2^2)
// fills, with k the bin of e (internal edges from the configuration, the
// first bin from 0 and the last one unbounded),
//   lm_<f>_<a>_<b>_<k>   the pair mass;
//   le_<f>_<a>_<b>_<k>   for pairs in 80 < m < 100 GeV the sums of w, w e,
//                        w e^2, w d_a^2, w d_b^2 (entries 1-5), d_a the
//                        error of the leg in region a (for a = b leg 1).
// The output is written to an attempt-qualified temporary, checked and
// renamed; an existing output is never overwritten (exit code 17).

#include "h4l/io.h"
#include "h4l/root_io.h"

#include <TError.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TTree.h>

#include <cmath>
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

bool g_read_error = false;
std::string g_read_message;
void record_errors(int level, Bool_t abort, const char* location, const char* message) {
  if (level >= kError) {
    g_read_error = true;
    g_read_message = std::string(location ? location : "?") + ": " + (message ? message : "");
  }
  DefaultErrorHandler(level, abort, location, message);
}

double upper_edge(const json& range) {
  return range.at(1).is_null() ? std::numeric_limits<double>::infinity() : range.at(1).get<double>();
}

struct Region {
  double eta_low = 0, eta_high = 0, err_low = 0, err_high = 0;
  bool contains(double abs_eta, double rel_err) const {
    return abs_eta >= eta_low && abs_eta < eta_high && rel_err >= err_low && rel_err < err_high;
  }
};

struct Flavour {
  std::string tag;
  std::vector<Region> regions;
  std::vector<double> e_edges;  // internal edges
  int region(double abs_eta, double rel_err) const {
    for (std::size_t k = 0; k < regions.size(); ++k)
      if (regions[k].contains(abs_eta, rel_err)) return static_cast<int>(k);
    return -1;
  }
  int e_bin(double e) const {
    int k = 0;
    while (k < static_cast<int>(e_edges.size()) && e >= e_edges[k]) ++k;
    return k;
  }
};

Flavour read_flavour(const json& node, const std::string& tag) {
  Flavour f;
  f.tag = tag;
  for (const auto& r : node.at("regions")) {
    Region region;
    region.eta_low = r.at("abs_eta").at(0).get<double>();
    region.eta_high = upper_edge(r.at("abs_eta"));
    if (r.contains("rel_err")) {
      region.err_low = r.at("rel_err").at(0).get<double>();
      region.err_high = upper_edge(r.at("rel_err"));
    } else {
      region.err_low = 0;
      region.err_high = std::numeric_limits<double>::infinity();
    }
    if (!(region.eta_high > region.eta_low) || !(region.err_high > region.err_low))
      throw std::runtime_error("empty lambda region for " + tag);
    f.regions.push_back(region);
  }
  for (const auto& value : node.at("pair_error_edges")) f.e_edges.push_back(value.get<double>());
  for (std::size_t k = 0; k < f.e_edges.size(); ++k)
    if (!(f.e_edges[k] > 0) || (k > 0 && !(f.e_edges[k] > f.e_edges[k - 1])))
      throw std::runtime_error("pair_error_edges must be positive and increasing for " + tag);
  return f;
}
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
    std::cerr << "usage: lambda_histograms --job JOB.json --out OUT.root\n";
    return 64;
  }
  gROOT->SetBatch(true);
  TH1::AddDirectory(false);
  SetErrorHandler(record_errors);
  try {
    const json job = h4l::read_json(job_path);
    const json& lambda = job.at("lambda_config");
    const bool weighted = job.at("weighted").get<bool>();
    const double mass_low = lambda.at("mass_range").at(0).get<double>(), mass_high = lambda.at("mass_range").at(1).get<double>();
    const int mass_bins = lambda.at("mass_bins").get<int>();
    const std::map<int, Flavour> flavours = {{13, read_flavour(lambda.at("muon"), "mm")},
                                             {11, read_flavour(lambda.at("electron"), "ee")}};
    const h4l::fs::path target(out_path);
    if (!target.is_absolute()) throw std::runtime_error("--out must be an absolute path");
    if (h4l::fs::exists(target)) {
      std::cerr << "ERROR: output exists: " << out_path << "\n";
      return 17;
    }
    h4l::ensure_directory(target.parent_path());
    std::map<std::string, std::unique_ptr<TH1D>> hists;
    auto histogram = [&](const std::string& name, int bins, double low, double high) -> TH1D& {
      auto& h = hists[name];
      if (!h) {
        h = std::make_unique<TH1D>(name.c_str(), name.c_str(), bins, low, high);
        h->Sumw2();
      }
      return *h;
    };
    long long rows = 0, used = 0;
    std::set<std::string> seen;
    for (const auto& input : job.at("inputs")) {
      const std::string path = input.at("path").get<std::string>();
      const long long expected = input.at("full_pairs").get<long long>();
      if (!seen.insert(path).second) throw std::runtime_error("duplicate input " + path);
      std::string open_error;
      auto file = h4l::open_input(path, 5, open_error);
      if (!file) throw std::runtime_error("cannot open " + path + ": " + open_error);
      g_read_error = false;
      auto* tree = dynamic_cast<TTree*>(file->Get("FullPairs"));
      if (!tree) throw std::runtime_error("no FullPairs in " + path);
      if (tree->GetEntries() != expected) throw std::runtime_error(path + ": FullPairs entries differ from the job");
      Char_t flavour = 0;
      Float_t m = 0, w = 0, aeta[2], rel_err[2];
      for (auto [name, where] : std::initializer_list<std::pair<const char*, void*>>{
               {"flavour", &flavour}, {"m", &m}, {"w", &w}, {"abs_eta_bin", aeta}, {"rel_err", rel_err}})
        if (tree->SetBranchAddress(name, where) < 0) throw std::runtime_error(std::string("cannot address ") + name);
      for (Long64_t row = 0; row < expected; ++row) {
        if (tree->GetEntry(row) <= 0) throw std::runtime_error("cannot read FullPairs row of " + path);
        if (g_read_error) throw std::runtime_error("ROOT error while reading " + path + ": " + g_read_message);
        ++rows;
        const auto found = flavours.find(flavour);
        if (found == flavours.end()) throw std::runtime_error("unknown flavour in " + path);
        const Flavour& f = found->second;
        int r1 = f.region(aeta[0], rel_err[0]), r2 = f.region(aeta[1], rel_err[1]);
        if (r1 < 0 || r2 < 0) continue;
        double d1 = rel_err[0], d2 = rel_err[1];
        if (r1 > r2) {
          std::swap(r1, r2);
          std::swap(d1, d2);
        }
        const double weight = weighted ? w : 1.0;
        const double e = 0.5 * std::hypot(d1, d2);
        const std::string suffix = f.tag + "_" + std::to_string(r1) + "_" + std::to_string(r2) + "_" + std::to_string(f.e_bin(e));
        histogram("lm_" + suffix, mass_bins, mass_low, mass_high).Fill(m, weight);
        if (m > 80 && m < 100) {
          TH1D& stats = histogram("le_" + suffix, 5, 0.5, 5.5);
          stats.Fill(1, weight);
          stats.Fill(2, weight * e);
          stats.Fill(3, weight * e * e);
          stats.Fill(4, weight * d1 * d1);
          stats.Fill(5, weight * d2 * d2);
        }
        ++used;
      }
      if (g_read_error) throw std::runtime_error("ROOT error while reading " + path + ": " + g_read_message);
    }
    const std::string temporary = h4l::attempt_path(target, "lambda").string() + ".tmp";
    std::vector<std::string> names;
    {
      std::unique_ptr<TFile> output(TFile::Open(temporary.c_str(), "CREATE"));
      if (!output || output->IsZombie()) throw std::runtime_error("cannot create " + temporary);
      output->cd();
      for (auto& [name, h] : hists) {
        if (h->Write() <= 0) throw std::runtime_error("cannot write " + name);
        names.push_back(name);
      }
      TNamed job_text("job", h4l::dump_json(job).c_str());
      if (job_text.Write() <= 0) throw std::runtime_error("cannot write the job");
      names.push_back("job");
      output->Close();
    }
    if (g_read_error) throw std::runtime_error("ROOT error while writing: " + g_read_message);
    h4l::validate_root_output(temporary, {}, names);
    if (h4l::fs::exists(target)) throw std::runtime_error("output appeared meanwhile: " + out_path);
    h4l::fs::rename(temporary, target);
    std::cout << "[lambda] " << out_path << ": " << rows << " pairs, " << used << " in lambda regions, " << hists.size()
              << " histograms\n";
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
