// Lepton calibration helpers shared by the calibration programs and every
// consumer of the calibration payload.
//
// Factorized per-lepton model of one flavour (user decision 2026-09-24):
//   u(pT, |eta|) = ln(1 + s) = a[e] + b_R(pT),
//   v(pT, |eta|) = r^2       = c[e] + d_R(pT),
// e the fine |eta| bin (piecewise constant), R the coarse |eta| region, and
// b_R, d_R linear in pT between the region's node positions (the mean pT of
// the Z legs of each pT bin in the region), constant beyond the outermost
// nodes, zero at the reference pT bin.  |eta| is abs(eta) for muons and
// abs(eta_SC) for electrons.  Bins are given by lower edges; internal
// boundaries are half-open [low, high) and the terminal bin of each axis is
// unbounded.  A value below the first edge has no bin (-1).
//
// Application convention (per lepton, with the pT the analysis reconstructs,
// before any correction):
//   data: pT -> pT / (1 + s) = pT exp(-u);
//   MC:   pT -> pT (1 + r N), r = sqrt(max(v, 0)), N the per-lepton deviate
//         of h4l/hash.h (object_normal with stream 0).
#pragma once

#include "h4l/io.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace h4l {

inline int lower_edge_bin(const std::vector<double>& lower, double value) {
  if (lower.empty() || !(value >= lower.front())) return -1;
  int bin = 0;
  while (bin + 1 < static_cast<int>(lower.size()) && value >= lower[bin + 1]) ++bin;
  return bin;
}

inline std::vector<double> increasing_edges(const json& node, const char* key) {
  std::vector<double> edges;
  for (const auto& value : node.at(key)) edges.push_back(value.get<double>());
  if (edges.empty()) throw std::runtime_error(std::string("empty ") + key);
  for (std::size_t k = 0; k + 1 < edges.size(); ++k)
    if (!(edges[k] < edges[k + 1])) throw std::runtime_error(std::string(key) + " are not increasing");
  return edges;
}

// The binning of the factorized model of one flavour.
struct CalibModel {
  std::vector<double> eta_lower, region_lower, pt_lower;
  bool eta_sc = false;
  int reference_pt_bin = 0;

  int n_eta() const { return static_cast<int>(eta_lower.size()); }
  int n_regions() const { return static_cast<int>(region_lower.size()); }
  int n_pt() const { return static_cast<int>(pt_lower.size()); }
  int eta_bin(double abs_eta) const { return lower_edge_bin(eta_lower, abs_eta); }
  int region(double abs_eta) const { return lower_edge_bin(region_lower, abs_eta); }
  int pt_bin(double pt) const { return lower_edge_bin(pt_lower, pt); }
  // The (pT bin, region) bin of family-B categories: k = ipt * n_regions + region.
  int pt_region_bin(double pt, double abs_eta) const {
    const int ip = pt_bin(pt), ir = region(abs_eta);
    return ip < 0 || ir < 0 ? -1 : ip * n_regions() + ir;
  }
  bool operator==(const CalibModel& other) const {
    return eta_lower == other.eta_lower && region_lower == other.region_lower && pt_lower == other.pt_lower &&
           eta_sc == other.eta_sc && reference_pt_bin == other.reference_pt_bin;
  }
  static CalibModel from_json(const json& node) {
    CalibModel model;
    model.eta_lower = increasing_edges(node, "eta_edges");
    model.region_lower = increasing_edges(node, "region_edges");
    model.pt_lower = increasing_edges(node, "pt_edges");
    const std::string variable = node.at("eta_variable").get<std::string>();
    if (variable != "abs_eta" && variable != "abs_eta_sc") throw std::runtime_error("unknown eta_variable " + variable);
    model.eta_sc = variable == "abs_eta_sc";
    model.reference_pt_bin = node.at("reference_pt_bin").get<int>();
    if (model.reference_pt_bin < 0 || model.reference_pt_bin >= model.n_pt())
      throw std::runtime_error("reference_pt_bin outside the pT bins");
    if (model.eta_lower.front() != 0.0 || model.region_lower.front() != 0.0)
      throw std::runtime_error("the eta and region edges must start at 0");
    return model;
  }
};

// Named lepton identification working points on the stage-2 flag bits.
// Muon bits: 0 isGlobal, 1 isTracker, 3 isPFcand, 10 nStations>0.
// Electron bits: 0/1/2 mvaFall17V2noIso WPL/WP90/WP80, 3/4/5 Iso WPL/WP90/WP80.
inline bool muon_id(const std::string& name, unsigned flags) {
  const bool global = flags & 1U, tracker = flags >> 1 & 1U, pf = flags >> 3 & 1U, stations = flags >> 10 & 1U;
  if (name == "an_tight") return pf && (global || (tracker && stations));
  if (name == "an_loose") return global || (tracker && stations);
  throw std::runtime_error("unknown muon id " + name);
}
// The AN-16-442 tight muon (3.2.1): a loose muon that is a PF muon or, above
// high_pt GeV, passes the tracker high-pT ID (stage-2 bit 8, highPtId >= 1).
inline bool an_tight_muon(unsigned flags, double pt, double high_pt = 200.0) {
  const bool loose = (flags & 1U) || ((flags >> 1 & 1U) && (flags >> 10 & 1U));
  return loose && ((flags >> 3 & 1U) || (pt > high_pt && (flags >> 8 & 1U)));
}
// Whether an FSR photon at dr from a lepton lies inside the lepton's isolation
// cone and outside its veto (AN-16-442 3.3 step 6): muons dr > 0.01;
// electrons |eta_SC| < 1.479 or dr > 0.08.
inline bool fsr_in_isolation(bool muon, double abs_eta_sc, double dr, double cone = 0.3) {
  if (dr >= cone) return false;
  return muon ? dr > 0.01 : (abs_eta_sc < 1.479 || dr > 0.08);
}
inline bool electron_id(const std::string& name, unsigned flags) {
  if (name == "mvaFall17V2noIso_WPL") return flags & 1U;
  if (name == "mvaFall17V2noIso_WP90") return flags >> 1 & 1U;
  if (name == "mvaFall17V2noIso_WP80") return flags >> 2 & 1U;
  if (name == "mvaFall17V2Iso_WPL") return flags >> 3 & 1U;
  if (name == "mvaFall17V2Iso_WP90") return flags >> 4 & 1U;
  if (name == "mvaFall17V2Iso_WP80") return flags >> 5 & 1U;
  if (name == "none") return true;
  throw std::runtime_error("unknown electron id " + name);
}

// The factorized payload of one flavour.
struct FactorizedPayload {
  CalibModel model;
  std::vector<double> a, c;                        // per fine |eta| bin
  std::vector<std::vector<double>> node_pt, b, d;  // [region][pT bin]
  bool empty() const { return a.empty(); }

  static double along_pt(const std::vector<double>& nodes, const std::vector<double>& values, double pt) {
    const int n = static_cast<int>(nodes.size());
    if (n == 1 || pt <= nodes.front()) return values.front();
    if (pt >= nodes.back()) return values.back();
    int i = 0;
    while (i + 2 < n && pt >= nodes[i + 1]) ++i;
    const double t = (pt - nodes[i]) / (nodes[i + 1] - nodes[i]);
    return (1 - t) * values[i] + t * values[i + 1];
  }
  double u(double pt, double abs_eta) const {
    const int e = model.eta_bin(abs_eta), r = model.region(abs_eta);
    if (e < 0 || r < 0) throw std::runtime_error("lepton outside the calibration bins");
    return a[e] + along_pt(node_pt[r], b[r], pt);
  }
  double v(double pt, double abs_eta) const {
    const int e = model.eta_bin(abs_eta), r = model.region(abs_eta);
    if (e < 0 || r < 0) throw std::runtime_error("lepton outside the calibration bins");
    return c[e] + along_pt(node_pt[r], d[r], pt);
  }
  double scale_at(double pt, double abs_eta) const { return std::expm1(u(pt, abs_eta)); }
  double smear_at(double pt, double abs_eta) const { return std::sqrt(std::max(v(pt, abs_eta), 0.0)); }
  double corrected_data_pt(double pt, double abs_eta) const { return pt * std::exp(-u(pt, abs_eta)); }
  double smeared_mc_pt(double pt, double abs_eta, double normal) const {
    return pt * (1.0 + smear_at(pt, abs_eta) * normal);
  }
  static FactorizedPayload from_json(const json& node) {
    FactorizedPayload p;
    p.model = CalibModel::from_json(node.at("model"));
    auto vector = [](const json& values) {
      std::vector<double> out;
      for (const auto& value : values) {
        const double x = value.get<double>();
        if (!std::isfinite(x)) throw std::runtime_error("non-finite payload value");
        out.push_back(x);
      }
      return out;
    };
    p.a = vector(node.at("a"));
    p.c = vector(node.at("c"));
    for (const auto& row : node.at("node_pt")) p.node_pt.push_back(vector(row));
    for (const auto& row : node.at("b")) p.b.push_back(vector(row));
    for (const auto& row : node.at("d")) p.d.push_back(vector(row));
    const std::size_t n_eta = p.model.n_eta(), n_regions = p.model.n_regions(), n_pt = p.model.n_pt();
    if (p.a.size() != n_eta || p.c.size() != n_eta || p.node_pt.size() != n_regions || p.b.size() != n_regions ||
        p.d.size() != n_regions)
      throw std::runtime_error("payload sizes differ from its model");
    for (std::size_t r = 0; r < n_regions; ++r) {
      if (p.node_pt[r].size() != n_pt || p.b[r].size() != n_pt || p.d[r].size() != n_pt)
        throw std::runtime_error("payload pT sizes differ from its model");
      for (std::size_t k = 0; k + 1 < n_pt; ++k)
        if (!(p.node_pt[r][k] < p.node_pt[r][k + 1])) throw std::runtime_error("payload pT nodes are not increasing");
      if (p.b[r][p.model.reference_pt_bin] != 0.0 || p.d[r][p.model.reference_pt_bin] != 0.0)
        throw std::runtime_error("payload pT terms are not zero at the reference bin");
    }
    for (double x : p.a)
      if (!(std::fabs(x) < 0.5)) throw std::runtime_error("undefined scale in the payload");
    return p;
  }
};
}  // namespace h4l
