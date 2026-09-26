// template_fit: event-level MC template fits of the Z -> ll mass in the
// leg-bin categories of the lepton calibration.
//
// Usage: template_fit --job JOB.json --out OUT.json
//
// The job names the histogram files of the data role (zpeak_histograms
// outputs with the delta-smeared md histograms) and of the MC role (with
// TemplatePairs trees), the category family ("A": fine |eta| pairs, "B":
// (pT, region) pairs), the flavour tag ("mm" or "ee") and the categories
// "<i>_<j>".  For every category the data histogram md (the data smeared
// once by the common known relative delta with frozen pair deviates) is
// fitted, in the fit binning, with the template of the category's MC pairs e
// moved to
//   mu_e = k m_e (1 + sqrt(delta^2 + D) eps_e),
// eps_e the frozen pair deviate, each spread by a Gaussian kernel of fixed
// relative width eta (T_b = sum_e w_e [Phi((x_{b+1} - mu_e) / (eta mu_e))
// - Phi((x_b - mu_e) / (eta mu_e))]).  The deviates carry the D-dependence,
// so the smoothing does not change with D and cannot favour a larger smear,
// while the fixed kernel keeps the likelihood smooth (without it, the MC
// noise gives wiggles and false local minima).  D is the residual extra
// relative pair variance of the data relative to the MC, down to -delta^2;
// the kernel adds the known variance eta^2, so E = D + eta^2 is reported as
// the residual pair variance.  The expected content is
// nu_b = N T_b / sum_window T with a free normalization N; the likelihood is
// Poisson, for a weighted data-role histogram (closure tests, response
// passes) a scaled Poisson in effective counts with ONE scale for the whole
// window, s = sum n / sum Var(n) over the window (a scale estimated per group
// from the observed data is biased and explodes in sparse groups: a group
// holding only a low-weight ZZ/ggZZ pair gets s = 1/w ~ 1e4 and dominates the
// likelihood), with one Barlow-Beeston-lite nuisance per likelihood group for
// the MC statistical uncertainty Var(T_g) = sum_e (w_e P_eg)^2, P_eg the
// share of event e in group g (the shares of one event in several bins of a
// group are fully correlated), profiled analytically.  MIGRAD and HESSE
// (Minuit2) give ln k, D and their covariance.
//
// The window is seeded from the data: the mode of the md histogram smoothed
// over the configured bins and searched in mode_search, then
// [mode - below, mode + above], clipped to window_clip and snapped to fit-bin
// edges.  The likelihood groups merge consecutive fit bins until they are at
// least min_group_width wide (several kernel widths, so that the kernel
// correlates neighbouring groups little and the per-group Barlow-Beeston
// nuisances carry the MC statistical uncertainty) and the MC events, placed
// at the starting values without the kernel, hold min_group_mc_entries
// effective entries (the kernel-smoothed contents would count each event
// about 2 sqrt(pi) eta m / fit_bin_width ~ 9 times).  Categories below the
// statistics thresholds are listed as skipped.
// Each result carries the curves (data, fitted template, template at k = 1
// and D = 0) for the fit galleries; fit quality is judged on those plots.

#include "h4l/hash.h"
#include "h4l/io.h"

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <TError.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <thread>
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

std::string hex64(std::uint64_t value) {
  char buffer[17];
  std::snprintf(buffer, sizeof(buffer), "%016llx", static_cast<unsigned long long>(value));
  return buffer;
}

struct Settings {
  double mode_low = 70, mode_high = 110, below = 20, above = 15, clip_low = 60, clip_high = 120;
  int smooth_bins = 45;
  double fit_bin_width = 0.1, delta = 0.005, kernel = 0.003, d_max = 2.5e-3, lnk_limit = 0.08;
  double min_data = 500, min_mc_effective = 500, min_group_mc = 20, min_group_width = 1.2;
  bool curves = true;
  int slices = 0;            // points per side of the diagnostic likelihood slices (0: none)
  double slice_sigmas = 5;   // half-width of the slices in fitted errors
};

struct Event {
  float m, w, eps;
};

std::unique_ptr<TFile> open_retrying(const std::string& path) {
  for (int attempt = 1; attempt <= 5; ++attempt) {
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
    if (file && !file->IsZombie() && !file->TestBit(TFile::kRecovered)) return file;
    std::this_thread::sleep_for(std::chrono::seconds(5 * attempt));
  }
  throw std::runtime_error("cannot open " + path);
}

std::unique_ptr<TH1D> summed(const std::vector<std::unique_ptr<TFile>>& files, const std::string& name) {
  std::unique_ptr<TH1D> total;
  for (const auto& file : files) {
    TObject* object = file->Get(name.c_str());
    if (!object) continue;
    auto* histogram = dynamic_cast<TH1D*>(object);
    if (!histogram) throw std::runtime_error(name + " is not a TH1D in " + file->GetName());
    if (!total) {
      total.reset(static_cast<TH1D*>(histogram->Clone((name + "_sum").c_str())));
      total->SetDirectory(nullptr);
    } else if (!total->Add(histogram)) {
      throw std::runtime_error("cannot add " + name);
    }
  }
  return total;
}

double smoothed_mode(const TH1D& h, double low, double high, int width) {
  const int n = h.GetNbinsX(), half = width / 2;
  std::vector<double> prefix(n + 1, 0.0);
  for (int b = 1; b <= n; ++b) prefix[b] = prefix[b - 1] + h.GetBinContent(b);
  int best = -1;
  double best_value = 0;
  for (int b = 1; b <= n; ++b) {
    const double x = h.GetBinCenter(b);
    if (x < low || x > high) continue;
    const int lo = std::max(1, b - half), hi = std::min(n, b + half);
    const double value = (prefix[hi] - prefix[lo - 1]) / (hi - lo + 1);
    if (best < 0 || value > best_value) {
      best = b;
      best_value = value;
    }
  }
  return best < 0 || !(best_value > 0) ? std::nan("") : h.GetBinCenter(best);
}

// Standard normal CDF, tabulated in [-6, 6] with linear interpolation (the
// template evaluates it tens of millions of times per likelihood call).
struct NormalCdf {
  static constexpr double kLimit = 6.0, kStep = 1e-3;
  std::vector<double> table;
  NormalCdf() {
    const int n = static_cast<int>(2 * kLimit / kStep) + 1;
    table.resize(n + 1);
    for (int i = 0; i <= n; ++i) table[i] = 0.5 * std::erfc(-(-kLimit + i * kStep) / std::sqrt(2.0));
  }
  double operator()(double z) const {
    if (z <= -kLimit) return 0.0;
    if (z >= kLimit) return 1.0;
    const double x = (z + kLimit) / kStep;
    const int i = static_cast<int>(x);
    const double f = x - i;
    return table[i] + f * (table[i + 1] - table[i]);
  }
};
const NormalCdf kCdf;

struct Category {
  double lo = 0, width = 0;
  int bins = 0;                     // fine fit bins of the template
  std::vector<int> group;           // likelihood group of each fine bin
  int groups = 0;
  std::vector<double> n, variance;  // data per likelihood group
  bool weighted = false;
  double data_scale = 1.0;          // effective counts per unit weight of the weighted data (window)
  const std::vector<Event>* events = nullptr;
  double delta2 = 0;

  double kernel = 0.003;

  // Template contents and MC variances per likelihood group at (k, D); the
  // shares of one event in the bins of a group are summed before squaring.
  void grouped(double k, double d, std::vector<double>& value, std::vector<double>& variance_out) const {
    value.assign(groups, 0.0);
    variance_out.assign(groups, 0.0);
    const double width_scale = std::sqrt(std::max(delta2 + d, 0.0));
    const double hi = lo + bins * width;
    for (const Event& e : *events) {
      const double mu = k * e.m * (1.0 + width_scale * e.eps);
      const double sigma = kernel * mu, reach = 4.5 * sigma;  // Phi(-4.5) = 3.4e-6
      if (mu + reach < lo || mu - reach > hi) continue;
      const int first = std::max(0, static_cast<int>(std::floor((mu - reach - lo) / width)));
      const int last = std::min(bins - 1, static_cast<int>(std::floor((mu + reach - lo) / width)));
      const double inv = 1.0 / sigma;
      double below = kCdf((lo + first * width - mu) * inv), share = 0;
      int g = group[first];
      for (int b = first; b <= last; ++b) {
        const double above = kCdf((lo + (b + 1) * width - mu) * inv);
        if (group[b] != g) {
          value[g] += share;
          variance_out[g] += share * share;
          g = group[b];
          share = 0;
        }
        share += (above - below) * e.w;
        below = above;
      }
      value[g] += share;
      variance_out[g] += share * share;
    }
  }

  double nll(const double* p, std::vector<double>* nu_out = nullptr) const {
    std::vector<double> value, variance_t;
    grouped(std::exp(p[0]), p[1], value, variance_t);
    double total = 0;
    for (double v : value) total += v;
    if (!(total > 0)) return 1e30;
    const double norm = p[2];
    double sum = 0;
    if (nu_out) nu_out->assign(groups, 0.0);
    for (int b = 0; b < groups; ++b) {
      double nu = norm * value[b] / total;
      const double data = n[b];
      const double scale = weighted ? data_scale : 1.0;
      if (!(nu > 0)) {
        if (data <= 0) continue;
        return 1e30;
      }
      const double rel2 = value[b] > 0 ? variance_t[b] / (value[b] * value[b]) : 0.0;
      double theta = 0;
      if (rel2 > 0) {
        const double a = 1.0 / rel2, nu_s = nu * scale, n_s = data * scale;
        const double bq = a + nu_s, cq = nu_s - n_s, disc = bq * bq - 4.0 * a * cq;
        theta = disc > 0 ? (-bq + std::sqrt(disc)) / (2.0 * a) : 0.0;
        theta = std::max(theta, -0.9);
        sum += 0.5 * a * theta * theta;
      }
      nu *= 1.0 + theta;
      if (nu_out) (*nu_out)[b] = nu;
      const double nu_s = nu * scale, n_s = data * scale;
      sum += nu_s - n_s;
      if (n_s > 0) sum += n_s * std::log(n_s / nu_s);
    }
    return sum;
  }
};

json fit_category(const TH1D& data_h, const std::vector<Event>& events, const Settings& s, json& entry) {
  const double mode = smoothed_mode(data_h, s.mode_low, s.mode_high, s.smooth_bins);
  entry["data_mode"] = mode;
  if (!std::isfinite(mode)) return json("no data peak");
  // Fit binning: an integer merge of the histogram bins.
  const double hist_width = data_h.GetXaxis()->GetBinWidth(1);
  const int merge = static_cast<int>(std::lround(s.fit_bin_width / hist_width));
  if (merge < 1 || std::fabs(merge * hist_width - s.fit_bin_width) > 1e-9)
    throw std::runtime_error("fit_bin_width is not a multiple of the histogram bin width");
  const double x0 = data_h.GetXaxis()->GetXmin();
  const double low = std::max(s.clip_low, mode - s.below), high = std::min(s.clip_high, mode + s.above);
  const int first_fit = static_cast<int>(std::ceil((low - x0) / s.fit_bin_width - 1e-9));
  const int last_fit = static_cast<int>(std::floor((high - x0) / s.fit_bin_width + 1e-9)) - 1;
  Category c;
  c.lo = x0 + first_fit * s.fit_bin_width;
  c.width = s.fit_bin_width;
  c.bins = last_fit - first_fit + 1;
  c.events = &events;
  c.delta2 = s.delta * s.delta;
  c.kernel = s.kernel;
  std::vector<double> fine_n, fine_variance;
  double n_data = 0, v_data = 0;
  for (int b = 0; b < c.bins; ++b) {
    double content = 0, variance = 0;
    for (int q = 0; q < merge; ++q) {
      const int bin = (first_fit + b) * merge + q + 1;
      content += data_h.GetBinContent(bin);
      variance += data_h.GetBinError(bin) * data_h.GetBinError(bin);
    }
    fine_n.push_back(content);
    fine_variance.push_back(variance);
    n_data += content;
    v_data += variance;
    if (std::fabs(content - std::round(content)) > 1e-6 || std::fabs(variance - content) > 1e-6 * std::max(1.0, content))
      c.weighted = true;
  }
  double mc_sum = 0, mc_w2 = 0;
  for (const Event& e : events)
    if (e.m >= c.lo && e.m < c.lo + c.bins * c.width) {
      mc_sum += e.w;
      mc_w2 += static_cast<double>(e.w) * e.w;
    }
  const double data_eff = v_data > 0 ? n_data * n_data / v_data : 0, mc_eff = mc_w2 > 0 ? mc_sum * mc_sum / mc_w2 : 0;
  entry["window"] = {c.lo, c.lo + c.bins * c.width};
  entry["data_entries"] = n_data;
  entry["data_effective_entries"] = data_eff;
  entry["data_weighted"] = c.weighted;
  entry["mc_events"] = events.size();
  entry["mc_effective_entries"] = mc_eff;
  if (data_eff < s.min_data || mc_eff < s.min_mc_effective) return json("low statistics");

  // Starting scale from the smeared-MC mode.
  TH1D mc_h("mc_mode", "", data_h.GetNbinsX(), data_h.GetXaxis()->GetXmin(), data_h.GetXaxis()->GetXmax());
  mc_h.SetDirectory(nullptr);
  for (const Event& e : events) mc_h.Fill(e.m * (1.0 + s.delta * e.eps), e.w);
  const double mc_mode = smoothed_mode(mc_h, s.mode_low, s.mode_high, s.smooth_bins);
  const double lnk0 = std::isfinite(mc_mode) ? std::clamp(std::log(mode / mc_mode), -0.05, 0.05) : 0.0;
  c.data_scale = n_data > 0 && v_data > 0 ? n_data / v_data : 1.0;
  // Likelihood groups: consecutive fine bins merged until the group is at
  // least min_group_width wide and the MC events at the starting values
  // (without the kernel) hold min_group_mc effective entries (fixed during
  // the fit; the template itself stays on the fine grid).
  {
    std::vector<double> t0(c.bins, 0.0), v0(c.bins, 0.0);
    const double k0 = std::exp(lnk0);
    for (const Event& e : events) {
      const int b = static_cast<int>(std::floor((k0 * e.m * (1.0 + s.delta * e.eps) - c.lo) / c.width));
      if (b < 0 || b >= c.bins) continue;
      t0[b] += e.w;
      v0[b] += static_cast<double>(e.w) * e.w;
    }
    c.group.assign(c.bins, 0);
    double sum_w = 0, sum_w2 = 0;
    int current = 0, start = 0;
    for (int b = 0; b < c.bins; ++b) {
      c.group[b] = current;
      sum_w += t0[b];
      sum_w2 += v0[b];
      const bool wide = (b - start + 1) * c.width >= s.min_group_width - 1e-9;
      if (wide && sum_w2 > 0 && sum_w * sum_w / sum_w2 >= s.min_group_mc && b + 1 < c.bins) {
        ++current;
        start = b + 1;
        sum_w = sum_w2 = 0;
      }
    }
    // A short last group joins its predecessor.
    if (current > 0 && !(sum_w2 > 0 && sum_w * sum_w / sum_w2 >= s.min_group_mc))
      for (int b = 0; b < c.bins; ++b)
        if (c.group[b] == current) c.group[b] = current - 1;
    c.groups = *std::max_element(c.group.begin(), c.group.end()) + 1;
    c.n.assign(c.groups, 0.0);
    c.variance.assign(c.groups, 0.0);
    for (int b = 0; b < c.bins; ++b) {
      c.n[c.group[b]] += fine_n[b];
      c.variance[c.group[b]] += fine_variance[b];
    }
  }
  entry["likelihood_groups"] = c.groups;
  auto function = [&](const double* p) { return c.nll(p); };
  json attempts = json::array();
  std::vector<double> best(3, 0.0);
  int status = -1, cov_status = -1, strategy_used = -1;
  double edm = -1, min_value = 0;
  std::vector<std::vector<double>> cov(3, std::vector<double>(3, 0.0));
  const double d_low = -0.98 * c.delta2;
  for (int strategy : {1, 2}) {
    std::unique_ptr<ROOT::Math::Minimizer> minimizer(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    if (!minimizer) throw std::runtime_error("cannot create the Minuit2 minimizer");
    ROOT::Math::Functor functor(function, 3);
    minimizer->SetFunction(functor);
    minimizer->SetErrorDef(0.5);
    minimizer->SetStrategy(strategy);
    minimizer->SetMaxFunctionCalls(20000);
    minimizer->SetTolerance(0.01);
    minimizer->SetPrintLevel(-1);
    minimizer->SetLimitedVariable(0, "lnk", lnk0, 2e-4, -s.lnk_limit, s.lnk_limit);
    minimizer->SetLimitedVariable(1, "D", 0.0, 0.2 * c.delta2, d_low, s.d_max);
    minimizer->SetLimitedVariable(2, "norm", n_data, 0.01 * n_data + 1.0, 0.0, 3.0 * n_data + 10.0);
    minimizer->Minimize();
    minimizer->Hesse();
    status = minimizer->Status();
    cov_status = minimizer->CovMatrixStatus();
    strategy_used = strategy;
    edm = minimizer->Edm();
    min_value = minimizer->MinValue();
    best.assign(minimizer->X(), minimizer->X() + 3);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) cov[i][j] = minimizer->CovMatrix(i, j);
    attempts.push_back({{"strategy", strategy}, {"status", status}, {"cov_status", cov_status}, {"edm", edm}});
    if (status == 0 && cov_status >= 2) break;
  }
  json fit;
  fit["status"] = status;
  fit["cov_status"] = cov_status;
  fit["strategy"] = strategy_used;
  fit["edm"] = edm;
  fit["nll"] = min_value;
  fit["attempts"] = attempts;
  const double e_lnk = std::sqrt(std::max(cov[0][0], 0.0)), e_d = std::sqrt(std::max(cov[1][1], 0.0));
  const double kernel_variance = s.kernel * s.kernel;
  fit["lnk"] = {{"value", best[0]}, {"error", e_lnk}};
  fit["D"] = {{"value", best[1]}, {"error", e_d}};
  fit["E"] = {{"value", best[1] + kernel_variance}, {"error", e_d}, {"kernel_variance", kernel_variance}};
  fit["norm"] = {{"value", best[2]}, {"error", std::sqrt(std::max(cov[2][2], 0.0))}};
  fit["corr_lnk_E"] = e_lnk > 0 && e_d > 0 ? cov[0][1] / (e_lnk * e_d) : 0.0;
  const bool at_limit = std::fabs(std::fabs(best[0]) - s.lnk_limit) < 1e-6 || best[1] <= d_low * 0.999 || best[1] >= s.d_max * 0.999;
  fit["at_limit"] = at_limit;
  fit["ok"] = status == 0 && cov_status >= 2 && !at_limit;
  std::vector<double> nu, nu0;
  c.nll(best.data(), &nu);
  const double start[3] = {0.0, 0.0, best[2]};
  c.nll(start, &nu0);
  // Pearson chi2 in effective counts (the observed variance of a sparse
  // weighted group is not its expected variance).
  double chi2 = 0;
  int used = 0;
  const double chi2_scale = c.weighted ? c.data_scale : 1.0;
  for (int b = 0; b < c.groups; ++b)
    if (nu[b] > 0) {
      chi2 += chi2_scale * (c.n[b] - nu[b]) * (c.n[b] - nu[b]) / nu[b];
      ++used;
    }
  fit["chi2"] = chi2;
  fit["ndf"] = used - 3;
  if (s.slices > 0) {
    // Likelihood slices through the minimum in ln k and in D (diagnostics).
    json slices = json::object();
    for (int par : {0, 1}) {
      const double error = par == 0 ? e_lnk : e_d;
      const double half = s.slice_sigmas * (error > 0 ? error : (par == 0 ? 1e-4 : c.delta2));
      json xs = json::array(), ys = json::array();
      for (int q = -s.slices; q <= s.slices; ++q) {
        double p[3] = {best[0], best[1], best[2]};
        p[par] += half * q / s.slices;
        if (par == 1 && p[1] < d_low) continue;
        xs.push_back(p[par]);
        ys.push_back(c.nll(p) - min_value);
      }
      slices[par == 0 ? "lnk" : "D"] = {{"x", xs}, {"dnll", ys}};
    }
    fit["slices"] = slices;
  }
  if (s.curves) {
    std::vector<double> x_lo(c.groups, 1e30), x_hi(c.groups, -1e30), ey(c.groups);
    for (int b = 0; b < c.bins; ++b) {
      x_lo[c.group[b]] = std::min(x_lo[c.group[b]], c.lo + b * c.width);
      x_hi[c.group[b]] = std::max(x_hi[c.group[b]], c.lo + (b + 1) * c.width);
    }
    for (int b = 0; b < c.groups; ++b) ey[b] = std::sqrt(c.variance[b]);
    fit["curve"] = {{"x_lo", x_lo}, {"x_hi", x_hi}, {"y", c.n}, {"ey", ey}, {"model", nu}, {"template_before", nu0}};
  }
  return fit;
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
    std::cerr << "usage: template_fit --job JOB.json --out OUT.json\n";
    return 64;
  }
  gROOT->SetBatch(true);
  TH1::AddDirectory(false);
  try {
    const json job = h4l::read_json(job_path);
    const json& fit = job.at("calibration_config").at("template_fit");
    const json variant = job.value("variant", json::object());
    Settings s;
    s.mode_low = fit.at("mode_search").at(0).get<double>();
    s.mode_high = fit.at("mode_search").at(1).get<double>();
    s.smooth_bins = fit.at("smooth_bins").get<int>();
    s.below = variant.value("below", fit.at("below").get<double>());
    s.above = variant.value("above", fit.at("above").get<double>());
    s.clip_low = fit.at("window_clip").at(0).get<double>();
    s.clip_high = fit.at("window_clip").at(1).get<double>();
    s.fit_bin_width = variant.value("fit_bin_width", fit.at("fit_bin_width").get<double>());
    s.delta = fit.at("common_delta").get<double>();
    s.kernel = fit.at("kernel_rel_sigma").get<double>();
    s.d_max = fit.at("d_max").get<double>();
    s.lnk_limit = fit.at("lnk_limit").get<double>();
    s.min_data = fit.at("min_data_events").get<double>();
    s.min_mc_effective = fit.at("min_mc_effective").get<double>();
    s.min_group_mc = fit.at("min_group_mc_entries").get<double>();
    s.min_group_width = fit.at("min_group_width").get<double>();
    s.curves = job.value("curves", true);
    s.slices = job.value("slices", 0);
    s.slice_sigmas = job.value("slice_sigmas", 5.0);
    const std::string flavour = job.at("flavour").get<std::string>();
    if (flavour != "mm" && flavour != "ee") throw std::runtime_error("flavour must be mm or ee");
    const int pdg = flavour == "mm" ? 13 : 11;
    const std::string family = job.at("family").get<std::string>();
    if (family != "A" && family != "B") throw std::runtime_error("family must be A or B");
    const std::string column = family == "A" ? "a" : "b";

    SetErrorHandler(record_errors);
    std::vector<std::unique_ptr<TFile>> data_files;
    std::set<std::string> seen;
    for (const auto& path_node : job.at("data")) {
      const std::string path = path_node.get<std::string>();
      if (!seen.insert(path).second) throw std::runtime_error("duplicate input " + path);
      data_files.push_back(open_retrying(path));
    }
    // MC template pairs of this flavour, bucketed by category.
    std::map<std::string, std::vector<Event>> buckets;
    long long mc_rows = 0;
    for (const auto& path_node : job.at("mc")) {
      const std::string path = path_node.get<std::string>();
      if (!seen.insert(path).second) throw std::runtime_error("duplicate input " + path);
      auto file = open_retrying(path);
      g_read_error = false;
      auto* tree = dynamic_cast<TTree*>(file->Get("TemplatePairs"));
      if (!tree) throw std::runtime_error("no TemplatePairs in " + path);
      Char_t t_flavour = 0;
      Short_t t_i = 0, t_j = 0;
      Float_t t_m = 0, t_w = 0, t_eps = 0;
      const std::string name_i = column + "_i", name_j = column + "_j";
      for (auto [name, where] : std::initializer_list<std::pair<const char*, void*>>{
               {"flavour", &t_flavour}, {name_i.c_str(), &t_i}, {name_j.c_str(), &t_j}, {"m", &t_m}, {"w", &t_w},
               {"eps", &t_eps}})
        if (tree->SetBranchAddress(name, where) < 0) throw std::runtime_error(std::string("cannot address ") + name);
      const Long64_t n = tree->GetEntries();
      for (Long64_t row = 0; row < n; ++row) {
        if (tree->GetEntry(row) <= 0) throw std::runtime_error("cannot read TemplatePairs row of " + path);
        if (t_flavour != pdg || t_i < 0) continue;
        buckets[std::to_string(t_i) + "_" + std::to_string(t_j)].push_back({t_m, t_w, t_eps});
        ++mc_rows;
      }
      if (g_read_error) throw std::runtime_error("ROOT error while reading " + path + ": " + g_read_message);
    }
    json results = json::array();
    for (const auto& category_node : job.at("categories")) {
      const std::string category = category_node.get<std::string>();
      g_read_error = false;
      auto data = summed(data_files, "md" + family + "_" + flavour + "_" + category);
      if (g_read_error) throw std::runtime_error("ROOT error while reading " + category + ": " + g_read_message);
      json entry = {{"category", category}, {"ok", false}};
      const auto bucket = buckets.find(category);
      if (!data || bucket == buckets.end()) {
        entry["skipped"] = !data ? "no data histogram" : "no MC template pairs";
        results.push_back(entry);
        continue;
      }
      const json result = fit_category(*data, bucket->second, s, entry);
      if (result.is_string()) {
        entry["skipped"] = result;
      } else {
        entry["fit"] = result;
        entry["ok"] = result.at("ok");
        std::cout << "[template] " << family << "_" << flavour << "_" << category << ": ln k " << result["lnk"]["value"].get<double>()
                  << " +- " << result["lnk"]["error"].get<double>() << ", E " << result["E"]["value"].get<double>() << " +- "
                  << result["E"]["error"].get<double>() << (entry["ok"].get<bool>() ? "" : "  [NOT OK]") << "\n";
      }
      results.push_back(entry);
    }
    json report = {{"schema", "h4l_v3_template_fits/3"},
                   {"flavour", flavour},
                   {"family", family},
                   {"settings", fit},
                   {"variant", variant},
                   {"data", job.at("data")},
                   {"mc", job.at("mc")},
                   {"mc_template_rows", mc_rows},
                   {"program", job.value("program", json())},
                   {"job_fnv1a64", hex64(h4l::fnv1a64(job.dump()))},
                   {"results", results},
                   {"finished_utc", h4l::utc_now()}};
    h4l::publish_json(out_path, report);
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
