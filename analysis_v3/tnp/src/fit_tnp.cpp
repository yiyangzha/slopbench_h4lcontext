// fit_tnp: simultaneous pass/fail fits of the Z tag-and-probe.
//
// Usage: fit_tnp --job JOB.json --out OUT.json
//        fit_tnp --task TASK.json --out-json OUT.json   (the same, as a Condor task of run_task.sh)
//
// The job names the histogram files of the data and of the MC
// (tnp_histograms outputs, each role the sum of its files), the flavour tag,
// the step, the probe bins ("all" = inclusive) and the number of |eta| bins
// per pT row.  The binning and the inspection practice follow
// efficiency_tnp_an; the line shapes follow the user decision of 2026-09-24:
//   * pass and fail spectra of a bin (the configured window, 60-120 GeV by
//     default, the narrowest of 0.5, 1, ..., 3 GeV bins with at least 40
//     entries per bin of the total spectrum) are fitted simultaneously with
//     the efficiency as a parameter: N_pass = eps N_sig, N_fail = (1 - eps)
//     N_sig, independent background yields and shapes per category;
//   * nominal signal: the template of the prompt-prompt tag-probe pairs of
//     the DY MC (tpass_/tfail_ histograms of the same bin and category;
//     below template.min_effective_entries the bin's |eta| column with the
//     neighbouring pT rows added nearest first, then the bin's pT row, then
//     the inclusive one; smoothed with an adaptive Gaussian kernel whose
//     width follows the local template density, template.smoothing)
//     convolved with a Gaussian (shift and width per category; one Gaussian
//     for both when a category has fewer than shared_resolution_below
//     entries).  The fail spectra contain genuine signal with degraded mass
//     (FSR photons inside the isolation cone, mismeasured probes), which the
//     template carries;
//   * nominal background: CMSShape, erfc((alpha - m) beta) exp(-gamma (m -
//     m_Z)), per category, parameter ranges of the configuration (the EGM
//     tag-and-probe ranges: alpha 50-80 GeV, beta 0.01-0.06 / GeV, gamma
//     0.005-1 / GeV); a category with fewer than small_category_below
//     entries keeps the background shape of the inclusive fit (its yield
//     floats);
//   * alternatives (fit-model systematic): alt_signal = an analytic signal, a
//     stand-alone double-sided Crystal Ball with a shift and width per
//     category and shared tails (starting widths from the width model
//     sigma_incl W(bin) / W(incl), W = (1 + a <|eta|>)(1 + b <pT> / 100
//     GeV)), the fail one mixed with a Gaussian low-mass component (mean and
//     width ranges of the configuration, free fraction), CMSShape background;
//     alt_background = the template signal with a cubic Bernstein
//     background;
//   * two starting points (the background shapes of the inclusive fit and
//     default shapes), the converged fit with the lower NLL kept;
//     convergence: Minuit status 0 or 1 (covariance forced positive), or
//     status 3 with EDM < max_edm;
//   * the efficiency uncertainty is the MINOS (profile-likelihood) interval
//     (HESSE fails when parameters sit at their bounds; for weighted MC the
//     interval and HESSE are scaled by the effective-weight factor
//     sqrt(sum w^2 / sum w) of the spectra: the numerical sandwich was
//     unstable with weakly constrained shape nuisances), else HESSE, else a
//     MINOS interval that Minuit flags (a side at eps = 1) but that brackets
//     the minimum (error_source records which); an uncertainty below
//     the binomial floor sqrt(eps (1 - eps) / N_sig,eff) is impossible for a
//     fit: it is raised to the floor and flagged;
//   * background shape parameters at a bound count as "at limit" only when
//     that category's background exceeds 2 % of its entries;
//   * a fit is rejected when a background has a narrow maximum near m_Z
//     (within 12 GeV, falling below 70 % of it within 8 GeV on both sides: a
//     background patching the signal); broad combinatorial maxima are
//     physical;
//   * per-bin overrides (the "overrides" list of the configuration) record the
//     manual fixes found by inspecting the fit galleries: window, starting or
//     fixed parameter values, minimum rebinned width;
//   * a category with fewer than min_category_entries raw entries is reported
//     with its counting efficiency only.
// Weighted MC gets the effective-weight covariance (h4l/binned_fit.h).  Each result
// has the efficiency and its error per model, the counting efficiency in
// 80-100 GeV, the fit status, parameters at non-physical limits, the template
// sources and the curves of every model (bin errors = sqrt of the variances); for
// the MC also the generator-truth efficiency of prompt probes, a diagnostic
// that never enters a scale factor.

#include "h4l/binned_fit.h"
#include "h4l/io.h"
#include "h4l/zpeak_model.h"

#include <TError.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <thread>
#include <vector>

namespace {
using h4l::json;

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
    if (!histogram) throw std::runtime_error(name + " is not a TH1D");
    if (!total) {
      total.reset(static_cast<TH1D*>(histogram->Clone((name + "_sum").c_str())));
      total->SetDirectory(nullptr);
    } else if (!total->Add(histogram)) {
      throw std::runtime_error("cannot add " + name);
    }
  }
  return total;
}

struct Spectrum {
  double lo = 60, width = 1;
  int bins = 0;
  h4l::FitHistogram pass, fail;
  double raw_pass = 0, raw_fail = 0, sum_w = 0, sum_w2 = 0;
};

// Rebin to the narrowest width in {0.5, 1, 1.5, ..., 3} GeV (not below min_width) with >= 40 entries per bin of the total.
Spectrum make_spectrum(const TH1D* pass, const TH1D* fail, double low, double high, double min_width) {
  Spectrum s;
  const TH1D* ref = pass ? pass : fail;
  const double storage = ref->GetXaxis()->GetBinWidth(1);
  double total = 0;
  for (const TH1D* h : {pass, fail})
    if (h)
      for (int b = 1; b <= h->GetNbinsX(); ++b)
        if (h->GetBinCenter(b) >= low && h->GetBinCenter(b) < high) total += h->GetBinContent(b);
  double width = 1.0;
  for (double candidate : {0.5, 1.0, 1.5, 2.0, 2.5, 3.0}) {
    width = candidate;
    if (candidate + 1e-9 < min_width) continue;
    if (total / ((high - low) / candidate) >= 40) break;
  }
  const int merge = static_cast<int>(std::lround(width / storage));
  s.lo = low;
  s.width = merge * storage;
  s.bins = static_cast<int>(std::lround((high - low) / s.width));
  for (auto [h, target] : {std::pair{pass, &s.pass}, std::pair{fail, &s.fail}}) {
    target->content.assign(s.bins, 0.0);
    target->variance.assign(s.bins, 0.0);
    if (!h) continue;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      const double x = h->GetBinCenter(b);
      if (x < low || x >= high) continue;
      const int k = std::min(s.bins - 1, static_cast<int>((x - low) / s.width));
      target->content[k] += h->GetBinContent(b);
      target->variance[k] += h->GetBinError(b) * h->GetBinError(b);
    }
  }
  for (double v : s.pass.content) s.raw_pass += v;
  for (double v : s.fail.content) s.raw_fail += v;
  for (const auto* h : {&s.pass, &s.fail})
    for (std::size_t b = 0; b < h->content.size(); ++b) {
      s.sum_w += h->content[b];
      s.sum_w2 += h->variance[b];
    }
  return s;
}

// A smoothed signal template (contents on a uniform grid).
struct Template {
  double lo = 0, width = 0, entries = 0, pilot_width = 0, peak_width = 0, max_width = 0;
  std::vector<double> content;
  std::string source;
};

struct Smoothing {
  double pilot_scale = 0.9, alpha = 0.5, min_width = 0.0, max_width = 4.0, peak_center = 91.1876, peak_half_window = 10.0;
  double min_kernel_entries = 0.0;
};

double effective_entries(const TH1D* h) {
  double sw = 0, sw2 = 0;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    sw += h->GetBinContent(b);
    sw2 += h->GetBinError(b) * h->GetBinError(b);
  }
  return sw2 > 0 ? sw * sw / sw2 : 0.0;
}

// Adaptive (sample-point) kernel smoothing of a template histogram (Abramson 1982).  A pilot density with the
// fixed bandwidth h0 = pilot_scale min(sigma, IQR / 1.34) n^(-1/5) (Silverman's rule) of the Z peak: sigma, IQR
// and the effective entries n of the template within peak_half_window of peak_center (the global spread of a
// template with a large low-mass continuum or FSR hump would oversmooth its peak); then every bin's content
// spread with its own Gaussian width h_k = h0 (f0(x_k) / g)^(-alpha), g the geometric mean of the pilot density
// over the entries, clamped to [min_width, max_width].  The width is narrow where the template is dense (the
// peak of a high-statistics template is not broadened, which the fitted Gaussian could not undo) and wide in
// its sparse tails.  Every width is at least the half-width of the smallest symmetric interval around the bin
// holding min_kernel_entries effective entries (a nearest-neighbour floor: the relative statistical fluctuation
// of the smoothed template stays below about 1/sqrt(min_kernel_entries) wherever the template is sparse, while
// dense templates keep their narrow kernels).  (Separate pilot bandwidths for the peak and the continuum were
// tried and rejected: the bandwidth step at the region boundary distorted the templates; a stronger Abramson
// sensitivity alpha = 1 oversmoothed the tails of dense pass templates.)  The contents keep their sum; negative
// (NLO-weight) sums per bin are spread like the others and the result is clamped at zero.
Template make_template(const TH1D* h, const Smoothing& cfg, const std::string& source) {
  Template t;
  t.lo = h->GetXaxis()->GetXmin();
  t.width = h->GetXaxis()->GetBinWidth(1);
  t.source = source;
  t.entries = effective_entries(h);
  const int n = h->GetNbinsX();
  std::vector<double> raw(n), x(n);
  double total = 0, peak_total = 0, mean = 0, mean2 = 0;
  for (int b = 1; b <= n; ++b) {
    raw[b - 1] = h->GetBinContent(b);
    x[b - 1] = h->GetBinCenter(b);
    const double w = std::max(raw[b - 1], 0.0);
    total += w;
    if (std::fabs(x[b - 1] - cfg.peak_center) >= cfg.peak_half_window) continue;
    peak_total += w;
    mean += w * x[b - 1];
    mean2 += w * x[b - 1] * x[b - 1];
  }
  if (peak_total <= 0) throw std::runtime_error("template " + source + " has no entries near the Z peak");
  mean /= peak_total;
  const double sigma = std::sqrt(std::max(mean2 / peak_total - mean * mean, 0.0));
  auto quantile = [&](double q) {
    double cumulative = 0;
    for (int k = 0; k < n; ++k) {
      if (std::fabs(x[k] - cfg.peak_center) >= cfg.peak_half_window) continue;
      const double w = std::max(raw[k], 0.0);
      if (w > 0 && cumulative + w >= q * peak_total) return x[k] - 0.5 * t.width + t.width * (q * peak_total - cumulative) / w;
      cumulative += w;
    }
    return cfg.peak_center + cfg.peak_half_window;
  };
  const double iqr = quantile(0.75) - quantile(0.25);
  const double spread = iqr > 0 ? std::min(sigma, iqr / 1.34) : sigma;
  const double h0 = std::clamp(cfg.pilot_scale * spread * std::pow(std::max(t.entries * peak_total / total, 1.0), -0.2),
                               cfg.min_width, cfg.max_width);
  auto smooth = [&](const std::vector<double>& widths) {
    std::vector<double> out(n, 0.0);
    for (int k = 0; k < n; ++k) {
      if (raw[k] == 0) continue;
      if (widths[k] < 0.25 * t.width) {
        out[k] += raw[k];
        continue;
      }
      const int reach = static_cast<int>(std::ceil(5.0 * widths[k] / t.width));
      const int j0 = std::max(0, k - reach), j1 = std::min(n - 1, k + reach);
      double norm = 0;
      for (int j = j0; j <= j1; ++j) norm += std::exp(-0.5 * std::pow((j - k) * t.width / widths[k], 2));
      for (int j = j0; j <= j1; ++j) out[j] += raw[k] * std::exp(-0.5 * std::pow((j - k) * t.width / widths[k], 2)) / norm;
    }
    return out;
  };
  const std::vector<double> pilot = smooth(std::vector<double>(n, h0));
  double log_g = 0, weight = 0;
  for (int k = 0; k < n; ++k)
    if (raw[k] > 0 && pilot[k] > 0) {
      log_g += raw[k] * std::log(pilot[k]);
      weight += raw[k];
    }
  const double g = std::exp(log_g / weight);
  std::vector<double> widths(n, cfg.max_width);
  int peak = 0;
  // Effective entries per bin (the contents scaled to the effective statistics of the template).
  const double per_content = total > 0 ? t.entries / total : 0.0;
  for (int k = 0; k < n; ++k) {
    if (pilot[k] > 0) widths[k] = std::clamp(h0 * std::pow(pilot[k] / g, -cfg.alpha), cfg.min_width, cfg.max_width);
    if (cfg.min_kernel_entries > 0) {
      double count = std::max(raw[k], 0.0) * per_content;
      int half = 0;
      while (count < cfg.min_kernel_entries && half < n) {
        ++half;
        if (k - half >= 0) count += std::max(raw[k - half], 0.0) * per_content;
        if (k + half < n) count += std::max(raw[k + half], 0.0) * per_content;
      }
      widths[k] = std::min(std::max(widths[k], (half + 0.5) * t.width), cfg.max_width);
    }
    if (pilot[k] > pilot[peak]) peak = k;
  }
  t.content = smooth(widths);
  for (double& c : t.content) c = std::max(c, 0.0);
  t.pilot_width = h0;
  t.peak_width = widths[peak];
  t.max_width = *std::max_element(widths.begin(), widths.end());
  return t;
}

struct Seed {
  double shift = 0, sigma = 1.2, alpha_l = 1.0, n_l = 4.0, alpha_r = 1.5, n_r = 8.0, mean_pt = 40, mean_eta = 1.0;
  double pass_bkg[3] = {60.0, 0.04, 0.05}, fail_bkg[3] = {60.0, 0.04, 0.05};
  bool valid = false;
};

struct Settings {
  double z_mass = 91.1876, z_width = 2.4952, max_edm = 0.01, width_a = 0.76, width_b = 1.0, shared_below = 300,
         small_category = 100;
  // CMSShape ranges [low, high] of alpha, beta, gamma and the fail low-mass Gaussian of the DSCB alternative.
  double cms_range[3][2] = {{50.0, 80.0}, {0.01, 0.06}, {0.005, 1.0}};
  double hump_mean[2] = {65.0, 88.0}, hump_width[2] = {3.0, 12.0};
  // Optional pass-like admixture of the fail template (fraction of the fail signal with the pass template's shape).
  bool fail_pass_like = false;
  // The DSCB alternative: a low-mass Gaussian in the pass as well.
  bool pass_hump = false;
  double pass_like_range[2] = {-0.3, 0.9};
};

// The manual fix of one fit (configuration "overrides").
struct Override {
  bool active = false;
  double window_low = -1, window_high = -1, min_width = 0.5;
  std::map<std::string, double> start, fix;
  std::string note;
};

Override find_override(const json& list, const std::string& flavour, const std::string& step, const std::string& role,
                       const std::string& bin, const std::string& model) {
  Override o;
  for (const auto& item : list) {
    if (item.at("flavour") != flavour || item.at("step") != step || item.at("role") != role || item.at("bin") != bin) continue;
    if (item.contains("models")) {
      bool listed = false;
      for (const auto& m : item.at("models")) listed = listed || m == model;
      if (!listed) continue;
    }
    o.active = true;
    if (item.contains("window")) {
      o.window_low = item.at("window").at(0).get<double>();
      o.window_high = item.at("window").at(1).get<double>();
    }
    o.min_width = item.value("min_width", 0.5);
    if (item.contains("start"))
      for (const auto& [k, v] : item.at("start").items()) o.start[k] = v.get<double>();
    if (item.contains("fix"))
      for (const auto& [k, v] : item.at("fix").items()) o.fix[k] = v.get<double>();
    o.note = item.value("note", "");
  }
  return o;
}

// Parameter layout: 0 shift_pass, 1 sigma_pass, 2 shift_fail, 3 sigma_fail, 4 alpha_l, 5 n_l, 6 alpha_r, 7 n_r
// (DSCB only), 8 efficiency, 9 n_sig, 10 b_pass, 11 b_fail, 12-14 pass background shape, 15-17 fail background
// shape, 18-20 the fail low-mass Gaussian of the DSCB alternative (mean, width, fraction).  Parameters that may end on a physical bound (no extra smearing, no background, a flat background, a
// turn-on below the window, a vanishing Bernstein coefficient) are not reported as "at limit".
json fit_bin(const Spectrum& s, const std::string& model, const Template* tpass, const Template* tfail, const Seed& seed,
             double sigma_start, bool inclusive, bool is_mc, const Settings& set, const Override& ov, Seed* result_seed,
             bool curves) {
  const bool dscb = model == "alt_signal";
  const std::string bkg = model == "alt_background" ? "bernstein3" : "cmsshape";
  const bool cms = bkg == "cmsshape";
  const double step = dscb ? 0.05 : 0.125;
  h4l::ZPeakModel shape_pass(s.lo, s.width, s.bins, set.z_mass, set.z_width, bkg, step, 30.0);
  h4l::ZPeakModel shape_fail(s.lo, s.width, s.bins, set.z_mass, set.z_width, bkg, step, 30.0);
  if (dscb) {
    // The pass DSCB of the alternative has its own low-mass Gaussian (the radiative continuum below the peak)
    // when configured; without it the pass background absorbed that continuum (v6 inspection).
    shape_pass.kernel_type = set.pass_hump ? h4l::ZPeakModel::Kernel::StandaloneDCBHump : h4l::ZPeakModel::Kernel::StandaloneDCB;
    shape_fail.kernel_type = h4l::ZPeakModel::Kernel::StandaloneDCBHump;
  } else {
    shape_pass.kernel_type = shape_fail.kernel_type = h4l::ZPeakModel::Kernel::Template;
    shape_pass.set_template(tpass->lo, tpass->width, tpass->content);
    shape_fail.set_template(tfail->lo, tfail->width, tfail->content);
    if (set.fail_pass_like) shape_fail.set_second_template(tpass->content);
  }
  const double n = s.raw_pass + s.raw_fail;
  const double eff_probes = s.sum_w2 > 0 ? s.sum_w * s.sum_w / s.sum_w2 : 0.0;
  const bool shared_resolution = s.raw_fail < set.shared_below || s.raw_pass < set.shared_below;
  const bool small_pass = !inclusive && s.raw_pass < set.small_category, small_fail = !inclusive && s.raw_fail < set.small_category;
  bool fix_n = false, fix_alpha = false, fix_core = false;
  if (dscb && !inclusive && seed.valid) {
    fix_n = eff_probes < 1000;
    fix_alpha = eff_probes < 300;
    fix_core = eff_probes < 100;
  }
  double sigma_low = 0.02, sigma_high = 4.0, sigma0 = is_mc ? 0.1 : 0.5;
  if (dscb) {
    sigma_low = std::max(0.3, sigma_start / 3);
    sigma_high = std::min(8.0, sigma_start * 3);
    sigma0 = std::clamp(sigma_start, sigma_low, sigma_high);
  } else if (!inclusive && seed.valid) {
    sigma0 = std::clamp(seed.sigma, sigma_low, sigma_high);
  }
  const double eps0 = std::clamp(s.raw_pass / std::max(n, 1e-9), 0.02, 0.98);
  auto bkg_parameters = [&](const char* prefix, const double* start, bool fixed) {
    std::vector<h4l::FitParameter> q;
    const std::string p(prefix);
    if (cms) {
      const auto& r = set.cms_range;
      q.push_back({p + "_alpha", std::clamp(start[0], r[0][0], r[0][1]), 1.0, r[0][0], r[0][1], fixed});
      q.push_back({p + "_beta", std::clamp(start[1], r[1][0], r[1][1]), 0.002, r[1][0], r[1][1], fixed});
      q.push_back({p + "_gamma", std::clamp(start[2], r[2][0], r[2][1]), 0.01, r[2][0], r[2][1], fixed});
    } else {
      for (int k = 1; k <= 3; ++k) q.push_back({p + "_c" + std::to_string(k), 1.0, 0.1, 0.0, 20.0, false});
    }
    return q;
  };
  auto parameters = [&](const double* pass_bkg, const double* fail_bkg) {
    std::vector<h4l::FitParameter> p = {
        {"shift_pass", dscb ? seed.shift : 0.0, 0.05, -3.0, 3.0, fix_core},
        {"sigma_pass", sigma0, 0.05, sigma_low, sigma_high, fix_core},
        {"shift_fail", dscb ? seed.shift : 0.0, 0.05, -3.0, 3.0, fix_core || shared_resolution},
        {"sigma_fail", dscb ? std::clamp(1.1 * sigma0, sigma_low, sigma_high) : sigma0, 0.05, sigma_low, sigma_high,
         fix_core || shared_resolution},
        {"alpha_l", seed.valid ? seed.alpha_l : 1.0, 0.05, 0.2, 10.0, !dscb || fix_alpha},
        {"n_l", seed.valid ? seed.n_l : 4.0, 0.2, 1.01, 80.0, !dscb || fix_n},
        {"alpha_r", seed.valid ? seed.alpha_r : 1.5, 0.05, 0.2, 10.0, !dscb || fix_alpha},
        {"n_r", seed.valid ? seed.n_r : 8.0, 0.2, 1.01, 80.0, !dscb || fix_n},
        {"efficiency", eps0, 0.01, 0.0, 1.0, false},
        {"n_sig", 0.95 * n, 0.01 * n + 1, 0.0, 3.0 * n + 10, false},
        {"b_pass", 0.05 * s.raw_pass + 1, 0.02 * s.raw_pass + 1, 0.0, 3.0 * s.raw_pass + 10, false},
        {"b_fail", 0.05 * s.raw_fail + 1, 0.02 * s.raw_fail + 1, 0.0, 3.0 * s.raw_fail + 10, false}};
    for (auto& q : bkg_parameters("pass", pass_bkg, cms && small_pass && seed.valid)) p.push_back(q);
    for (auto& q : bkg_parameters("fail", fail_bkg, cms && small_fail && seed.valid)) p.push_back(q);
    p.push_back({"hump_mean_fail", std::clamp(75.0, set.hump_mean[0], set.hump_mean[1]), 1.0, set.hump_mean[0], set.hump_mean[1],
                 !dscb});
    p.push_back({"hump_width_fail", std::clamp(5.0, set.hump_width[0], set.hump_width[1]), 0.2, set.hump_width[0],
                 set.hump_width[1], !dscb});
    if (dscb) p.push_back({"hump_fraction_fail", 0.1, 0.02, 0.0, 0.95, false});
    else p.push_back({"pass_like_fraction_fail", 0.0, 0.02, set.pass_like_range[0], set.pass_like_range[1], !set.fail_pass_like});
    const bool pass_hump = dscb && set.pass_hump;
    p.push_back({"hump_mean_pass", std::clamp(75.0, set.hump_mean[0], set.hump_mean[1]), 1.0, set.hump_mean[0], set.hump_mean[1],
                 !pass_hump});
    p.push_back({"hump_width_pass", std::clamp(8.0, set.hump_width[0], set.hump_width[1]), 0.2, set.hump_width[0], set.hump_width[1],
                 !pass_hump});
    p.push_back({"hump_fraction_pass", pass_hump ? 0.05 : 0.0, 0.01, 0.0, 0.5, !pass_hump});
    for (auto& q : p) {
      if (ov.start.count(q.name)) q.value = std::clamp(ov.start.at(q.name), q.low, q.high);
      if (ov.fix.count(q.name)) {
        q.value = ov.fix.at(q.name);
        q.fixed = true;
      }
    }
    return p;
  };
  const std::vector<h4l::FitHistogram> data = {s.pass, s.fail};
  auto evaluate = [&](const double* x, std::vector<std::vector<double>>& expected, std::vector<double>* bkg_pass = nullptr,
                      std::vector<double>* bkg_fail = nullptr) {
    const double fs = shared_resolution ? x[0] : x[2], fw = shared_resolution ? x[1] : x[3];
    const double pass_par[14] = {x[0], x[1], x[4], x[5], x[6], x[7], x[8] * x[9], x[10], x[12], x[13], x[14],
                                 x[21], x[22], x[23]};
    const double fail_par[14] = {fs, fw, x[4], x[5], x[6], x[7], (1.0 - x[8]) * x[9], x[11], x[15], x[16], x[17],
                                 x[18], x[19], x[20]};
    shape_pass.expected(pass_par, expected[0], bkg_pass);
    shape_fail.expected(fail_par, expected[1], bkg_fail);
  };
  auto model_fn = [&](const double* x, std::vector<std::vector<double>>& expected) { evaluate(x, expected); };
  auto converged = [&](const h4l::FitOutcome& o) {
    return (o.status == 0 || o.status == 1 || (o.status == 3 && o.edm < set.max_edm)) && std::isfinite(o.nll);
  };
  // Starting points: the background shapes of the inclusive fit (when available) and the default shapes.
  const double default_bkg[3] = {60.0, 0.04, 0.05};
  std::vector<std::vector<h4l::FitParameter>> starts;
  if (seed.valid && cms) starts.push_back(parameters(seed.pass_bkg, seed.fail_bkg));
  starts.push_back(parameters(default_bkg, default_bkg));
  std::vector<h4l::FitParameter> p = starts.front();
  h4l::FitOutcome outcome;
  bool have = false;
  int start_used = -1;
  for (std::size_t k = 0; k < starts.size(); ++k) {
    const h4l::FitOutcome o = h4l::binned_fit(data, starts[k], model_fn, 8, h4l::WeightedErrors::Effective);
    const bool better = !have || (converged(o) && !converged(outcome)) ||
                        (converged(o) == converged(outcome) && o.nll < outcome.nll - 1e-6);
    if (better) {
      outcome = o;
      p = starts[k];
      have = true;
      start_used = static_cast<int>(k);
    }
  }
  // DSCB alternative: a width resting on its bound is refitted with the bounds opened to [0.3, 8] GeV.
  bool relaxed = false;
  if (dscb) {
    for (int w : {1, 3}) {
      if (p[w].fixed) continue;
      const double range = p[w].high - p[w].low;
      const bool at_bound = outcome.values[w] - p[w].low < 1e-3 * range || p[w].high - outcome.values[w] < 1e-3 * range;
      if (!at_bound || (p[w].low <= 0.3 + 1e-9 && p[w].high >= 8.0 - 1e-9)) continue;
      std::vector<h4l::FitParameter> q = p;
      for (int v : {1, 3}) {
        q[v].low = 0.3;
        q[v].high = 8.0;
      }
      for (std::size_t i = 0; i < q.size(); ++i)
        if (!q[i].fixed) q[i].value = std::clamp(outcome.values[i], q[i].low, q[i].high);
      const h4l::FitOutcome retry = h4l::binned_fit(data, q, model_fn, 8, h4l::WeightedErrors::Effective);
      if (converged(retry) && retry.nll < outcome.nll) {
        outcome = retry;
        p = q;
        relaxed = true;
      }
      break;
    }
  }
  // The efficiency uncertainty: MINOS, else HESSE, else a MINOS interval flagged invalid by Minuit (one side at
  // the physical bound eps = 1) that still brackets the minimum; floored at the binomial error of the fitted signal.
  const double eps = outcome.values[8];
  double error = std::nan("");
  std::string error_source = "none";
  if (outcome.minos_ok) {
    error = 0.5 * (outcome.minos_high - outcome.minos_low);
    error_source = "minos";
  } else if (std::isfinite(outcome.errors[8]) && outcome.errors[8] > 0) {
    error = outcome.errors[8];
    error_source = "hesse";
  } else if (std::isfinite(outcome.minos_low) && std::isfinite(outcome.minos_high) && outcome.minos_low < 0 &&
             outcome.minos_high > 0) {
    error = 0.5 * (outcome.minos_high - outcome.minos_low);
    error_source = "minos_partial";
  }
  const bool error_valid = std::isfinite(error) && error > 0;
  const double n_sig_eff = outcome.values[9] * (n > 0 ? eff_probes / n : 0.0);
  const double floor = n_sig_eff > 0 ? std::sqrt(std::max(eps * (1 - eps), 0.0) / n_sig_eff) : 0.0;
  bool floor_applied = false;
  if (error_valid && error < 0.9 * floor) {
    error = floor;
    floor_applied = true;
  }
  json out;
  out["model"] = model;
  out["background_model"] = bkg;
  out["status"] = outcome.status;
  out["cov_status"] = outcome.cov_status;
  out["edm"] = outcome.edm;
  out["nll"] = outcome.nll;
  out["efficiency"] = eps;
  out["efficiency_error"] = error;
  out["efficiency_error_hesse"] = outcome.errors[8];
  out["efficiency_minos"] = {{"ok", outcome.minos_ok}, {"low", outcome.minos_low}, {"high", outcome.minos_high}};
  out["binomial_floor"] = floor;
  out["error_floor_applied"] = floor_applied;
  out["error_source"] = error_source;
  out["effective_probes"] = eff_probes;
  out["shared_resolution"] = shared_resolution;
  out["small_category_background_fixed"] = {{"pass", cms && small_pass && seed.valid}, {"fail", cms && small_fail && seed.valid}};
  out["sigma_start"] = sigma_start;
  out["start_used"] = start_used;
  out["width_bound_relaxed"] = relaxed;
  out["bin_width"] = s.width;
  out["window"] = {s.lo, s.lo + s.width * s.bins};
  if (!dscb) {
    auto describe = [](const Template* t) {
      return json{{"source", t->source}, {"entries", t->entries}, {"pilot_width", t->pilot_width}, {"peak_width", t->peak_width},
                  {"max_width", t->max_width}};
    };
    out["templates"] = {{"pass", describe(tpass)}, {"fail", describe(tfail)}};
  }
  if (ov.active) out["override"] = ov.note;
  // Parameters at a non-physical limit (a yield within half an event of a bound counts as at the bound).
  std::set<std::string> physical_low = {"b_pass", "b_fail", "pass_gamma", "fail_gamma", "pass_alpha", "fail_alpha", "pass_beta",
                                        "fail_beta", "pass_c1", "pass_c2", "pass_c3", "fail_c1", "fail_c2", "fail_c3",
                                        "hump_fraction_fail", "hump_fraction_pass"};
  if (!dscb) {
    physical_low.insert("sigma_pass");
    physical_low.insert("sigma_fail");
  }
  // Background shape parameters matter only where the category has a background (> 2 % of its entries).
  const bool pass_background = outcome.values[10] > 0.02 * s.raw_pass, fail_background = outcome.values[11] > 0.02 * s.raw_fail;
  json at_limit = json::array();
  for (std::size_t i = 0; i < p.size(); ++i) {
    if (p[i].fixed) continue;
    if ((i >= 12 && i <= 14 && !pass_background) || (i >= 15 && i <= 17 && !fail_background)) continue;
    const double tolerance = (i >= 9 && i <= 11) ? 0.5 : 1e-4 * (p[i].high - p[i].low);
    const bool low = outcome.values[i] - p[i].low < tolerance, high = p[i].high - outcome.values[i] < tolerance;
    if ((low && !physical_low.count(p[i].name)) || high) at_limit.push_back(p[i].name);
  }
  out["at_limit"] = at_limit;
  json parameters_out;
  for (std::size_t i = 0; i < p.size(); ++i) parameters_out[p[i].name] = {outcome.values[i], outcome.errors[i]};
  out["parameters"] = parameters_out;
  // A background with a narrow maximum near m_Z patches the signal.
  std::vector<std::vector<double>> expected(2), background(2);
  evaluate(outcome.values.data(), expected, &background[0], &background[1]);
  bool narrow_background_peak = false;
  for (int c = 0; c < 2; ++c) {
    const auto peak = std::max_element(background[c].begin(), background[c].end());
    const int b = static_cast<int>(peak - background[c].begin());
    const double x = s.lo + (b + 0.5) * s.width;
    if (!(*peak > 0) || b == 0 || b + 1 >= s.bins || std::fabs(x - set.z_mass) > 12.0) continue;
    const int reach = static_cast<int>(std::lround(8.0 / s.width));
    const double left = background[c][std::max(0, b - reach)], right = background[c][std::min(s.bins - 1, b + reach)];
    if (left < 0.7 * *peak && right < 0.7 * *peak) narrow_background_peak = true;
  }
  out["background_narrow_peak_near_z"] = narrow_background_peak;
  // Diagnostic chi2 with the variance floored at the expectation times the mean weight; never an acceptance criterion.
  const double mean_weight = s.sum_w > 0 ? s.sum_w2 / s.sum_w : 1.0;
  double chi2 = 0;
  int used = 0;
  for (int c = 0; c < 2; ++c)
    for (int b = 0; b < s.bins; ++b) {
      const double variance = std::max(data[c].variance[b], expected[c][b] * mean_weight);
      if (variance > 0) {
        chi2 += std::pow(data[c].content[b] - expected[c][b], 2) / variance;
        ++used;
      }
    }
  out["chi2"] = chi2;
  out["ndf"] = used - static_cast<int>(std::count_if(p.begin(), p.end(), [](const auto& q) { return !q.fixed; }));
  out["error_valid"] = error_valid;
  out["converged"] = converged(outcome);
  out["ok"] = converged(outcome) && outcome.minos_ok && at_limit.empty() && !narrow_background_peak && error_valid && !floor_applied;
  out["usable"] = converged(outcome) && !narrow_background_peak && error_valid;
  if (curves) {
    json x = json::array();
    std::vector<double> pass_error(s.bins), fail_error(s.bins);
    for (int b = 0; b < s.bins; ++b) {
      x.push_back(s.lo + (b + 0.5) * s.width);
      // The bin errors are the square roots of the (sum of squared weights) variances.
      pass_error[b] = std::sqrt(std::max(s.pass.variance[b], 0.0));
      fail_error[b] = std::sqrt(std::max(s.fail.variance[b], 0.0));
    }
    out["curve"] = {{"x", x}, {"pass", s.pass.content}, {"fail", s.fail.content},
                    {"pass_error", pass_error}, {"fail_error", fail_error},
                    {"model_pass", expected[0]}, {"model_fail", expected[1]},
                    {"background_pass", background[0]}, {"background_fail", background[1]}};
  }
  if (result_seed) {
    result_seed->shift = outcome.values[0];
    result_seed->sigma = outcome.values[1];
    result_seed->alpha_l = outcome.values[4];
    result_seed->n_l = outcome.values[5];
    result_seed->alpha_r = outcome.values[6];
    result_seed->n_r = outcome.values[7];
    if (cms)
      for (int k = 0; k < 3; ++k) {
        result_seed->pass_bkg[k] = outcome.values[12 + k];
        result_seed->fail_bkg[k] = outcome.values[15 + k];
      }
    result_seed->valid = converged(outcome);
  }
  return out;
}

double counting_efficiency(const TH1D* pass, const TH1D* fail, double lo, double hi, double& error) {
  auto integral = [&](const TH1D* h, double& variance) {
    double sum = 0;
    variance = 0;
    if (!h) return 0.0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      const double x = h->GetBinCenter(b);
      if (x < lo || x >= hi) continue;
      sum += h->GetBinContent(b);
      variance += h->GetBinError(b) * h->GetBinError(b);
    }
    return sum;
  };
  double vp = 0, vf = 0;
  const double np = integral(pass, vp), nf = integral(fail, vf), n = np + nf;
  if (!(n > 0)) {
    error = 0;
    return std::nan("");
  }
  const double eps = np / n;
  error = std::sqrt((1 - eps) * (1 - eps) * vp + eps * eps * vf) / n;
  return eps;
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
    if (key == "--job" || key == "--task") job_path = argv[index + 1];
    else if (key == "--out" || key == "--out-json") out_path = argv[index + 1];
    else {
      std::cerr << "unknown argument " << key << "\n";
      return 64;
    }
  }
  if (job_path.empty() || out_path.empty()) {
    std::cerr << "usage: fit_tnp --job JOB.json --out OUT.json\n";
    return 64;
  }
  gROOT->SetBatch(true);
  TH1::AddDirectory(false);
  try {
    const json job = h4l::read_json(job_path);
    const json& fit_cfg = job.at("tnp_config").at("fit");
    const std::string tag = job.at("flavour").get<std::string>(), step = job.at("step").get<std::string>();
    const int n_eta = job.at("n_eta").get<int>();
    const std::vector<double> pt_edges = job.at("pt_edges").get<std::vector<double>>();
    const int n_pt = static_cast<int>(pt_edges.size());
    const double low = fit_cfg.at("window").at(0).get<double>(), high = fit_cfg.at("window").at(1).get<double>();
    const int min_entries = fit_cfg.at("min_category_entries").get<int>();
    const json& template_cfg = fit_cfg.at("template");
    const double template_min = template_cfg.at("min_effective_entries").get<double>();
    // The |eta|-column merges, tried in order: the neighbouring pT rows whose lower pT edge lies within each factor
    // of the bin's (a single number is one factor).
    std::vector<double> merge_pt_factors;
    if (template_cfg.at("merge_pt_factor").is_array())
      for (const auto& f : template_cfg.at("merge_pt_factor")) merge_pt_factors.push_back(f.get<double>());
    else
      merge_pt_factors.push_back(template_cfg.at("merge_pt_factor").get<double>());
    // A template must carry at least this many times the effective statistics of the data spectrum it describes
    // (a template as sparse as the data puts its statistical fluctuations into the model).
    const double template_to_data = template_cfg.value("min_template_to_data_ratio", 0.0);
    Smoothing smoothing;
    smoothing.pilot_scale = template_cfg.at("smoothing").at("pilot_scale").get<double>();
    smoothing.alpha = template_cfg.at("smoothing").at("alpha").get<double>();
    smoothing.min_width = template_cfg.at("smoothing").at("min_width").get<double>();
    smoothing.max_width = template_cfg.at("smoothing").at("max_width").get<double>();
    smoothing.peak_half_window = template_cfg.at("smoothing").at("peak_half_window").get<double>();
    smoothing.peak_center = fit_cfg.at("z_mass").get<double>();
    smoothing.min_kernel_entries = template_cfg.at("smoothing").value("min_kernel_entries", 0.0);
    Settings set;
    const char* cms_names[3] = {"alpha", "beta", "gamma"};
    for (int k = 0; k < 3; ++k)
      for (int e = 0; e < 2; ++e) set.cms_range[k][e] = fit_cfg.at("cmsshape").at(cms_names[k]).at(e).get<double>();
    for (int e = 0; e < 2; ++e) {
      set.hump_mean[e] = fit_cfg.at("alt_signal_hump").at("mean").at(e).get<double>();
      set.hump_width[e] = fit_cfg.at("alt_signal_hump").at("width").at(e).get<double>();
    }
    set.pass_hump = fit_cfg.at("alt_signal_hump").value("pass", false);
    set.z_mass = fit_cfg.at("z_mass").get<double>();
    set.z_width = fit_cfg.at("z_width").get<double>();
    set.max_edm = fit_cfg.at("max_edm").get<double>();
    set.width_a = fit_cfg.at("start_width").at(tag).at("a").get<double>();
    set.width_b = fit_cfg.at("start_width").at(tag).at("b").get<double>();
    set.shared_below = fit_cfg.at("shared_resolution_below").get<double>();
    set.small_category = fit_cfg.at("small_category_below").get<double>();
    if (template_cfg.contains("fail_pass_like")) {
      set.fail_pass_like = template_cfg.at("fail_pass_like").at("enabled").get<bool>();
      set.pass_like_range[0] = template_cfg.at("fail_pass_like").at("range").at(0).get<double>();
      set.pass_like_range[1] = template_cfg.at("fail_pass_like").at("range").at(1).get<double>();
    }
    const json overrides = fit_cfg.value("overrides", json::array());
    std::set<std::string> seen;
    auto open_all = [&](const char* role) {
      std::vector<std::unique_ptr<TFile>> files;
      for (const auto& path : job.at(role)) {
        if (!seen.insert(path.get<std::string>()).second) throw std::runtime_error("duplicate input " + path.get<std::string>());
        files.push_back(open_retrying(path.get<std::string>()));
      }
      return files;
    };
    const auto data_files = open_all("data");
    const auto mc_files = open_all("mc");
    // Signal template of a bin and category: the bin's own when populated enough (min_effective_entries), else
    // the bin's |eta| column with the neighbouring pT rows added nearest first (the resolution, which the fitted
    // Gaussian can only broaden, depends mostly on |eta|) as long as their lower pT edges stay within each
    // merge_pt_factor of the bin's (tried in order), else the bin's pT row, else the inclusive one.
    auto template_of = [&](const std::string& kind, const std::string& bin) {
      const std::string prefix = kind + "_" + tag + "_" + step + "_";
      if (bin == "all") {
        auto h = summed(mc_files, prefix + "all");
        if (!h) throw std::runtime_error("no signal template " + prefix + "all");
        return make_template(h.get(), smoothing, "all");
      }
      // The data spectrum of this bin and category (same for the data and the MC fits of the bin).
      auto data_spectrum = summed(data_files, std::string(kind == "tpass" ? "pass_" : "fail_") + tag + "_" + step + "_" + bin);
      const double needed = std::max(template_min, template_to_data * (data_spectrum ? effective_entries(data_spectrum.get()) : 0.0));
      const int index = std::stoi(bin), row = index / n_eta, col = index % n_eta;
      std::unique_ptr<TH1D> merged = summed(mc_files, prefix + bin);
      if (merged && effective_entries(merged.get()) >= needed) return make_template(merged.get(), smoothing, bin);
      for (double factor : merge_pt_factors) {
        std::unique_ptr<TH1D> column = summed(mc_files, prefix + bin);
        int lo_row = row, hi_row = row;
        for (int distance = 1; distance < n_pt; ++distance) {
          for (int r : {row - distance, row + distance}) {
            if (r < 0 || r >= n_pt) continue;
            if (pt_edges[r] > factor * pt_edges[row] || pt_edges[r] * factor < pt_edges[row]) continue;
            auto h = summed(mc_files, prefix + std::to_string(r * n_eta + col));
            if (!h) continue;
            if (!column) column = std::move(h);
            else if (!column->Add(h.get())) throw std::runtime_error("cannot add template " + prefix + std::to_string(r * n_eta + col));
            lo_row = std::min(lo_row, r);
            hi_row = std::max(hi_row, r);
            if (effective_entries(column.get()) >= needed)
              return make_template(column.get(), smoothing,
                                   "eta" + std::to_string(col) + "_pt" + std::to_string(lo_row) + "-" + std::to_string(hi_row));
          }
        }
      }
      for (const std::string& label : {"pt" + std::to_string(row), std::string("all")}) {
        auto h = summed(mc_files, prefix + label);
        if (!h) continue;
        if (effective_entries(h.get()) >= needed || label == "all") return make_template(h.get(), smoothing, label);
      }
      throw std::runtime_error("no signal template " + prefix + bin);
    };
    json results = json::array();
    for (const auto& [role, files] : {std::pair{"data", &data_files}, std::pair{"mc", &mc_files}}) {
      Seed seeds[3];
      const bool is_mc = std::string(role) == "mc";
      for (const auto& bin_node : job.at("bins")) {
        const std::string bin = bin_node.get<std::string>();
        const std::string suffix = tag + "_" + step + "_" + bin;
        auto pass = summed(*files, "pass_" + suffix), fail = summed(*files, "fail_" + suffix);
        auto pstat = summed(*files, "pstat_" + tag + "_" + bin);
        json entry = {{"role", role}, {"flavour", tag}, {"step", step}, {"bin", bin}};
        double count_error = 0;
        entry["counting"] = {{"efficiency", counting_efficiency(pass.get(), fail.get(), 80, 100, count_error)},
                             {"error", count_error}};
        double mean_pt = 40, mean_eta = 1.0;
        if (pstat && pstat->GetBinContent(1) > 0) {
          mean_pt = pstat->GetBinContent(2) / pstat->GetBinContent(1);
          mean_eta = pstat->GetBinContent(3) / pstat->GetBinContent(1);
        }
        entry["mean_pt"] = mean_pt;
        entry["mean_abs_eta"] = mean_eta;
        if (is_mc) {
          auto truth = summed(*files, "truth_" + suffix);
          if (truth) {
            const double f = truth->GetBinContent(1), p = truth->GetBinContent(2);
            entry["truth_diagnostic"] = {
                {"efficiency", p + f > 0 ? p / (p + f) : std::nan("")},
                {"error", p + f > 0 ? std::sqrt(std::pow(truth->GetBinError(2) * f, 2) + std::pow(truth->GetBinError(1) * p, 2)) /
                                          std::pow(p + f, 2)
                                    : 0.0}};
          }
        }
        if (!pass && !fail) {
          entry["status"] = "empty";
          results.push_back(entry);
          continue;
        }
        const Template tpass = template_of("tpass", bin), tfail = template_of("tfail", bin);
        const bool inclusive = bin == "all";
        json fits;
        const char* models[] = {"nominal", "alt_signal", "alt_background"};
        bool counting_only = false;
        for (int k = 0; k < 3; ++k) {
          const Override ov = find_override(overrides, tag, step, role, bin, models[k]);
          const double w_low = ov.window_low > 0 ? ov.window_low : low, w_high = ov.window_high > 0 ? ov.window_high : high;
          const Spectrum s = make_spectrum(pass.get(), fail.get(), w_low, w_high, ov.min_width);
          if (k == 0) {
            entry["raw_pass"] = s.raw_pass;
            entry["raw_fail"] = s.raw_fail;
            if (s.raw_pass < min_entries || s.raw_fail < min_entries) {
              counting_only = true;
              break;
            }
          }
          const Seed seed = inclusive ? Seed{} : seeds[k];
          double sigma_start = seed.sigma;
          if (k == 1 && !inclusive && seed.valid) {
            const double w_bin = (1 + set.width_a * mean_eta) * (1 + set.width_b * mean_pt / 100.0);
            const double w_incl = (1 + set.width_a * seed.mean_eta) * (1 + set.width_b * seed.mean_pt / 100.0);
            sigma_start = std::clamp(seed.sigma * w_bin / w_incl, 0.3, 8.0);
          }
          Seed* store = inclusive ? &seeds[k] : nullptr;
          fits[models[k]] = fit_bin(s, models[k], &tpass, &tfail, seed, sigma_start, inclusive, is_mc, set, ov, store,
                                    job.value("curves", true));
          if (store) {
            store->mean_pt = mean_pt;
            store->mean_eta = mean_eta;
          }
        }
        if (counting_only) {
          entry["status"] = "counting_only";
          results.push_back(entry);
          continue;
        }
        entry["fits"] = fits;
        entry["status"] = fits["nominal"]["ok"].get<bool>() ? "ok" : (fits["nominal"]["usable"].get<bool>() ? "constrained" : "rejected");
        std::cout << "[tnp] " << role << " " << suffix << ": eps " << fits["nominal"]["efficiency"].get<double>() << " +- "
                  << fits["nominal"]["efficiency_error"] << " (" << entry["status"].get<std::string>() << ")" << std::endl;
        results.push_back(entry);
      }
    }
    json report = {{"schema", "h4l_v3_tnp_fits/3"}, {"flavour", tag}, {"step", step}, {"data", job.at("data")},
                   {"mc", job.at("mc")}, {"program", job.value("program", json())},
                   {"frozen_program", job.value("frozen_program", json())}, {"task_id", job.value("task_id", json())},
                   {"bins", job.at("bins")}, {"results", results},
                   {"finished_utc", h4l::utc_now()}};
    h4l::publish_json(out_path, report);
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
