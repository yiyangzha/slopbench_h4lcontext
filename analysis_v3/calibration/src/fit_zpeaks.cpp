// fit_zpeaks: BW (x) double-sided Crystal Ball + background fits of the
// Z -> ll mass histograms of the lepton calibration (AN-16-442 3.1.4, 5.3).
//
// Usage: fit_zpeaks --job JOB.json --out OUT.json
//
// The job names the histogram files of the data role and of the MC role
// (zpeak_histograms outputs; each role is the sum of its files), the flavour
// tag ("mm" or "ee"), the list of category histograms to fit and the fit
// variant.  For every category:
//   1. The window is seeded from the data role: the mode of the histogram
//      smoothed over the configured number of bins, searched in mode_search;
//      the window is [mode - below, mode + above], snapped to bin edges.  The
//      MC role is fitted in the same window.
//   2. The MC role is fitted with every parameter free: BW(m_Z, Gamma_Z
//      fixed) (x) DCB(delta, sigma, alphaL, nL, alphaR, nR) plus the
//      background model (h4l/zpeak_model.h), extended binned Poisson
//      likelihood with Minuit2 and HESSE; weighted histograms get the
//      sandwich covariance (h4l/binned_fit.h).
//   3. The data role is fitted with the DCB tails fixed to the MC values
//      (variant data_tails = "fixed_from_mc") or free ("free").
// Every fit reports its status, covariance status, parameters with errors,
// yields, the Pearson chi2 in the window and, when curves are requested, the
// histogram and the fitted model per bin.  Categories below the statistics
// thresholds or without a histogram in both roles are listed as skipped.

#include "h4l/binned_fit.h"
#include "h4l/io.h"
#include "h4l/zpeak_model.h"

#include <TError.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
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

struct Settings {
  double mode_low = 70, mode_high = 110, below = 25, above = 20, clip_low = 53, clip_high = 125, z_mass = 91.1876,
         z_width = 2.4952;
  int smooth_bins = 9;
  double grid_step = 0.05, kernel_half = 30;
  double min_data = 300, min_mc_effective = 300;
  std::string background = "exponential", data_tails = "fixed_from_mc";
  bool curves = true;
};

// Sum of one histogram over several files; nullptr if no file has it.
std::unique_ptr<TH1D> summed(const std::vector<std::unique_ptr<TFile>>& files, const std::string& name) {
  std::unique_ptr<TH1D> total;
  for (const auto& file : files) {
    auto* histogram = dynamic_cast<TH1D*>(file->Get(name.c_str()));
    if (!histogram) continue;
    if (!total) {
      total.reset(static_cast<TH1D*>(histogram->Clone((name + "_sum").c_str())));
      total->SetDirectory(nullptr);
    } else if (!total->Add(histogram)) {
      throw std::runtime_error("cannot add " + name);
    }
  }
  return total;
}

// Mode of the histogram smoothed with a running mean over `width` bins,
// searched in [low, high]; also the half-maximum width of the smoothed peak.
std::pair<double, double> smoothed_mode(const TH1D& h, double low, double high, int width) {
  const int n = h.GetNbinsX(), half = width / 2;
  std::vector<double> smooth(n + 2, 0.0);
  for (int b = 1; b <= n; ++b) {
    double sum = 0;
    int count = 0;
    for (int k = std::max(1, b - half); k <= std::min(n, b + half); ++k) {
      sum += h.GetBinContent(k);
      ++count;
    }
    smooth[b] = sum / count;
  }
  int best = -1;
  for (int b = 1; b <= n; ++b) {
    const double x = h.GetBinCenter(b);
    if (x < low || x > high) continue;
    if (best < 0 || smooth[b] > smooth[best]) best = b;
  }
  if (best < 0 || !(smooth[best] > 0)) return {std::nan(""), std::nan("")};
  const double half_max = 0.5 * smooth[best];
  int left = best, right = best;
  while (left > 1 && smooth[left] > half_max) --left;
  while (right < n && smooth[right] > half_max) ++right;
  return {h.GetBinCenter(best), h.GetBinCenter(right) - h.GetBinCenter(left)};
}

struct Window {
  int first = 0, last = 0;
  double lo = 0, width = 0;
  int bins() const { return last - first + 1; }
};

Window make_window(const TH1D& h, double low, double high) {
  Window w;
  w.first = h.GetXaxis()->FindFixBin(low + 1e-9);
  w.last = h.GetXaxis()->FindFixBin(high - 1e-9);
  w.first = std::max(w.first, 1);
  w.last = std::min(w.last, h.GetNbinsX());
  w.lo = h.GetXaxis()->GetBinLowEdge(w.first);
  w.width = h.GetXaxis()->GetBinWidth(w.first);
  return w;
}

h4l::FitHistogram window_data(const TH1D& h, const Window& w) {
  h4l::FitHistogram data;
  for (int b = w.first; b <= w.last; ++b) {
    data.content.push_back(h.GetBinContent(b));
    data.variance.push_back(h.GetBinError(b) * h.GetBinError(b));
  }
  return data;
}

double effective_entries(const h4l::FitHistogram& h) {
  double sw = 0, sw2 = 0;
  for (std::size_t b = 0; b < h.content.size(); ++b) {
    sw += h.content[b];
    sw2 += h.variance[b];
  }
  return sw2 > 0 ? sw * sw / sw2 : 0;
}

double total(const h4l::FitHistogram& h) {
  double sum = 0;
  for (double value : h.content) sum += value;
  return sum;
}

struct Tails {
  double alpha_l = 1.2, n_l = 3, alpha_r = 1.5, n_r = 5;
};

const char* kNames[] = {"delta", "sigma", "alpha_l", "n_l", "alpha_r", "n_r", "n_sig", "n_bkg", "bkg1", "bkg2"};

// One extended binned fit of one window; returns the result block.
json fit_one(const h4l::FitHistogram& data, const Window& w, const Settings& s, double delta0, double sigma0,
             const Tails& tails, bool fix_tails, Tails* fitted_tails) {
  const h4l::ZPeakModel model(w.lo, w.width, w.bins(), s.z_mass, s.z_width, s.background, s.grid_step, s.kernel_half);
  const double n = total(data);
  std::vector<h4l::FitParameter> p = {
      {"delta", std::clamp(delta0, -9.0, 9.0), 0.05, -10.0, 10.0, false},
      {"sigma", std::clamp(sigma0, 0.4, 7.0), 0.05, 0.2, 10.0, false},
      {"alpha_l", tails.alpha_l, 0.05, 0.2, 10.0, fix_tails},
      {"n_l", tails.n_l, 0.2, 1.01, 80.0, fix_tails},
      {"alpha_r", tails.alpha_r, 0.05, 0.2, 10.0, fix_tails},
      {"n_r", tails.n_r, 0.2, 1.01, 80.0, fix_tails},
      {"n_sig", 0.95 * n, 0.01 * n + 1.0, 0.0, 3.0 * n + 10.0, false},
      {"n_bkg", 0.05 * n, 0.01 * n + 1.0, 0.0, 3.0 * n + 10.0, model.background == h4l::ZPeakModel::Background::None},
      {"bkg1", model.background == h4l::ZPeakModel::Background::Exponential ? -0.03 : 0.0, 0.01, -1.0, 1.0,
       model.background == h4l::ZPeakModel::Background::None},
      {"bkg2", 0.0, 0.01, -1.0, 1.0, model.background != h4l::ZPeakModel::Background::Chebychev2}};
  if (model.background == h4l::ZPeakModel::Background::None) p[7].value = 0.0;
  const std::vector<h4l::FitHistogram> histograms = {data};
  const auto outcome = h4l::binned_fit(histograms, p, [&](const double* x, std::vector<std::vector<double>>& expected) {
    model.expected(x, expected[0]);
  });
  json out;
  out["weighted"] = outcome.weighted;
  out["window"] = {w.lo, w.lo + w.width * w.bins()};
  out["entries"] = n;
  out["effective_entries"] = effective_entries(data);
  out["status"] = outcome.status;
  out["cov_status"] = outcome.cov_status;
  out["strategy"] = outcome.strategy;
  out["ok"] = outcome.ok();
  out["edm"] = outcome.edm;
  out["nll"] = outcome.nll;
  for (std::size_t i = 0; i < p.size(); ++i)
    out[kNames[i]] = {{"value", outcome.values[i]}, {"error", outcome.errors[i]}, {"fixed", p[i].fixed}};
  out["peak"] = s.z_mass + outcome.values[0];
  const double e0 = outcome.errors[0], e1 = outcome.errors[1];
  out["corr_delta_sigma"] = e0 > 0 && e1 > 0 ? outcome.covariance[0][1] / (e0 * e1) : 0.0;
  // Pearson chi2 with the histogram variance, and the curves.
  std::vector<double> nu, nu_bkg;
  model.expected(outcome.values.data(), nu, &nu_bkg);
  double chi2 = 0;
  int bins_used = 0, floating = 0;
  for (const auto& parameter : p) floating += parameter.fixed ? 0 : 1;
  json x = json::array(), y = json::array(), ey = json::array(), f_total = json::array(), f_bkg = json::array();
  for (int b = 0; b < w.bins(); ++b) {
    if (data.variance[b] > 0) {
      chi2 += (data.content[b] - nu[b]) * (data.content[b] - nu[b]) / data.variance[b];
      ++bins_used;
    }
    if (s.curves) {
      x.push_back(w.lo + (b + 0.5) * w.width);
      y.push_back(data.content[b]);
      ey.push_back(std::sqrt(data.variance[b]));
      f_total.push_back(nu[b]);
      f_bkg.push_back(nu_bkg[b]);
    }
  }
  out["chi2"] = chi2;
  out["ndf"] = bins_used - floating;
  if (s.curves) out["curve"] = {{"x", x}, {"y", y}, {"ey", ey}, {"model", f_total}, {"background", f_bkg}};
  if (fitted_tails) *fitted_tails = Tails{outcome.values[2], outcome.values[3], outcome.values[4], outcome.values[5]};
  return out;
}
}  // namespace

int main(int argc, char** argv) {
  std::string job_path, out_path;
  for (int index = 1; index + 1 < argc; index += 2) {
    const std::string key = argv[index];
    if (key == "--job") job_path = argv[index + 1];
    else if (key == "--out") out_path = argv[index + 1];
  }
  if (job_path.empty() || out_path.empty()) {
    std::cerr << "usage: fit_zpeaks --job JOB.json --out OUT.json\n";
    return 64;
  }
  gROOT->SetBatch(true);
  try {
    const json job = h4l::read_json(job_path);
    const json& fit = job.at("calibration_config").at("validation_fit");
    const json variant = job.value("variant", json::object());
    Settings s;
    s.mode_low = fit.at("mode_search").at(0).get<double>();
    s.mode_high = fit.at("mode_search").at(1).get<double>();
    s.smooth_bins = fit.at("smooth_bins").get<int>();
    s.below = variant.value("below", fit.at("below").get<double>());
    s.above = variant.value("above", fit.at("above").get<double>());
    s.clip_low = fit.at("window_clip").at(0).get<double>();
    s.clip_high = fit.at("window_clip").at(1).get<double>();
    s.min_data = fit.at("min_data_events").get<double>();
    s.min_mc_effective = fit.at("min_mc_effective").get<double>();
    s.z_mass = fit.at("z_mass").get<double>();
    s.z_width = fit.at("z_width").get<double>();
    s.grid_step = fit.at("grid_step").get<double>();
    s.kernel_half = fit.at("kernel_half_width").get<double>();
    s.background = variant.value("background", fit.at("background").get<std::string>());
    s.data_tails = variant.value("data_tails", fit.at("data_tails").get<std::string>());
    s.curves = job.value("curves", true);
    if (s.data_tails != "fixed_from_mc" && s.data_tails != "free") throw std::runtime_error("unknown data_tails " + s.data_tails);
    h4l::ZPeakModel::parse_background(s.background);

    SetErrorHandler(record_errors);
    std::vector<std::unique_ptr<TFile>> data_files, mc_files;
    for (const auto& [role, files] : {std::pair{"data", &data_files}, std::pair{"mc", &mc_files}})
      for (const auto& path : job.at(role)) {
        std::unique_ptr<TFile> file(TFile::Open(path.get<std::string>().c_str(), "READ"));
        if (!file || file->IsZombie() || file->TestBit(TFile::kRecovered))
          throw std::runtime_error("cannot open " + path.get<std::string>());
        files->push_back(std::move(file));
      }
    json results = json::array();
    for (const auto& category_node : job.at("categories")) {
      const std::string category = category_node.get<std::string>();
      g_read_error = false;
      auto data = summed(data_files, category), mc = summed(mc_files, category);
      if (g_read_error) throw std::runtime_error("ROOT error while reading " + category + ": " + g_read_message);
      json entry = {{"category", category}, {"ok", false}};
      if (!data || !mc) {
        entry["skipped"] = !data ? "no data histogram" : "no MC histogram";
        results.push_back(entry);
        continue;
      }
      const auto [mode, fwhm] = smoothed_mode(*data, s.mode_low, s.mode_high, s.smooth_bins);
      entry["data_mode"] = mode;
      entry["data_fwhm"] = fwhm;
      if (!std::isfinite(mode)) {
        entry["skipped"] = "no data peak";
        results.push_back(entry);
        continue;
      }
      // The MC role may be at the finer template binning: rebin exactly.
      if (mc->GetNbinsX() != data->GetNbinsX()) {
        const int ratio = mc->GetNbinsX() / data->GetNbinsX();
        if (ratio < 1 || ratio * data->GetNbinsX() != mc->GetNbinsX() ||
            std::fabs(mc->GetXaxis()->GetXmin() - data->GetXaxis()->GetXmin()) > 1e-9 ||
            std::fabs(mc->GetXaxis()->GetXmax() - data->GetXaxis()->GetXmax()) > 1e-9)
          throw std::runtime_error(category + ": MC binning is not an integer refinement of the data binning");
        mc->Rebin(ratio);
      }
      const Window w = make_window(*data, std::max(s.clip_low, mode - s.below), std::min(s.clip_high, mode + s.above));
      const auto data_window = window_data(*data, w), mc_window = window_data(*mc, w);
      entry["data_entries"] = total(data_window);
      entry["data_effective_entries"] = effective_entries(data_window);
      entry["mc_effective_entries"] = effective_entries(mc_window);
      if (effective_entries(data_window) < s.min_data || effective_entries(mc_window) < s.min_mc_effective) {
        entry["skipped"] = "low statistics";
        results.push_back(entry);
        continue;
      }
      // Initial resolution from the half-maximum width of a Voigt profile.
      const double g = s.z_width;
      const double core = std::pow(fwhm - 0.5346 * g, 2) - 0.2166 * g * g;
      const double sigma0 = core > 0 && fwhm > 0.5346 * g ? std::sqrt(core) / 2.3548 : 1.0;
      const auto [mc_mode, mc_fwhm] = smoothed_mode(*mc, s.mode_low, s.mode_high, s.smooth_bins);
      (void)mc_fwhm;
      Tails tails;
      entry["mc"] = fit_one(mc_window, w, s, (std::isfinite(mc_mode) ? mc_mode : mode) - s.z_mass, sigma0, Tails{}, false, &tails);
      entry["data"] = fit_one(data_window, w, s, mode - s.z_mass, sigma0, tails, s.data_tails == "fixed_from_mc", nullptr);
      entry["ok"] = entry["mc"]["ok"].get<bool>() && entry["data"]["ok"].get<bool>();
      std::cout << "[fit] " << category << ": data " << entry["data"]["peak"].get<double>() << " +- "
                << entry["data"]["delta"]["error"].get<double>() << " sigma " << entry["data"]["sigma"]["value"].get<double>()
                << ", MC " << entry["mc"]["peak"].get<double>() << " +- " << entry["mc"]["delta"]["error"].get<double>()
                << " sigma " << entry["mc"]["sigma"]["value"].get<double>() << (entry["ok"].get<bool>() ? "" : "  [FAILED]")
                << "\n";
      results.push_back(entry);
    }
    json report = {{"schema", "h4l_v3_zpeak_fits/1"},
                   {"flavour", job.at("flavour")},
                   {"variant", {{"below", s.below}, {"above", s.above}, {"background", s.background}, {"data_tails", s.data_tails}}},
                   {"data", job.at("data")},
                   {"mc", job.at("mc")},
                   {"results", results},
                   {"finished_utc", h4l::utc_now()}};
    h4l::publish_json(out_path, report);
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
