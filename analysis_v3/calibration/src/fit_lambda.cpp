// fit_lambda: the per-event mass-uncertainty correction factors lambda of
// AN-16-442 5.3.1 from Z -> ll pairs binned in the predicted uncertainty.
//
// Usage: fit_lambda --job JOB.json --out OUT.json
//
// The job names the lambda_histograms outputs of the data role and of the MC
// role (each role is the sum of its files), the flavour tag ("mm" or "ee")
// and the lambda configuration (lambda_ul16_v2.json).  Every configured fit
// of the flavour, in order, is done for the MC role and then for the data
// role:
//   * the pair class is (target, target) for mode "same" and (reference,
//     target) for mode "reference"; its categories are the e bins with at
//     least min_category_entries pairs in the window;
//   * all categories are fitted simultaneously (extended binned Poisson
//     likelihood, h4l/binned_fit.h) with BW(m_Z, Gamma_Z) (x) DCB plus an
//     exponential (h4l/zpeak_model.h): the DCB tails (alphaL, nL, alphaR,
//     nR) are shared, every category has its own shift, signal and
//     background yields and slope, and the DCB sigma of category k is
//       sigma_k = (m_Z / 2) sqrt(lambda_a^2 <d_a^2>_k + lambda_b^2 <d_b^2>_k),
//     with <d^2> the mean squared relative momentum error of the legs in
//     80-100 GeV; the target lambda is the free parameter, a reference
//     leg's lambda is fixed to its own fit of the same role; the data fit
//     takes the MC tails when data_tails = "fixed_from_mc";
//   * closure (AN Figure 18): every category is refitted alone with sigma
//     free and the tails of the simultaneous fit, giving the measured sigma
//     against the predicted sigma before (lambda = 1) and after the
//     correction.
// The report lists every fit's status, parameters, errors and per-category
// curves.

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
#include <map>
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

std::unique_ptr<TH1D> summed(const std::vector<std::unique_ptr<TFile>>& files, const std::string& name) {
  std::unique_ptr<TH1D> total;
  for (const auto& file : files) {
    auto* h = dynamic_cast<TH1D*>(file->Get(name.c_str()));
    if (!h) continue;
    if (!total) {
      total.reset(static_cast<TH1D*>(h->Clone((name + "_sum").c_str())));
      total->SetDirectory(nullptr);
    } else {
      total->Add(h);
    }
  }
  return total;
}

struct Category {
  int bin = 0;
  h4l::FitHistogram hist;
  double entries = 0, mean_e = 0, mean_da2 = 0, mean_db2 = 0;
};

struct Settings {
  double lo = 60, hi = 120, z_mass = 91.1876, z_width = 2.4952, grid_step = 0.125, kernel_half = 25, min_entries = 150;
  std::string background = "exponential", data_tails = "fixed_from_mc";
};

struct Tails {
  double alpha_l = 1.5, n_l = 3, alpha_r = 1.8, n_r = 5;
};

// Window of the histograms: bins whose centres lie in [lo, hi).
h4l::FitHistogram window_of(const TH1D& h, const Settings& s, double& low_edge, double& width, int& bins) {
  h4l::FitHistogram out;
  width = h.GetXaxis()->GetBinWidth(1);
  low_edge = 0;
  bins = 0;
  for (int b = 1; b <= h.GetNbinsX(); ++b) {
    const double centre = h.GetXaxis()->GetBinCenter(b);
    if (centre < s.lo || centre >= s.hi) continue;
    if (bins == 0) low_edge = h.GetXaxis()->GetBinLowEdge(b);
    out.content.push_back(h.GetBinContent(b));
    out.variance.push_back(h.GetBinError(b) * h.GetBinError(b));
    ++bins;
  }
  return out;
}

double total(const h4l::FitHistogram& h) {
  double sum = 0;
  for (double v : h.content) sum += v;
  return sum;
}

struct ClassFit {
  bool ok = false;
  double lambda = 1, lambda_error = 0;
  Tails tails;
  json report;
};

// One simultaneous fit of the categories of a pair class.
ClassFit fit_class(const std::vector<Category>& cats, const h4l::ZPeakModel& model, const Settings& s, bool target_is_a,
                   bool same, double lambda_ref, const Tails* fixed_tails, bool curves) {
  const int n_cat = static_cast<int>(cats.size());
  std::vector<h4l::FitParameter> pars = {
      {"lambda", 1.3, 0.02, 0.3, 5.0, false},
      {"alpha_l", fixed_tails ? fixed_tails->alpha_l : 1.5, 0.05, 0.2, 10.0, fixed_tails != nullptr},
      {"n_l", fixed_tails ? fixed_tails->n_l : 3.0, 0.2, 1.01, 80.0, fixed_tails != nullptr},
      {"alpha_r", fixed_tails ? fixed_tails->alpha_r : 1.8, 0.05, 0.2, 10.0, fixed_tails != nullptr},
      {"n_r", fixed_tails ? fixed_tails->n_r : 5.0, 0.2, 1.01, 80.0, fixed_tails != nullptr}};
  for (int k = 0; k < n_cat; ++k) {
    const double n = std::max(total(cats[k].hist), 1.0);
    const std::string t = "_" + std::to_string(cats[k].bin);
    pars.push_back({"delta" + t, -0.2, 0.02, -5.0, 5.0, false});
    pars.push_back({"n_sig" + t, 0.97 * n, 0.01 * n + 1.0, 0.0, 3.0 * n + 10.0, false});
    pars.push_back({"n_bkg" + t, 0.03 * n, 0.01 * n + 1.0, 0.0, 3.0 * n + 10.0, false});
    pars.push_back({"tau" + t, -0.03, 0.01, -1.0, 1.0, false});
  }
  auto sigma_of = [&](const double* p, const Category& c) {
    const double la = same ? p[0] : (target_is_a ? p[0] : lambda_ref);
    const double lb = same ? p[0] : (target_is_a ? lambda_ref : p[0]);
    return 0.5 * s.z_mass * std::sqrt(la * la * c.mean_da2 + lb * lb * c.mean_db2);
  };
  auto category_parameters = [&](const double* p, int k, double out[10]) {
    const double* q = p + 5 + 4 * k;
    out[0] = q[0];
    out[1] = sigma_of(p, cats[k]);
    out[2] = p[1];
    out[3] = p[2];
    out[4] = p[3];
    out[5] = p[4];
    out[6] = q[1];
    out[7] = q[2];
    out[8] = q[3];
    out[9] = 0.0;
  };
  std::vector<h4l::FitHistogram> data;
  for (const auto& c : cats) data.push_back(c.hist);
  h4l::FitModel fit_model = [&](const double* p, std::vector<std::vector<double>>& expected) {
    double q[10];
    for (int k = 0; k < n_cat; ++k) {
      category_parameters(p, k, q);
      model.expected(q, expected[k]);
    }
  };
  // One attempt from the given starting values.  A free nuisance parameter that ends at (or next to) a limit
  // leaves MIGRAD/HESSE unreliable (a DCB tail exponent near its upper limit is degenerate with its alpha; a
  // background yield at zero): fix every such parameter at its fitted value and refit, up to three rounds (the
  // staged freezing of the tag-and-probe fits).  "At a limit": within 1 % of the range for shapes, within half
  // an event for the background yields; lambda is never fixed.
  auto attempt = [&](std::vector<h4l::FitParameter> start, json& fixed) {
    h4l::FitOutcome o = h4l::binned_fit(data, start, fit_model);
    for (int round = 0; round < 3 && !o.ok(); ++round) {
      bool changed = false;
      for (std::size_t i = 1; i < start.size(); ++i) {
        if (start[i].fixed) continue;
        const double span = start[i].high - start[i].low, v = o.values[i];
        const bool yield = start[i].name.rfind("n_bkg", 0) == 0 || start[i].name.rfind("n_sig", 0) == 0;
        const double tolerance = yield ? 0.5 : 1e-2 * span;
        if (v - start[i].low < tolerance || start[i].high - v < tolerance) {
          start[i].value = v;
          start[i].fixed = true;
          fixed.push_back(start[i].name);
          changed = true;
        }
      }
      if (!changed) break;
      for (std::size_t i = 0; i < start.size(); ++i)
        if (!start[i].fixed) start[i].value = std::clamp(o.values[i], start[i].low, start[i].high);
      o = h4l::binned_fit(data, start, fit_model);
    }
    return std::pair{o, start};
  };
  json fixed_at_limit = json::array();
  auto [outcome, final_pars] = attempt(pars, fixed_at_limit);
  // Free tails that did not converge from the default start: alternative starting points of the shared tails,
  // the converged attempt with the lowest NLL kept.
  int tail_start = 0;
  if (!outcome.ok() && !fixed_tails) {
    const double starts[][4] = {{1.0, 5.0, 1.5, 10.0}, {2.0, 2.0, 2.5, 3.0}, {1.2, 10.0, 3.0, 20.0}};
    for (int k = 0; k < 3; ++k) {
      std::vector<h4l::FitParameter> alt = pars;
      for (int i = 0; i < 4; ++i) alt[1 + i].value = starts[k][i];
      json fixed = json::array();
      auto [o, q] = attempt(alt, fixed);
      if (o.ok() && (!outcome.ok() || o.nll < outcome.nll - 1e-6)) {
        outcome = o;
        final_pars = q;
        fixed_at_limit = fixed;
        tail_start = k + 1;
      }
    }
  }
  pars = final_pars;
  ClassFit result;
  result.ok = outcome.ok();
  result.lambda = outcome.values[0];
  result.lambda_error = outcome.errors[0];
  result.tails = {outcome.values[1], outcome.values[2], outcome.values[3], outcome.values[4]};
  json& r = result.report;
  r["status"] = outcome.status;
  r["cov_status"] = outcome.cov_status;
  r["strategy"] = outcome.strategy;
  r["ok"] = outcome.ok();
  r["weighted"] = outcome.weighted;
  r["nll"] = outcome.nll;
  r["edm"] = outcome.edm;
  r["fixed_at_limit"] = fixed_at_limit;
  r["tail_start"] = tail_start;
  r["lambda"] = {{"value", outcome.values[0]}, {"error", outcome.errors[0]}};
  r["tails"] = {{"alpha_l", {{"value", outcome.values[1]}, {"error", outcome.errors[1]}, {"fixed", pars[1].fixed}}},
                {"n_l", {{"value", outcome.values[2]}, {"error", outcome.errors[2]}, {"fixed", pars[2].fixed}}},
                {"alpha_r", {{"value", outcome.values[3]}, {"error", outcome.errors[3]}, {"fixed", pars[3].fixed}}},
                {"n_r", {{"value", outcome.values[4]}, {"error", outcome.errors[4]}, {"fixed", pars[4].fixed}}}};
  double chi2 = 0;
  int bins_used = 0;
  json categories = json::array();
  std::vector<std::vector<double>> expected(n_cat);
  fit_model(outcome.values.data(), expected);
  for (int k = 0; k < n_cat; ++k) {
    const Category& c = cats[k];
    double q[10];
    category_parameters(outcome.values.data(), k, q);
    const double predicted_before = 0.5 * s.z_mass * std::sqrt(c.mean_da2 + c.mean_db2);
    json cat = {{"bin", c.bin}, {"entries", c.entries}, {"mean_e", c.mean_e}, {"mean_da2", c.mean_da2}, {"mean_db2", c.mean_db2},
                {"predicted_sigma_before", predicted_before}, {"predicted_sigma_after", q[1]},
                {"delta", {{"value", outcome.values[5 + 4 * k]}, {"error", outcome.errors[5 + 4 * k]}}},
                {"n_sig", outcome.values[6 + 4 * k]}, {"n_bkg", outcome.values[7 + 4 * k]}, {"tau", outcome.values[8 + 4 * k]}};
    double chi2_k = 0;
    int bins_k = 0;
    std::vector<double> x, y, ey, f_total, f_bkg;
    std::vector<double> nu_bkg;
    model.expected(q, expected[k], &nu_bkg);
    // Pearson chi2 with the variance floored at the expectation times the mean weight, so that
    // sparsely populated bins of mixed-weight MC do not dominate.
    double sum_n = 0, sum_v = 0;
    for (std::size_t b = 0; b < c.hist.content.size(); ++b) {
      sum_n += c.hist.content[b];
      sum_v += c.hist.variance[b];
    }
    const double mean_weight = sum_n > 0 ? sum_v / sum_n : 1.0;
    for (std::size_t b = 0; b < c.hist.content.size(); ++b) {
      const double v = c.hist.variance[b];
      const double variance = std::max(v, expected[k][b] * mean_weight);
      if (variance > 0) {
        chi2_k += std::pow(c.hist.content[b] - expected[k][b], 2) / variance;
        ++bins_k;
      }
      if (curves) {
        x.push_back(model.lo + (b + 0.5) * model.bin_width);
        y.push_back(c.hist.content[b]);
        ey.push_back(std::sqrt(std::max(v, 0.0)));
        f_total.push_back(expected[k][b]);
        f_bkg.push_back(nu_bkg[b]);
      }
    }
    chi2 += chi2_k;
    bins_used += bins_k;
    cat["chi2"] = chi2_k;
    cat["bins"] = bins_k;
    if (curves) cat["curve"] = {{"x", x}, {"y", y}, {"ey", ey}, {"model", f_total}, {"background", f_bkg}};
    categories.push_back(cat);
  }
  int floating = 0;
  for (const auto& p : pars) floating += p.fixed ? 0 : 1;
  r["chi2"] = chi2;
  r["ndf"] = bins_used - floating;
  r["categories"] = categories;
  return result;
}

// Closure: one category alone, sigma free, tails fixed.
json closure_fit(const Category& c, const h4l::ZPeakModel& model, const Tails& tails, double sigma0) {
  const double n = std::max(total(c.hist), 1.0);
  std::vector<h4l::FitParameter> pars = {
      {"delta", -0.2, 0.02, -5.0, 5.0, false},        {"sigma", std::clamp(sigma0, 0.3, 6.0), 0.02, 0.1, 10.0, false},
      {"alpha_l", tails.alpha_l, 0.05, 0.2, 10.0, true}, {"n_l", tails.n_l, 0.2, 1.01, 80.0, true},
      {"alpha_r", tails.alpha_r, 0.05, 0.2, 10.0, true}, {"n_r", tails.n_r, 0.2, 1.01, 80.0, true},
      {"n_sig", 0.97 * n, 0.01 * n + 1.0, 0.0, 3.0 * n + 10.0, false},
      {"n_bkg", 0.03 * n, 0.01 * n + 1.0, 0.0, 3.0 * n + 10.0, false},
      {"tau", -0.03, 0.01, -1.0, 1.0, false},         {"unused", 0.0, 0.01, -1.0, 1.0, true}};
  std::vector<h4l::FitHistogram> data = {c.hist};
  h4l::FitModel fit_model = [&](const double* p, std::vector<std::vector<double>>& expected) {
    model.expected(p, expected[0]);
  };
  const h4l::FitOutcome outcome = h4l::binned_fit(data, pars, fit_model);
  return {{"ok", outcome.ok()}, {"sigma", outcome.values[1]}, {"sigma_error", outcome.errors[1]},
          {"delta", outcome.values[0]}, {"status", outcome.status}, {"cov_status", outcome.cov_status}};
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
    std::cerr << "usage: fit_lambda --job JOB.json --out OUT.json\n";
    return 64;
  }
  gROOT->SetBatch(true);
  TH1::AddDirectory(false);
  try {
    const json job = h4l::read_json(job_path);
    const json& lambda = job.at("lambda_config");
    const std::string tag = job.at("flavour").get<std::string>();
    const json& flavour = lambda.at(tag == "mm" ? "muon" : "electron");
    const json& fit = lambda.at("fit");
    Settings s;
    s.lo = fit.at("window").at(0).get<double>();
    s.hi = fit.at("window").at(1).get<double>();
    s.z_mass = fit.at("z_mass").get<double>();
    s.z_width = fit.at("z_width").get<double>();
    s.grid_step = fit.at("grid_step").get<double>();
    s.kernel_half = fit.at("kernel_half_width").get<double>();
    s.background = fit.at("background").get<std::string>();
    s.min_entries = fit.at("min_category_entries").get<double>();
    s.data_tails = fit.at("data_tails").get<std::string>();
    if (s.data_tails != "fixed_from_mc" && s.data_tails != "free") throw std::runtime_error("unknown data_tails " + s.data_tails);
    const bool curves = job.value("curves", true);
    const int n_bins = static_cast<int>(flavour.at("pair_error_edges").size()) + 1;

    SetErrorHandler(record_errors);
    std::map<std::string, std::vector<std::unique_ptr<TFile>>> files;
    for (const std::string role : {"data", "mc"})
      for (const auto& path : job.at(role)) {
        std::unique_ptr<TFile> file(TFile::Open(path.get<std::string>().c_str(), "READ"));
        if (!file || file->IsZombie() || file->TestBit(TFile::kRecovered))
          throw std::runtime_error("cannot open " + path.get<std::string>());
        files[role].push_back(std::move(file));
      }
    std::unique_ptr<h4l::ZPeakModel> model;
    std::map<std::pair<std::string, int>, double> fitted_lambda;  // (role, region) -> lambda
    json fits_report = json::array();
    for (const auto& spec : flavour.at("fits")) {
      const int target = spec.at("target").get<int>();
      const std::string mode = spec.at("mode").get<std::string>();
      const bool same = mode == "same";
      if (!same && mode != "reference") throw std::runtime_error("unknown lambda fit mode " + mode);
      const int reference = same ? target : spec.at("reference").get<int>();
      if (!same && reference == target) throw std::runtime_error("a reference fit needs another region");
      const int a = std::min(target, reference), b = std::max(target, reference);
      json entry = {{"target", target}, {"mode", mode}, {"class", {a, b}}};
      if (!same) entry["reference"] = reference;
      Tails mc_tails;
      bool have_mc_tails = false;
      for (const std::string role : {"mc", "data"}) {
        double lambda_ref = 1.0;
        if (!same) {
          const auto found = fitted_lambda.find({role, reference});
          if (found == fitted_lambda.end())
            throw std::runtime_error("reference region " + std::to_string(reference) + " has no lambda for " + role);
          lambda_ref = found->second;
        }
        std::vector<Category> cats;
        for (int k = 0; k < n_bins; ++k) {
          const std::string suffix = tag + "_" + std::to_string(a) + "_" + std::to_string(b) + "_" + std::to_string(k);
          g_read_error = false;
          auto h = summed(files[role], "lm_" + suffix);
          auto stats = summed(files[role], "le_" + suffix);
          if (g_read_error) throw std::runtime_error("ROOT error while reading " + suffix + ": " + g_read_message);
          if (!h || !stats || !(stats->GetBinContent(1) > 0)) continue;
          Category c;
          c.bin = k;
          double low = 0, width = 0;
          int bins = 0;
          c.hist = window_of(*h, s, low, width, bins);
          c.entries = total(c.hist);
          if (c.entries < s.min_entries) continue;
          const double w = stats->GetBinContent(1);
          c.mean_e = stats->GetBinContent(2) / w;
          c.mean_da2 = stats->GetBinContent(4) / w;
          c.mean_db2 = stats->GetBinContent(5) / w;
          if (!model) {
            model = std::make_unique<h4l::ZPeakModel>(low, width, bins, s.z_mass, s.z_width, s.background, s.grid_step,
                                                      s.kernel_half);
          } else if (std::fabs(model->lo - low) > 1e-9 || std::fabs(model->bin_width - width) > 1e-9 || model->bins != bins) {
            throw std::runtime_error("the lambda histograms do not share one binning");
          }
          cats.push_back(std::move(c));
        }
        json role_report = {{"categories_used", static_cast<int>(cats.size())}, {"lambda_reference", lambda_ref}};
        if (cats.size() < 2) {
          role_report["ok"] = false;
          role_report["skipped"] = "fewer than two categories above min_category_entries";
          entry[role] = role_report;
          std::cout << "[lambda] " << tag << " region " << target << " " << role << ": skipped (" << cats.size()
                    << " categories)\n";
          continue;
        }
        const bool fix_tails = role == "data" && s.data_tails == "fixed_from_mc";
        if (fix_tails && !have_mc_tails) {
          role_report["ok"] = false;
          role_report["skipped"] = "no MC tails";
          entry[role] = role_report;
          continue;
        }
        const bool target_is_a = target == a;
        ClassFit result = fit_class(cats, *model, s, target_is_a, same, lambda_ref, fix_tails ? &mc_tails : nullptr, curves);
        role_report.update(result.report);
        if (role == "mc" && result.ok) {
          mc_tails = result.tails;
          have_mc_tails = true;
        }
        if (result.ok) fitted_lambda[{role, target}] = result.lambda;
        // Closure fits per category with the tails of this fit.
        for (std::size_t k = 0; k < cats.size(); ++k) {
          const double sigma0 = role_report["categories"][k]["predicted_sigma_after"].get<double>();
          role_report["categories"][k]["closure"] = closure_fit(cats[k], *model, result.tails, sigma0);
        }
        std::cout << "[lambda] " << tag << " region " << target << " (" << mode << ") " << role << ": lambda "
                  << result.lambda << " +- " << result.lambda_error << ", " << cats.size() << " categories, chi2 "
                  << role_report["chi2"].get<double>() << "/" << role_report["ndf"].get<int>()
                  << (result.ok ? "" : "  [FAILED]") << "\n";
        entry[role] = role_report;
      }
      fits_report.push_back(entry);
    }
    json report = {{"schema", "h4l_v3_lambda_fits/1"},
                   {"flavour", tag},
                   {"regions", flavour.at("regions")},
                   {"pair_error_edges", flavour.at("pair_error_edges")},
                   {"settings", fit},
                   {"data", job.at("data")},
                   {"mc", job.at("mc")},
                   {"fits", fits_report},
                   {"finished_utc", h4l::utc_now()}};
    h4l::publish_json(out_path, report);
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
