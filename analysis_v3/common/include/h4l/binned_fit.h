// Extended binned Poisson likelihood fits with Minuit2.
//
// A model maps the parameter vector to the expected content of every bin of
// one or more histograms; the fit minimizes
//   NLL = sum_h sum_b [nu_hb - n_hb + n_hb ln(n_hb / nu_hb)]
// (the Poisson deviance / 2, zero for a saturated model) with MIGRAD, then
// runs HESSE.  For weighted histograms (n_hb = sum of weights, variance
// v_hb = sum of squared weights) the parameter covariance is the sandwich
// estimator C = V G V with V the HESSE covariance and
// G_ij = sum_hb v_hb (d nu_hb / d theta_i)(d nu_hb / d theta_j) / nu_hb^2,
// which reduces to V for unweighted histograms (v = n).  With minos_index
// >= 0 the MINOS (profile-likelihood) interval of that parameter is computed
// after a successful HESSE; for weighted histograms it is scaled by the ratio
// of the sandwich to the HESSE error of the parameter.  WeightedErrors::
// Effective instead scales the HESSE covariance and the MINOS interval by
// the effective-weight factor s^2 = sum v / sum n over every bin (exact for
// uniform weights; the numerical sandwich is unstable when nuisance
// parameters are weakly constrained or at their bounds).
#pragma once

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>

#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace h4l {

struct FitParameter {
  std::string name;
  double value = 0, step = 0.1, low = -1e30, high = 1e30;
  bool fixed = false;
};

struct FitHistogram {
  std::vector<double> content, variance;  // per bin
};

using FitModel = std::function<void(const double* parameters, std::vector<std::vector<double>>& expected)>;

enum class WeightedErrors { Sandwich, Effective };

struct FitOutcome {
  int status = -1, cov_status = -1, strategy = -1;
  double nll = 0, edm = 0;
  bool weighted = false;
  std::vector<double> values, errors;
  std::vector<std::vector<double>> covariance;  // full matrix over all parameters (zero rows for fixed ones)
  bool minos_ok = false;
  double minos_low = 0, minos_high = 0;  // signed interval ends (low <= 0 <= high)
  bool ok() const { return status == 0 && cov_status >= 2; }
};

inline bool histogram_weighted(const FitHistogram& h) {
  for (std::size_t b = 0; b < h.content.size(); ++b) {
    const double n = h.content[b];
    if (std::fabs(n - std::round(n)) > 1e-6 || std::fabs(h.variance[b] - n) > 1e-6 * std::max(1.0, n)) return true;
  }
  return false;
}

inline double poisson_deviance_half(const std::vector<FitHistogram>& data, const std::vector<std::vector<double>>& expected) {
  double nll = 0;
  for (std::size_t h = 0; h < data.size(); ++h)
    for (std::size_t b = 0; b < data[h].content.size(); ++b) {
      const double nu = expected[h][b], n = data[h].content[b];
      if (!(nu > 0) || !std::isfinite(nu)) {
        if (nu == 0 && n == 0) continue;
        return 1e30;
      }
      nll += nu - n;
      if (n > 0) nll += n * std::log(n / nu);
      else if (n < 0) nll -= n * std::log(nu);  // negative-weight bins keep the likelihood defined
    }
  return nll;
}

// Minimize with strategies 1 then 2 until MIGRAD and HESSE succeed.
inline FitOutcome binned_fit(const std::vector<FitHistogram>& data, const std::vector<FitParameter>& parameters,
                             const FitModel& model, int minos_index = -1,
                             WeightedErrors weighted_errors = WeightedErrors::Sandwich) {
  const std::size_t n_par = parameters.size();
  std::vector<std::vector<double>> expected(data.size());
  for (std::size_t h = 0; h < data.size(); ++h) expected[h].assign(data[h].content.size(), 0.0);
  auto nll = [&](const double* p) {
    model(p, expected);
    return poisson_deviance_half(data, expected);
  };
  FitOutcome outcome;
  outcome.weighted = false;
  for (const auto& h : data) outcome.weighted = outcome.weighted || histogram_weighted(h);
  for (int strategy : {1, 2}) {
    std::unique_ptr<ROOT::Math::Minimizer> minimizer(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    if (!minimizer) throw std::runtime_error("cannot create the Minuit2 minimizer");
    ROOT::Math::Functor functor(nll, static_cast<unsigned>(n_par));
    minimizer->SetFunction(functor);
    minimizer->SetErrorDef(0.5);
    minimizer->SetStrategy(strategy);
    minimizer->SetMaxFunctionCalls(200000);
    minimizer->SetMaxIterations(200000);
    minimizer->SetTolerance(0.01);
    minimizer->SetPrintLevel(-1);
    for (std::size_t i = 0; i < n_par; ++i) {
      const auto& p = parameters[i];
      if (p.fixed) minimizer->SetFixedVariable(static_cast<unsigned>(i), p.name, p.value);
      else minimizer->SetLimitedVariable(static_cast<unsigned>(i), p.name, p.value, p.step, p.low, p.high);
    }
    minimizer->Minimize();
    minimizer->Hesse();
    outcome.status = minimizer->Status();
    outcome.cov_status = minimizer->CovMatrixStatus();
    outcome.strategy = strategy;
    outcome.edm = minimizer->Edm();
    outcome.nll = minimizer->MinValue();
    const double* x = minimizer->X();
    outcome.values.assign(x, x + n_par);
    outcome.covariance.assign(n_par, std::vector<double>(n_par, 0.0));
    for (std::size_t i = 0; i < n_par; ++i)
      for (std::size_t j = 0; j < n_par; ++j)
        if (!parameters[i].fixed && !parameters[j].fixed)
          outcome.covariance[i][j] = minimizer->CovMatrix(static_cast<unsigned>(i), static_cast<unsigned>(j));
    const bool converged = outcome.status == 0 || (outcome.status == 3 && outcome.edm < 0.01);
    if (minos_index >= 0 && static_cast<std::size_t>(minos_index) < n_par && !parameters[minos_index].fixed && converged) {
      double low = 0, high = 0;
      outcome.minos_ok = minimizer->GetMinosError(static_cast<unsigned>(minos_index), low, high) && std::isfinite(low) &&
                         std::isfinite(high) && low <= 0 && high >= 0 && (high - low) > 0;
      outcome.minos_low = low;
      outcome.minos_high = high;
      // MINOS may move the minimum: keep the values of the best point.
      const double* xm = minimizer->X();
      if (minimizer->MinValue() < outcome.nll - 1e-6) {
        outcome.values.assign(xm, xm + n_par);
        outcome.nll = minimizer->MinValue();
      }
    }
    if (outcome.ok()) break;
  }
  const double hesse_error_before_sandwich =
      minos_index >= 0 && static_cast<std::size_t>(minos_index) < n_par && !outcome.covariance.empty()
          ? std::sqrt(std::max(outcome.covariance[minos_index][minos_index], 0.0))
          : 0.0;
  if (outcome.weighted && weighted_errors == WeightedErrors::Effective) {
    double sum_n = 0, sum_v = 0;
    for (const auto& h : data)
      for (std::size_t b = 0; b < h.content.size(); ++b) {
        sum_n += h.content[b];
        sum_v += h.variance[b];
      }
    const double s2 = sum_n > 0 ? sum_v / sum_n : 1.0;
    for (auto& row : outcome.covariance)
      for (double& c : row) c *= s2;
    outcome.errors.assign(n_par, 0.0);
    for (std::size_t i = 0; i < n_par; ++i) outcome.errors[i] = std::sqrt(std::max(outcome.covariance[i][i], 0.0));
    outcome.minos_low *= std::sqrt(s2);
    outcome.minos_high *= std::sqrt(s2);
    return outcome;
  }
  if (outcome.weighted && outcome.ok()) {
    // Sandwich covariance: Jacobian of the expected contents by central differences.
    std::vector<std::size_t> free;
    for (std::size_t i = 0; i < n_par; ++i)
      if (!parameters[i].fixed) free.push_back(i);
    std::vector<std::vector<std::vector<double>>> jacobian(free.size());
    std::vector<double> p = outcome.values;
    for (std::size_t f = 0; f < free.size(); ++f) {
      const std::size_t i = free[f];
      const double error = std::sqrt(std::max(outcome.covariance[i][i], 0.0));
      const double step = error > 0 ? 0.01 * error : 1e-6 * std::max(1.0, std::fabs(p[i]));
      std::vector<std::vector<double>> up = expected, down = expected;
      p[i] = outcome.values[i] + step;
      model(p.data(), up);
      p[i] = outcome.values[i] - step;
      model(p.data(), down);
      p[i] = outcome.values[i];
      jacobian[f].resize(data.size());
      for (std::size_t h = 0; h < data.size(); ++h) {
        jacobian[f][h].resize(data[h].content.size());
        for (std::size_t b = 0; b < data[h].content.size(); ++b) jacobian[f][h][b] = (up[h][b] - down[h][b]) / (2 * step);
      }
    }
    model(outcome.values.data(), expected);
    std::vector<std::vector<double>> g(free.size(), std::vector<double>(free.size(), 0.0));
    for (std::size_t h = 0; h < data.size(); ++h)
      for (std::size_t b = 0; b < data[h].content.size(); ++b) {
        const double nu = expected[h][b];
        if (!(nu > 0)) continue;
        const double factor = data[h].variance[b] / (nu * nu);
        for (std::size_t a = 0; a < free.size(); ++a)
          for (std::size_t c = 0; c < free.size(); ++c) g[a][c] += factor * jacobian[a][h][b] * jacobian[c][h][b];
      }
    std::vector<std::vector<double>> corrected(n_par, std::vector<double>(n_par, 0.0));
    for (std::size_t a = 0; a < free.size(); ++a)
      for (std::size_t c = 0; c < free.size(); ++c) {
        double sum = 0;
        for (std::size_t d = 0; d < free.size(); ++d)
          for (std::size_t e = 0; e < free.size(); ++e)
            sum += outcome.covariance[free[a]][free[d]] * g[d][e] * outcome.covariance[free[e]][free[c]];
        corrected[free[a]][free[c]] = sum;
      }
    outcome.covariance = corrected;
  }
  outcome.errors.assign(n_par, 0.0);
  for (std::size_t i = 0; i < n_par; ++i) outcome.errors[i] = std::sqrt(std::max(outcome.covariance[i][i], 0.0));
  if (outcome.weighted && outcome.minos_ok && minos_index >= 0) {
    // Scale the profile interval by the sandwich / HESSE ratio of the parameter.
    const double hesse = hesse_error_before_sandwich;
    if (hesse > 0 && outcome.errors[minos_index] > 0) {
      const double ratio = outcome.errors[minos_index] / hesse;
      outcome.minos_low *= ratio;
      outcome.minos_high *= ratio;
    }
  }
  return outcome;
}
}  // namespace h4l
