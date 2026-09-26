// Z-peak line shape of the lepton calibration: a Breit-Wigner of fixed mass
// and width convolved with a double-sided Crystal Ball, plus a background.
//
// The convolution (BW (x) DCB)(m) = int BW(m - x) DCB(x) dx is computed on a
// uniform grid of step h that divides the histogram bin width: the window
// [lo, hi) is sampled at the sub-bin midpoints m_a, the kernel at x_i = i h
// for |x_i| <= kernel_half, and the Breit-Wigner, which has no free
// parameter, is tabulated once.  Bin contents are midpoint sums over the
// sub-bins, normalized over the window, so no FFT is needed.
//
// DCB(x) = exp(-t^2/2) for -alphaL <= t <= alphaR, t = (x - delta) / sigma,
// with power-law tails A (B - t)^-n beyond them (as RooCrystalBall).
// BW(m) = 1 / ((m - mZ)^2 + Gamma^2 / 4) (as RooBreitWigner).
// Parameters: 0 delta, 1 sigma, 2 alphaL, 3 nL, 4 alphaR, 5 nR, 6 n_sig,
// 7 n_bkg, 8 and 9 background shape (exponential: 8 = tau; chebychev2:
// 8 = c1, 9 = c2 on x = 2 (m - lo) / (hi - lo) - 1; bernstein2_falling:
// coefficients 1, 8, 8 * 9 of (1-x)^2, 2 x (1-x), x^2 on x = (m - lo) /
// (hi - lo), with 8 and 9 in [0, 1], so the polynomial falls
// monotonically).  With the kernel "double_gaussian" (the alternative
// resolution of the tag-and-probe) the resolution is f G(delta, s1) +
// (1 - f) G(delta, s2): 0 delta, 1 s1, 2 s2, 3 f (4, 5 unused).  With
// "standalone_dcb" the signal is the double-sided Crystal Ball itself at
// m_Z + delta, without the Breit-Wigner (the line shape of probes binned in
// pT, whose natural Breit-Wigner tails are sculpted away).  With the kernel
// "template" the signal is a binned MC template (set_template: lower edge,
// bin width, contents) convolved with a Gaussian of mean delta (0) and width
// sigma (1): s(m) = sum_j T_j [Phi((m - delta - lo_j) / sigma) -
// Phi((m - delta - hi_j) / sigma)], the sum restricted to the template bins
// within 7 sigma (2-5 unused).  With "standalone_dcb_hump" the stand-alone
// DCB (normalized over the window) is mixed with a Gaussian of mean p[11],
// width p[12] and fraction p[13] (the parameter array must then hold 14
// values): the low-mass component of fail spectra (FSR photons inside the
// isolation cone, mismeasured probes).  Backgrounds with three shape parameters
// (8, 9, 10; the parameter array must then hold 11 values): "cmsshape"
// erfc((alpha - m) beta) exp(-gamma (m - m_Z)) (8 alpha, 9 beta, 10 gamma, as
// RooCMSShape with its peak at m_Z) and "bernstein3" with the coefficients
// 1, 8, 9, 10 of the cubic Bernstein basis on x = (m - lo) / (hi - lo)
// (non-negative coefficients keep it non-negative).
#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace h4l {

inline double dcb_shape(double t, double alpha_l, double n_l, double alpha_r, double n_r) {
  if (t < -alpha_l) {
    const double a = std::pow(n_l / alpha_l, n_l) * std::exp(-0.5 * alpha_l * alpha_l), b = n_l / alpha_l - alpha_l;
    return a * std::pow(b - t, -n_l);
  }
  if (t > alpha_r) {
    const double a = std::pow(n_r / alpha_r, n_r) * std::exp(-0.5 * alpha_r * alpha_r), b = n_r / alpha_r - alpha_r;
    return a * std::pow(b + t, -n_r);
  }
  return std::exp(-0.5 * t * t);
}

struct ZPeakModel {
  enum class Background { Exponential, Chebychev2, Bernstein2Falling, CMSShape, Bernstein3, None };
  enum class Kernel { DoubleCrystalBall, DoubleGaussian, StandaloneDCB, StandaloneDCBHump, Template };
  double lo = 0, hi = 0, bin_width = 0, mz = 91.1876, gz = 2.4952, step = 0.05, kernel_half = 30.0;
  int bins = 0, sub = 0, kernel = 0;
  Background background = Background::Exponential;
  Kernel kernel_type = Kernel::DoubleCrystalBall;
  std::vector<double> bw;  // BW at lo + (j + 0.5) h - kernel_half for j = 0 .. sub-points + 2 kernel
  double template_lo = 0, template_width = 0;
  std::vector<double> template_content, template2_content;

  void set_template(double low, double width, const std::vector<double>& content) {
    if (!(width > 0) || content.empty()) throw std::runtime_error("empty signal template");
    template_lo = low;
    template_width = width;
    template_content = content;
  }

  // A second template on the same grid: the Template kernel then mixes the two (each normalized in the
  // window) with the fraction p[13] of the second, clamped at zero (p[13] may be negative).
  void set_second_template(const std::vector<double>& content) {
    if (content.size() != template_content.size()) throw std::runtime_error("the second template needs the grid of the first");
    template2_content = content;
  }

  static Background parse_background(const std::string& name) {
    if (name == "exponential") return Background::Exponential;
    if (name == "chebychev2") return Background::Chebychev2;
    if (name == "bernstein2_falling") return Background::Bernstein2Falling;
    if (name == "cmsshape") return Background::CMSShape;
    if (name == "bernstein3") return Background::Bernstein3;
    if (name == "none") return Background::None;
    throw std::runtime_error("unknown background model " + name);
  }

  ZPeakModel(double low, double width, int n_bins, double z_mass, double z_width, const std::string& bkg,
             double grid_step = 0.05, double half = 30.0)
      : lo(low), hi(low + width * n_bins), bin_width(width), mz(z_mass), gz(z_width), step(grid_step), kernel_half(half),
        bins(n_bins) {
    sub = static_cast<int>(std::lround(bin_width / step));
    if (sub < 1 || std::fabs(sub * step - bin_width) > 1e-9) throw std::runtime_error("grid step must divide the bin width");
    kernel = static_cast<int>(std::lround(kernel_half / step));
    background = parse_background(bkg);
    const int points = bins * sub;
    bw.resize(points + 2 * kernel + 1);
    for (int j = 0; j < static_cast<int>(bw.size()); ++j) {
      const double m = lo + (j - kernel + 0.5) * step;
      bw[j] = 1.0 / ((m - mz) * (m - mz) + 0.25 * gz * gz);
    }
  }

  // Signal density at the sub-bin midpoints (unnormalized).
  void signal_grid(const double* p, std::vector<double>& out) const {
    const int points = bins * sub;
    if (kernel_type == Kernel::Template) {
      const double sigma = p[1], delta = p[0], root2 = std::sqrt(2.0);
      auto convolve = [&](const std::vector<double>& content, std::vector<double>& target) {
        target.assign(points, 0.0);
        const int n = static_cast<int>(content.size());
        for (int a = 0; a < points; ++a) {
          const double x = lo + (a + 0.5) * step - delta;
          const int first = std::max(0, static_cast<int>(std::floor((x - 7.0 * sigma - template_lo) / template_width)));
          const int last = std::min(n - 1, static_cast<int>(std::floor((x + 7.0 * sigma - template_lo) / template_width)));
          if (last < first) continue;
          double previous = 0.5 * std::erfc(-(x - (template_lo + first * template_width)) / (sigma * root2));
          double sum = 0;
          for (int j = first; j <= last; ++j) {
            const double upper_edge = template_lo + (j + 1) * template_width;
            const double current = 0.5 * std::erfc(-(x - upper_edge) / (sigma * root2));
            sum += content[j] * (previous - current);
            previous = current;
          }
          target[a] = sum;
        }
      };
      convolve(template_content, out);
      if (!template2_content.empty() && p[13] != 0.0) {
        std::vector<double> second;
        convolve(template2_content, second);
        double first_sum = 0, second_sum = 0;
        for (int a = 0; a < points; ++a) {
          first_sum += out[a];
          second_sum += second[a];
        }
        for (int a = 0; a < points; ++a)
          out[a] = std::max((1.0 - p[13]) * (first_sum > 0 ? out[a] / first_sum : 0.0) +
                                p[13] * (second_sum > 0 ? second[a] / second_sum : 0.0),
                            0.0);
      }
      return;
    }
    if (kernel_type == Kernel::StandaloneDCB) {
      out.assign(points, 0.0);
      for (int a = 0; a < points; ++a) {
        const double m = lo + (a + 0.5) * step;
        out[a] = dcb_shape((m - mz - p[0]) / p[1], p[2], p[3], p[4], p[5]);
      }
      return;
    }
    if (kernel_type == Kernel::StandaloneDCBHump) {
      out.assign(points, 0.0);
      std::vector<double> hump(points, 0.0);
      double dcb_sum = 0, hump_sum = 0;
      for (int a = 0; a < points; ++a) {
        const double m = lo + (a + 0.5) * step, z = (m - p[11]) / p[12];
        out[a] = dcb_shape((m - mz - p[0]) / p[1], p[2], p[3], p[4], p[5]);
        hump[a] = std::exp(-0.5 * z * z);
        dcb_sum += out[a];
        hump_sum += hump[a];
      }
      for (int a = 0; a < points; ++a)
        out[a] = (1.0 - p[13]) * (dcb_sum > 0 ? out[a] / dcb_sum : 0.0) + p[13] * (hump_sum > 0 ? hump[a] / hump_sum : 0.0);
      return;
    }
    std::vector<double> k(2 * kernel + 1);
    if (kernel_type == Kernel::DoubleGaussian) {
      for (int i = -kernel; i <= kernel; ++i) {
        const double x = i * step - p[0];
        k[i + kernel] = p[3] * std::exp(-0.5 * x * x / (p[1] * p[1])) / p[1] +
                        (1.0 - p[3]) * std::exp(-0.5 * x * x / (p[2] * p[2])) / p[2];
      }
    } else {
      for (int i = -kernel; i <= kernel; ++i) k[i + kernel] = dcb_shape((i * step - p[0]) / p[1], p[2], p[3], p[4], p[5]);
    }
    out.assign(points, 0.0);
    // m_a - x_i = lo + (a - i + 0.5) h = bw index (a - i + kernel).
    for (int a = 0; a < points; ++a) {
      double sum = 0;
      const double* bw_row = bw.data() + a + kernel;
      for (int i = -kernel; i <= kernel; ++i) sum += k[i + kernel] * bw_row[-i];
      out[a] = sum;
    }
  }

  double background_density(const double* p, double m) const {
    switch (background) {
      case Background::Exponential: return std::exp(p[8] * (m - lo));
      case Background::Chebychev2: {
        const double x = 2.0 * (m - lo) / (hi - lo) - 1.0;
        return 1.0 + p[8] * x + p[9] * (2.0 * x * x - 1.0);
      }
      case Background::Bernstein2Falling: {
        const double x = (m - lo) / (hi - lo);
        return (1 - x) * (1 - x) + 2.0 * p[8] * x * (1 - x) + p[8] * p[9] * x * x;
      }
      case Background::CMSShape: return std::erfc((p[8] - m) * p[9]) * std::exp(-p[10] * (m - mz));
      case Background::Bernstein3: {
        const double x = (m - lo) / (hi - lo), y = 1 - x;
        return y * y * y + 3.0 * p[8] * x * y * y + 3.0 * p[9] * x * x * y + p[10] * x * x * x;
      }
      case Background::None: return 0.0;
    }
    return 0.0;
  }

  // Expected bin contents: n_sig and n_bkg are the window yields.
  void expected(const double* p, std::vector<double>& nu, std::vector<double>* nu_bkg = nullptr) const {
    std::vector<double> s;
    signal_grid(p, s);
    double s_total = 0, b_total = 0;
    std::vector<double> s_bin(bins, 0.0), b_bin(bins, 0.0);
    for (int b = 0; b < bins; ++b) {
      for (int a = b * sub; a < (b + 1) * sub; ++a) {
        s_bin[b] += s[a];
        b_bin[b] += background_density(p, lo + (a + 0.5) * step);
      }
      s_total += s_bin[b];
      b_total += b_bin[b];
    }
    nu.assign(bins, 0.0);
    if (nu_bkg) nu_bkg->assign(bins, 0.0);
    for (int b = 0; b < bins; ++b) {
      const double signal = s_total > 0 ? p[6] * s_bin[b] / s_total : 0.0;
      const double bkg = background == Background::None || !(b_total > 0) ? 0.0 : p[7] * b_bin[b] / b_total;
      nu[b] = signal + bkg;
      if (nu_bkg) (*nu_bkg)[b] = bkg;
    }
  }
};
}  // namespace h4l
