// Small kinematics helpers shared by the v3 programs.
#pragma once

#include <Math/Vector4D.h>
#include <TVector2.h>

#include <cmath>

namespace h4l {

using P4 = ROOT::Math::PtEtaPhiMVector;

inline double delta_phi(double a, double b) { return TVector2::Phi_mpi_pi(a - b); }

inline double delta_r(double eta1, double phi1, double eta2, double phi2) {
  return std::hypot(eta1 - eta2, delta_phi(phi1, phi2));
}

constexpr double kMuonMass = 0.1056584;
constexpr double kElectronMass = 0.000511;
constexpr double kZMass = 91.1876;
}  // namespace h4l
