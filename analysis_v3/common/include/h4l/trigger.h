// Matching of offline leptons to HLT trigger objects (NanoAOD TrigObj).
//
// * UL16 muons often carry near-duplicate HLT objects that differ in their
//   filter bits, so muon objects closer than `merge_dr` are merged into one
//   object (OR of the bits, highest pT, position of the highest-pT member).
// * A match ORs the bits of every object inside the cone and keeps the
//   highest pT and the smallest distance.
// * HLT electron objects sit at the supercluster position: eta = eta_SC and
//   phi = phi_track - q * 0.3 * B * r / (2 pT), with r the transverse radius
//   at which the electron reaches the calorimeter (barrel radius, or endcap
//   z / |sinh eta_SC|).
#pragma once

#include "h4l/kinematics.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace h4l {

struct TrigObject {
  float pt = 0, eta = 0, phi = 0;
  int bits = 0, members = 0;
};

struct TrigMatch {
  int bits = 0, objects = 0;
  float pt = -1, dr = -1;
};

struct Bending {
  double field_tesla = 3.8, barrel_radius_m = 1.29, endcap_z_m = 3.14;
};

// Objects of one id; muons are merged within merge_dr (use 0 for no merging).
template <typename F, typename I>
std::vector<TrigObject> trigger_objects(const F& pt, const F& eta, const F& phi, const I& id, const I& bits, int wanted_id,
                                        double merge_dr) {
  std::vector<std::size_t> order;
  for (std::size_t k = 0; k < pt.GetSize(); ++k)
    if (id[k] == wanted_id) order.push_back(k);
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return pt[a] > pt[b]; });
  std::vector<TrigObject> objects;
  for (std::size_t k : order) {
    bool merged = false;
    for (auto& object : objects) {
      if (merge_dr > 0 && delta_r(object.eta, object.phi, eta[k], phi[k]) < merge_dr) {
        object.bits |= bits[k];
        ++object.members;
        merged = true;
        break;
      }
    }
    if (!merged) objects.push_back({pt[k], eta[k], phi[k], bits[k], 1});
  }
  return objects;
}

inline TrigMatch match_trigger(const std::vector<TrigObject>& objects, double eta, double phi, double cone) {
  TrigMatch match;
  for (const auto& object : objects) {
    const double dr = delta_r(eta, phi, object.eta, object.phi);
    if (dr >= cone) continue;
    match.bits |= object.bits;
    match.pt = std::max(match.pt, object.pt);
    match.dr = match.objects == 0 ? static_cast<float>(dr) : std::min(match.dr, static_cast<float>(dr));
    ++match.objects;
  }
  return match;
}

// Predicted (eta, phi) of the HLT object of an electron.
inline std::pair<double, double> electron_hlt_position(double eta_sc, double phi, int charge, double pt, const Bending& b) {
  const double sinh_eta = std::fabs(std::sinh(eta_sc));
  const double radius = sinh_eta > 0 ? std::min(b.barrel_radius_m, b.endcap_z_m / sinh_eta) : b.barrel_radius_m;
  const double bend = 0.3 * b.field_tesla * radius / (2.0 * pt);
  return {eta_sc, TVector2::Phi_mpi_pi(phi - charge * bend)};
}
}  // namespace h4l
