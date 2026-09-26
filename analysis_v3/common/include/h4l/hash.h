// Stable 64-bit hashes and the seeds of the MC smearing.
//
// The pseudo-data have no event identifiers, so a seed is built from the
// manifest file key of the ORIGINAL input file, the ORIGINAL entry number in
// that file, the object collection and the object index.  Skims must carry
// these originals (file key, entry, object index) so that every stage smears
// the same lepton identically.  One standard-normal deviate is drawn per
// lepton; a smearing variation scales that same deviate by a different width,
// so the up/down templates are correlated with the nominal.  `stream`
// separates statistically independent uses (for example toys); it is 0 for
// the analysis itself.
#pragma once

#include <cmath>
#include <cstdint>
#include <string>
#include <string_view>

namespace h4l {

inline std::uint64_t fnv1a64(std::string_view text) {
  std::uint64_t hash = 14695981039346656037ULL;
  for (unsigned char character : text) {
    hash ^= character;
    hash *= 1099511628211ULL;
  }
  return hash;
}

inline std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9E3779B97F4A7C15ULL;
  value = (value ^ (value >> 30U)) * 0xBF58476D1CE4E5B9ULL;
  value = (value ^ (value >> 27U)) * 0x94D049BB133111EBULL;
  return value ^ (value >> 31U);
}

// Manifest file keys are 16-digit hexadecimal FNV-1a values.
inline std::uint64_t parse_file_key(const std::string& hex) { return std::stoull(hex, nullptr, 16); }

enum class Collection : std::uint64_t { Muon = 1, Electron = 2, Photon = 3, Jet = 4 };

inline std::uint64_t object_seed(std::uint64_t file_key, std::uint64_t entry, Collection collection,
                                 std::uint64_t object_index, std::uint64_t stream = 0) {
  std::uint64_t seed = splitmix64(file_key);
  seed = splitmix64(seed ^ entry);
  seed = splitmix64(seed ^ (static_cast<std::uint64_t>(collection) << 56U) ^ object_index);
  return splitmix64(seed ^ (stream * 0xD1B54A32D192ED03ULL));
}

// The per-object standard-normal deviate used by every smearing variation.
// Box-Muller on two splitmix64 uniforms, so that the value is the same with
// every compiler and standard library (std::normal_distribution is
// implementation-defined).
inline double object_normal(std::uint64_t file_key, std::uint64_t entry, Collection collection,
                            std::uint64_t object_index, std::uint64_t stream = 0) {
  const std::uint64_t seed = object_seed(file_key, entry, collection, object_index, stream);
  const std::uint64_t a = splitmix64(seed), b = splitmix64(seed ^ 0xA0761D6478BD642FULL);
  const double u1 = (static_cast<double>(a >> 11U) + 0.5) * 0x1.0p-53;  // in (0, 1)
  const double u2 = static_cast<double>(b >> 11U) * 0x1.0p-53;          // in [0, 1)
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
}

// The per-object uniform in (0, 1) of a stream (for example the efficiency thinning of the closure tests).
inline double object_uniform(std::uint64_t file_key, std::uint64_t entry, Collection collection,
                             std::uint64_t object_index, std::uint64_t stream) {
  const std::uint64_t seed = object_seed(file_key, entry, collection, object_index, stream);
  return (static_cast<double>(splitmix64(seed ^ 0x6A09E667F3BCC909ULL) >> 11U) + 0.5) * 0x1.0p-53;
}
}  // namespace h4l
