// Shared file-system helpers of the v3 programs: EOS-safe directory creation,
// attempt-qualified temporary paths and atomic publication of text outputs.
#pragma once

#include <nlohmann/json.hpp>

#include <chrono>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <unistd.h>

namespace h4l {
namespace fs = std::filesystem;
using json = nlohmann::json;

inline std::string utc_now() {
  const std::time_t now = std::time(nullptr);
  char buffer[32];
  std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&now));
  return buffer;
}

// EOS can answer a concurrent mkdir with a transient EEXIST or ENOENT; the
// directory counts as created once it verifiably exists.
inline void ensure_directory(const fs::path& directory, int attempts = 5) {
  for (int attempt = 1; attempt <= attempts; ++attempt) {
    std::error_code error;
    fs::create_directories(directory, error);
    if (fs::is_directory(directory)) return;
    if (attempt < attempts) std::this_thread::sleep_for(std::chrono::seconds(3));
  }
  throw std::runtime_error("cannot create directory " + directory.string());
}

// A sibling path that no other attempt can produce: host, PID, time, random.
inline fs::path attempt_path(const fs::path& target, const std::string& tag) {
  std::random_device entropy;
  const auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::system_clock::now().time_since_epoch()).count();
  char host[64] = {0};
  ::gethostname(host, sizeof(host) - 1);
  return fs::path(target.string() + ".partial." + tag + "." + host + "." + std::to_string(::getpid()) + "." +
                  std::to_string(nanoseconds) + "." + std::to_string(entropy()));
}

inline json read_json(const fs::path& path) {
  std::ifstream stream(path);
  if (!stream) throw std::runtime_error("cannot read " + path.string());
  return json::parse(stream);
}

// Write, flush and re-read the text, then rename it onto `target`; an
// existing target is never overwritten.
inline void publish_text(const fs::path& target, const std::string& text) {
  ensure_directory(target.parent_path());
  const fs::path temporary = attempt_path(target, "text");
  {
    std::ofstream stream(temporary);
    stream << text;
    stream.flush();
    if (!stream) throw std::runtime_error("cannot write " + temporary.string());
  }
  std::ifstream check(temporary);
  const std::string reread((std::istreambuf_iterator<char>(check)), std::istreambuf_iterator<char>());
  if (reread != text) throw std::runtime_error("re-read mismatch for " + temporary.string());
  if (fs::exists(target)) throw std::runtime_error("refusing to overwrite " + target.string() + "; kept " + temporary.string());
  fs::rename(temporary, target);
}

// Invalid UTF-8 (for example in a ROOT key name) is replaced, not thrown.
inline std::string dump_json(const json& value) {
  return value.dump(1, ' ', false, json::error_handler_t::replace) + "\n";
}

inline void publish_json(const fs::path& target, const json& value) { publish_text(target, dump_json(value)); }

// Plain write used for the temporary outputs a worker wrapper asks for: the
// wrapper, not the program, performs the final publication.
inline void write_json_file(const fs::path& path, const json& value) {
  ensure_directory(path.parent_path().empty() ? fs::path(".") : path.parent_path());
  std::ofstream stream(path);
  stream << dump_json(value);
  stream.flush();
  if (!stream) throw std::runtime_error("cannot write " + path.string());
  stream.close();
  // The re-read must parse and reproduce the same serialization.
  if (dump_json(read_json(path)) != dump_json(value)) throw std::runtime_error("re-read mismatch for " + path.string());
}
}  // namespace h4l
