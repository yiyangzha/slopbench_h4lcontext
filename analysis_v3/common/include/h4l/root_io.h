// ROOT helpers of the v3 programs: robust opening of inputs and the
// close / reopen / keys / entries validation of every ROOT output.
#pragma once

#include <TFile.h>
#include <TKey.h>
#include <TTree.h>

#include <chrono>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace h4l {

// Open an input for reading with a bounded number of attempts.  A file that
// ROOT had to recover is treated as unreadable: its content is not trusted.
inline std::unique_ptr<TFile> open_input(const std::string& path, int attempts, std::string& error) {
  for (int attempt = 1; attempt <= attempts; ++attempt) {
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
    if (file && !file->IsZombie() && !file->TestBit(TFile::kRecovered)) return file;
    error = !file ? "TFile::Open returned null" : (file->IsZombie() ? "zombie file" : "file needed recovery");
    if (attempt < attempts) std::this_thread::sleep_for(std::chrono::seconds(5 * attempt));
  }
  return nullptr;
}

// Reopen a written ROOT file and verify that every required tree has exactly
// the expected entries and that its first and last entries read.
inline void validate_root_output(const std::string& path, const std::map<std::string, long long>& trees,
                                 const std::vector<std::string>& objects) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie() || file->TestBit(TFile::kRecovered))
    throw std::runtime_error("output does not reopen cleanly: " + path);
  for (const auto& [name, entries] : trees) {
    auto* tree = dynamic_cast<TTree*>(file->Get(name.c_str()));
    if (!tree) throw std::runtime_error("output lacks tree " + name + ": " + path);
    if (tree->GetEntries() != entries)
      throw std::runtime_error("output tree " + name + " has " + std::to_string(tree->GetEntries()) + " entries, expected " +
                               std::to_string(entries) + ": " + path);
    if (entries > 0 && (tree->GetEntry(0) <= 0 || tree->GetEntry(entries - 1) <= 0))
      throw std::runtime_error("output tree " + name + " has unreadable entries: " + path);
  }
  for (const auto& name : objects)
    if (!file->Get(name.c_str())) throw std::runtime_error("output lacks object " + name + ": " + path);
}

// The same checks, but every entry of every required tree is read with all
// branches: used on the local output before publication, so that a write
// error in the middle of a file cannot hide behind correct entry counts.
inline void validate_root_output_full(const std::string& path, const std::map<std::string, long long>& trees) {
  validate_root_output(path, trees, {});
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) throw std::runtime_error("cannot reopen the output " + path);
  for (const auto& [name, entries] : trees) {
    auto* tree = dynamic_cast<TTree*>(file->Get(name.c_str()));
    if (!tree) throw std::runtime_error("output tree " + name + " missing on reopening: " + path);
    for (long long entry = 0; entry < entries; ++entry)
      if (tree->GetEntry(entry) <= 0)
        throw std::runtime_error("output tree " + name + " has an unreadable entry " + std::to_string(entry) + ": " + path);
  }
}
}  // namespace h4l
