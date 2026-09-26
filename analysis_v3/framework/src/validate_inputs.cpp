// validate_inputs: integrity and bookkeeping of the v3 input files.
//
// For every file of one task it checks that the file opens cleanly, that the
// expected keys and every required branch exist, and that the first and the
// last Events entries read.  For MC it also sums the Runs bookkeeping, reads
// the provenance record and sums genWeight over all entries (which reads
// every basket of that branch).  For the pseudo-data, which may never be
// skipped, it reads every entry of every required branch, and it reads only
// the schema and the luminosity from the provenance record, as the blinding
// rules allow.  A file whose check fails is re-checked from a fresh open a
// bounded number of times before it is reported.
//
// Usage: validate_inputs --task TASK.json --out-json OUTPUT.json
// Exit code 0: the report was written; MC files that stay unreadable or
// invalid are recorded in it.  Exit code 3: a pseudo-data shard is not valid,
// or no file of the task could be opened (a node problem is more likely than
// a whole task of bad files); no report is written, so the task is retried.

#include "h4l/hash.h"
#include "h4l/io.h"
#include "h4l/root_io.h"

#include <TBranch.h>
#include <TLeaf.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace {
using h4l::json;

std::string hex64(std::uint64_t value) {
  char buffer[17];
  std::snprintf(buffer, sizeof(buffer), "%016llx", static_cast<unsigned long long>(value));
  return buffer;
}

// Name and leaf type of every Events branch, in name order.
std::vector<std::string> branch_signature(TTree& tree) {
  std::vector<std::string> signature;
  TObjArray* branches = tree.GetListOfBranches();
  for (int index = 0; index < branches->GetEntriesFast(); ++index) {
    auto* branch = static_cast<TBranch*>(branches->At(index));
    auto* leaf = static_cast<TLeaf*>(branch->GetListOfLeaves()->At(0));
    signature.push_back(std::string(branch->GetName()) + ":" + (leaf ? leaf->GetTypeName() : "?"));
  }
  std::sort(signature.begin(), signature.end());
  return signature;
}

std::string join(const std::vector<std::string>& items) {
  std::ostringstream text;
  for (std::size_t index = 0; index < items.size(); ++index) text << (index ? "; " : "") << items[index];
  return text.str();
}

// One check of one file.  Throws only on programming errors; every problem
// with the file itself is reported in the record.
json check_file(const json& input, const json& config) {
  const std::string path = input.at("path").get<std::string>();
  const bool is_mc = config.at("kind").get<std::string>() == "mc";
  json record = {{"file_key", input.at("file_key")}, {"path", path}, {"relative", input.at("relative")}};
  std::vector<std::string> problems;
  std::vector<std::string> warnings;

  std::string open_error;
  auto file = h4l::open_input(path, config.value("open_attempts", 3), open_error);
  if (!file) {
    record["status"] = "unreadable";
    record["reason"] = open_error;
    return record;
  }
  record["size_bytes"] = file->GetSize();

  json keys = json::array();
  for (TObject* key : *file->GetListOfKeys()) keys.push_back(key->GetName());
  record["keys"] = keys;

  auto* events = dynamic_cast<TTree*>(file->Get("Events"));
  if (!events) {
    record["status"] = "invalid";
    record["reason"] = "no Events tree";
    return record;
  }
  const Long64_t entries = events->GetEntries();
  record["events_entries"] = entries;
  if (entries <= 0) problems.push_back("Events has no entries");

  const auto signature = branch_signature(*events);
  std::set<std::string> names;
  std::string joined;
  for (const auto& item : signature) {
    names.insert(item.substr(0, item.find(':')));
    joined += item + "\n";
  }
  record["branch_count"] = signature.size();
  record["branch_signature_fnv1a64"] = hex64(h4l::fnv1a64(joined));

  json missing = json::array();
  for (const auto& name : config.at("required_events"))
    if (!names.count(name.get<std::string>())) missing.push_back(name);
  record["missing_required"] = missing;
  if (!missing.empty()) problems.push_back(std::to_string(missing.size()) + " required branches missing");
  json recorded = json::object();
  for (const auto& name : config.at("recorded_events")) recorded[name.get<std::string>()] = names.count(name.get<std::string>()) > 0;
  record["recorded_present"] = recorded;

  if (entries > 0) {
    const bool first = events->GetEntry(0) > 0;
    const bool last = events->GetEntry(entries - 1) > 0;
    record["first_last_entry_ok"] = first && last;
    if (!(first && last)) problems.push_back("first or last Events entry unreadable");
  }

  // The pseudo-data can never be skipped: read every entry of every
  // required branch, so that a corrupt basket anywhere is found now.
  if (config.value("full_read", false) && entries > 0 && missing.empty()) {
    events->SetBranchStatus("*", false);
    for (const auto& name : config.at("required_events")) events->SetBranchStatus(name.get<std::string>().c_str(), true);
    Long64_t bad_entry = -1;
    for (Long64_t index = 0; index < entries; ++index) {
      if (events->GetEntry(index) <= 0) {
        bad_entry = index;
        break;
      }
    }
    events->SetBranchStatus("*", true);
    record["full_read_ok"] = bad_entry < 0;
    if (bad_entry >= 0) problems.push_back("unreadable entry " + std::to_string(bad_entry) + " in the required branches");
  }

  const std::string provenance_name = config.at("provenance_object").get<std::string>();
  auto* provenance = dynamic_cast<TNamed*>(file->Get(provenance_name.c_str()));
  json provenance_json;
  if (!provenance) {
    problems.push_back("no " + provenance_name);
  } else {
    try {
      provenance_json = json::parse(provenance->GetTitle());
    } catch (const std::exception&) {
      problems.push_back(provenance_name + " is not JSON");
    }
    if (!provenance_json.is_null() && !provenance_json.is_object()) problems.push_back(provenance_name + " is not a JSON object");
  }
  const std::string expected_schema = config.at("expected_schema").get<std::string>();

  if (is_mc) {
    json mc = json::object();
    auto* runs = dynamic_cast<TTree*>(file->Get("Runs"));
    if (!runs || runs->GetEntries() < 1) {
      problems.push_back("Runs tree missing or empty");
    } else {
      Long64_t count = 0;
      Double_t sumw = 0, sumw2 = 0;
      const bool addressed = runs->SetBranchAddress("genEventCount", &count) >= 0 &&
                             runs->SetBranchAddress("genEventSumw", &sumw) >= 0 &&
                             runs->SetBranchAddress("genEventSumw2", &sumw2) >= 0;
      if (!addressed) {
        problems.push_back("Runs lacks genEventCount/genEventSumw/genEventSumw2 of the expected types");
      } else {
        long long total_count = 0;
        double total_sumw = 0, total_sumw2 = 0;
        for (Long64_t index = 0; index < runs->GetEntries(); ++index) {
          if (runs->GetEntry(index) <= 0) {
            problems.push_back("unreadable Runs entry");
            break;
          }
          total_count += count;
          total_sumw += sumw;
          total_sumw2 += sumw2;
        }
        mc["runs_entries"] = runs->GetEntries();
        mc["genEventCount"] = total_count;
        mc["genEventSumw"] = total_sumw;
        mc["genEventSumw2"] = total_sumw2;
        if (!std::isfinite(total_sumw) || total_sumw <= 0) problems.push_back("non-positive or non-finite genEventSumw");
      }
      runs->ResetBranchAddresses();
      // Theory-variation normalizations before the preselection: NanoAOD
      // stores LHE{Scale,Pdf}Sumw[i] = sum(genWeight * w_i) / genEventSumw per
      // Runs entry; record the absolute sums so that files can be added.
      for (const std::string name : {"LHEScaleSumw", "LHEPdfSumw"}) {
        const std::string counter = "n" + name;
        if (!runs->GetBranch(name.c_str()) || !runs->GetBranch(counter.c_str())) continue;
        UInt_t length = 0;
        std::vector<Double_t> values(512, 0.0);
        Double_t sumw_entry = 0;
        const bool ok = runs->SetBranchAddress(counter.c_str(), &length) >= 0 &&
                        runs->SetBranchAddress(name.c_str(), values.data()) >= 0 &&
                        runs->SetBranchAddress("genEventSumw", &sumw_entry) >= 0;
        if (!ok) {
          problems.push_back(name + " has an unexpected type");
        } else {
          std::vector<double> absolute;
          bool consistent = true;
          for (Long64_t index = 0; index < runs->GetEntries(); ++index) {
            if (runs->GetEntry(index) <= 0 || length > values.size()) {
              consistent = false;
              break;
            }
            if (absolute.empty()) absolute.assign(length, 0.0);
            if (absolute.size() != length) {
              consistent = false;
              break;
            }
            for (UInt_t item = 0; item < length; ++item) absolute[item] += values[item] * sumw_entry;
          }
          if (!consistent) problems.push_back(name + " is unreadable or changes length between Runs entries");
          else mc[name + "_absolute"] = absolute;
        }
        runs->ResetBranchAddresses();
      }
    }
    if (provenance_json.is_object()) {
      mc["provenance_schema"] = provenance_json.value("schema", "");
      mc["provenance_sample"] = provenance_json.value("sample", "");
      mc["entries_before_selection"] = provenance_json.value("entries_before_selection", -1LL);
      mc["entries_after_selection"] = provenance_json.value("entries_after_selection", -1LL);
      if (mc["provenance_schema"] != expected_schema) problems.push_back("unexpected provenance schema");
      if (mc["provenance_sample"] != config.at("sample")) problems.push_back("provenance sample differs from the catalogue");
      if (mc["entries_before_selection"].get<long long>() <= 0) problems.push_back("provenance lacks entries_before_selection");
      if (mc["entries_after_selection"].get<long long>() != entries)
        problems.push_back("provenance entries_after_selection differs from the Events entries");
      const json finalization = provenance_json.value("finalization", json::object());
      if (finalization.is_object() && finalization.contains("entries_after_selection") &&
          finalization["entries_after_selection"] != entries)
        problems.push_back("provenance finalization.entries_after_selection differs from the Events entries");
      // The Runs sums must describe exactly the events this file was made from.
      if (mc.contains("genEventCount") && mc["genEventCount"].get<long long>() != mc["entries_before_selection"].get<long long>())
        problems.push_back("Runs genEventCount differs from entries_before_selection");
    }
    // Lengths of the per-event weight vectors in the first entry.
    if (entries > 0 && events->GetEntry(0) > 0) {
      for (const std::string counter : {"nLHEScaleWeight", "nLHEPdfWeight", "nPSWeight"}) {
        if (auto* leaf = events->GetLeaf(counter.c_str())) mc["first_entry_" + counter] = leaf->GetValue();
      }
    }
    // Sum genWeight over every entry; this reads every basket of the branch.
    if (names.count("genWeight") && entries > 0) {
      events->SetBranchStatus("*", false);
      events->SetBranchStatus("genWeight", true);
      Float_t weight = 0;
      if (events->SetBranchAddress("genWeight", &weight) < 0) {
        problems.push_back("genWeight has an unexpected type");
      } else {
        double sum = 0, sum2 = 0;
        long long negative = 0, non_finite = 0;
        bool readable = true;
        for (Long64_t index = 0; index < entries; ++index) {
          if (events->GetEntry(index) <= 0) {
            readable = false;
            break;
          }
          if (!std::isfinite(weight)) {
            ++non_finite;
            continue;
          }
          sum += weight;
          sum2 += static_cast<double>(weight) * weight;
          negative += weight < 0;
        }
        if (!readable) problems.push_back("unreadable genWeight basket");
        if (non_finite) problems.push_back(std::to_string(non_finite) + " non-finite genWeight values");
        if (readable && !non_finite) {
          mc["sum_genWeight"] = sum;
          mc["sum_genWeight2"] = sum2;
          mc["n_negative_genWeight"] = negative;
        }
      }
      events->ResetBranchAddresses();
      events->SetBranchStatus("*", true);
    }
    record["mc"] = mc;
  } else {
    json data = json::object();
    auto* runs = dynamic_cast<TTree*>(file->Get("Runs"));
    auto* lumis = dynamic_cast<TTree*>(file->Get("LuminosityBlocks"));
    data["runs_entries"] = runs ? runs->GetEntries() : -1;
    data["lumiblocks_entries"] = lumis ? lumis->GetEntries() : -1;
    if (provenance_json.is_object()) {
      // Only the schema and the luminosity are read (blinding rules).
      data["provenance_schema"] = provenance_json.value("schema", "");
      data["lumi_fb"] = provenance_json.contains("lumi_fb") ? provenance_json["lumi_fb"] : json();
      if (data["provenance_schema"] != expected_schema) problems.push_back("unexpected provenance schema");
      const json expected_lumi = config.value("expected_lumi_fb", json());
      if (!data["lumi_fb"].is_number() || !expected_lumi.is_number() ||
          std::abs(data["lumi_fb"].get<double>() - expected_lumi.get<double>()) > 1e-9)
        problems.push_back("provenance lumi_fb differs from the catalogue");
    }
    record["data"] = data;
  }

  record["status"] = problems.empty() ? "ok" : "invalid";
  if (!problems.empty()) record["reason"] = join(problems);
  record["warnings"] = warnings;
  return record;
}

// Re-check a failing file from a fresh open, so that a transient EOS read
// error is not recorded as a bad file.
json validate_file(const json& input, const json& config) {
  const auto started = std::chrono::steady_clock::now();
  const int attempts = config.value("check_attempts", 3);
  json record;
  json history = json::array();
  for (int attempt = 1; attempt <= attempts; ++attempt) {
    try {
      record = check_file(input, config);
    } catch (const std::exception& error) {
      record = {{"file_key", input.at("file_key")}, {"path", input.at("path")}, {"relative", input.at("relative")},
                {"status", "invalid"}, {"reason", std::string("exception: ") + error.what()}};
    }
    if (record.at("status") == "ok") break;
    history.push_back(record.value("reason", ""));
    if (attempt < attempts) std::this_thread::sleep_for(std::chrono::seconds(20 * attempt));
  }
  record["check_attempts"] = static_cast<int>(history.size()) + (record.at("status") == "ok" ? 1 : 0);
  if (!history.empty()) record["failed_attempt_reasons"] = history;
  record["seconds"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
  return record;
}
}  // namespace

int main(int argc, char** argv) {
  std::string task_path, output_path;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--task" && index + 1 < argc) task_path = argv[++index];
    else if (argument == "--out-json" && index + 1 < argc) output_path = argv[++index];
    else {
      std::cerr << "usage: validate_inputs --task TASK.json --out-json OUTPUT.json\n";
      return 2;
    }
  }
  if (task_path.empty() || output_path.empty()) {
    std::cerr << "usage: validate_inputs --task TASK.json --out-json OUTPUT.json\n";
    return 2;
  }
  gROOT->SetBatch(true);
  try {
    const json task = h4l::read_json(task_path);
    const json& config = task.at("config");
    const bool is_data = config.at("kind").get<std::string>() == "data";
    json report = {{"schema", "h4l_v3_validation_report/2"},
                   {"task_id", task.at("task_id")},
                   {"sample", config.at("sample")},
                   {"kind", config.at("kind")},
                   {"config_fnv1a64", hex64(h4l::fnv1a64(config.dump()))},
                   {"frozen_program", task.value("frozen_program", json())},
                   {"started_utc", h4l::utc_now()},
                   {"files", json::array()}};
    int opened = 0;
    std::vector<std::string> bad_data;
    for (const auto& input : task.at("inputs")) {
      json record = validate_file(input, config);
      const std::string status = record.at("status").get<std::string>();
      std::cout << "[validate] " << status << " " << input.at("relative").get<std::string>() << "\n";
      opened += status != "unreadable";
      if (is_data && status != "ok") bad_data.push_back(input.at("relative").get<std::string>() + ": " + record.value("reason", ""));
      report["files"].push_back(record);
    }
    if (!bad_data.empty()) {
      std::cerr << "ERROR: pseudo-data shards are not valid (never skippable): " << join(bad_data) << "\n";
      return 3;
    }
    if (opened == 0) {
      std::cerr << "ERROR: no file of the task could be opened; not publishing (node problem?)\n";
      return 3;
    }
    report["finished_utc"] = h4l::utc_now();
    h4l::write_json_file(output_path, report);
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
