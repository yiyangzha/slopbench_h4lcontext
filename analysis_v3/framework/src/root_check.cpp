// root_check: verify a published or about-to-be-published ROOT output.
//
// Usage: root_check FILE [TREE=ENTRIES ...] [OBJECT ...]
// The file must open cleanly (no recovery), every TREE must have exactly
// ENTRIES entries with readable first and last entries, and every OBJECT
// must exist.  Exit code 0 on success, 1 otherwise.

#include "h4l/root_io.h"

#include <TROOT.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "usage: root_check FILE [TREE=ENTRIES ...] [OBJECT ...]\n";
    return 2;
  }
  gROOT->SetBatch(true);
  std::map<std::string, long long> trees;
  std::vector<std::string> objects;
  try {
    for (int index = 2; index < argc; ++index) {
      const std::string item = argv[index];
      const auto equal = item.find('=');
      if (equal == std::string::npos) objects.push_back(item);
      else trees[item.substr(0, equal)] = std::stoll(item.substr(equal + 1));
    }
    h4l::validate_root_output(argv[1], trees, objects);
  } catch (const std::exception& error) {
    std::cerr << "root_check: " << error.what() << "\n";
    return 1;
  }
  std::cout << "root_check: ok " << argv[1] << "\n";
  return 0;
}
