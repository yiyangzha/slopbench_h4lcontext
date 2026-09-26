// skim_histograms: validation histograms of the stage-2 skims of one sample.
//
// Usage: skim_histograms OUTPUT.root SKIM.root [SKIM.root ...]
// Pairs rows whose event fired the analysis trigger OR (mask given by the
// environment variable H4L_TRIGGER_OR_MASK, decimal) fill, per flavour
// (mm = muon pairs, ee = electron pairs): the pair mass (and with FSR), both
// legs' pT and eta, and for pairs in 80-100 GeV the leg isolation, SIP,
// impact parameters, electron BDT scores and muon/electron ID flags.  The
// first pair of every event fills the event-level histograms (vertices,
// rho, MET).  Every histogram is weighted by genWeight only (1 for data);
// the luminosity normalization is applied by the plotting script.

#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: skim_histograms OUTPUT.root SKIM.root [SKIM.root ...]\n";
    return 2;
  }
  const char* mask_text = std::getenv("H4L_TRIGGER_OR_MASK");
  if (!mask_text) {
    std::cerr << "ERROR: set H4L_TRIGGER_OR_MASK\n";
    return 2;
  }
  const unsigned long mask = std::stoul(mask_text);
  gROOT->SetBatch(true);
  TH1::SetDefaultSumw2(true);
  std::map<std::string, std::unique_ptr<TH1D>> h;
  auto book = [&](const std::string& name, int bins, double low, double high) {
    h[name] = std::make_unique<TH1D>(name.c_str(), name.c_str(), bins, low, high);
    h[name]->SetDirectory(nullptr);
  };
  for (const std::string flavour : {"mm", "ee"}) {
    book("mass_" + flavour, 100, 40, 140);
    book("mass_fsr_" + flavour, 100, 40, 140);
    book("mass_zoom_" + flavour, 80, 81, 101);
    book("pair_pt_" + flavour, 50, 0, 200);
    for (const std::string leg : {"l1", "l2"}) {
      book(leg + "_pt_" + flavour, 60, 0, 150);
      book(leg + "_eta_" + flavour, 50, -2.5, 2.5);
    }
    book("iso03_" + flavour, 50, 0, 1.0);
    book("sip_" + flavour, 50, 0, 10);
    book("dxy_" + flavour, 50, -0.1, 0.1);
    book("dz_" + flavour, 50, -0.2, 0.2);
    book("pt_err_rel_" + flavour, 50, 0, 0.1);
  }
  book("mva_noiso_ee", 50, -1, 1);
  book("flags_mm", 16, -0.5, 15.5);
  book("flags_ee", 16, -0.5, 15.5);
  book("probes_mm", 1, 0, 1);
  book("probes_ee", 1, 0, 1);
  book("npv", 70, 0, 70);
  book("npv_good", 70, 0, 70);
  book("rho", 60, 0, 60);
  book("met", 50, 0, 200);
  book("events", 1, 0, 1);

  for (int index = 2; index < argc; ++index) {
    std::unique_ptr<TFile> file(TFile::Open(argv[index], "READ"));
    if (!file || file->IsZombie()) {
      std::cerr << "ERROR: cannot open " << argv[index] << "\n";
      return 1;
    }
    auto* tree = dynamic_cast<TTree*>(file->Get("Pairs"));
    if (!tree) {
      std::cerr << "ERROR: no Pairs in " << argv[index] << "\n";
      return 1;
    }
    ULong64_t key = 0, last_key = ~0ULL;
    Long64_t entry = 0, last_entry = -1;
    Float_t weight = 1, rho = 0, met = 0, mass = 0, mass_fsr = 0, pair_pt = 0;
    Int_t npv = 0, npv_good = 0, flavour = 0;
    UInt_t trig = 0;
    Float_t pt[2], eta[2], iso[2], sip[2], dxy[2], dz[2], pterr[2], mva[2];
    UShort_t flags[2];
    tree->SetBranchStatus("*", false);
    auto address = [&](const char* name, void* where) {
      tree->SetBranchStatus(name, true);
      if (tree->SetBranchAddress(name, where) < 0) throw std::runtime_error(std::string("cannot address ") + name);
    };
    try {
      address("file_key", &key);
      address("entry", &entry);
      address("weight", &weight);
      address("rho", &rho);
      address("met", &met);
      address("npv", &npv);
      address("npv_good", &npv_good);
      address("trig", &trig);
      address("flavour", &flavour);
      address("mass", &mass);
      address("mass_fsr", &mass_fsr);
      address("pair_pt", &pair_pt);
      const char* legs[2] = {"l1_", "l2_"};
      for (int l = 0; l < 2; ++l) {
        const std::string p = legs[l];
        address((p + "pt").c_str(), &pt[l]);
        address((p + "eta").c_str(), &eta[l]);
        address((p + "iso03").c_str(), &iso[l]);
        address((p + "sip").c_str(), &sip[l]);
        address((p + "dxy").c_str(), &dxy[l]);
        address((p + "dz").c_str(), &dz[l]);
        address((p + "pt_err").c_str(), &pterr[l]);
        address((p + "mva_noiso").c_str(), &mva[l]);
        address((p + "flags").c_str(), &flags[l]);
      }
    } catch (const std::exception& error) {
      std::cerr << "ERROR: " << error.what() << " in " << argv[index] << "\n";
      return 1;
    }
    const Long64_t rows = tree->GetEntries();
    for (Long64_t row = 0; row < rows; ++row) {
      if (tree->GetEntry(row) <= 0) {
        std::cerr << "ERROR: unreadable Pairs row " << row << " in " << argv[index] << "\n";
        return 1;
      }
      if ((trig & mask) == 0) continue;
      const std::string f = flavour == 13 ? "mm" : "ee";
      if (key != last_key || entry != last_entry) {
        h["npv"]->Fill(npv, weight);
        h["npv_good"]->Fill(npv_good, weight);
        h["rho"]->Fill(rho, weight);
        h["met"]->Fill(met, weight);
        h["events"]->Fill(0.5, weight);
        last_key = key;
        last_entry = entry;
      }
      h["mass_" + f]->Fill(mass, weight);
      h["mass_fsr_" + f]->Fill(mass_fsr, weight);
      h["mass_zoom_" + f]->Fill(mass, weight);
      h["pair_pt_" + f]->Fill(pair_pt, weight);
      h["l1_pt_" + f]->Fill(pt[0], weight);
      h["l2_pt_" + f]->Fill(pt[1], weight);
      h["l1_eta_" + f]->Fill(eta[0], weight);
      h["l2_eta_" + f]->Fill(eta[1], weight);
      if (mass > 80 && mass < 100) {
        for (int l = 0; l < 2; ++l) {
          h["iso03_" + f]->Fill(iso[l], weight);
          h["sip_" + f]->Fill(sip[l], weight);
          h["dxy_" + f]->Fill(dxy[l], weight);
          h["dz_" + f]->Fill(dz[l], weight);
          h["pt_err_rel_" + f]->Fill(pterr[l] / pt[l], weight);
          if (flavour == 11) h["mva_noiso_ee"]->Fill(mva[l], weight);
          h["probes_" + f]->Fill(0.5, weight);
          for (int bit = 0; bit < 16; ++bit)
            if (flags[l] >> bit & 1) h["flags_" + f]->Fill(bit, weight);
        }
      }
    }
  }
  std::unique_ptr<TFile> output(TFile::Open(argv[1], "RECREATE"));
  if (!output || output->IsZombie()) {
    std::cerr << "ERROR: cannot create " << argv[1] << "\n";
    return 1;
  }
  for (auto& [name, histogram] : h) histogram->Write();
  output->Close();
  std::cout << "[skim_histograms] " << argv[1] << " from " << argc - 2 << " files\n";
  return 0;
}
