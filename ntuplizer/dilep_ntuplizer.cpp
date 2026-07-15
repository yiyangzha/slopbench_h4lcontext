// NanoAOD -> dilepton control-region ntuplizer: output schema and selection only.

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kMZ = 91.1876;
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kTwoPi = 2.0 * kPi;

constexpr std::array<const char *, 11> kHltPaths = {
    "IsoMu24",
    "IsoMu27",
    "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
    "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
    "Ele32_WPTight_Gsf",
    "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
    "Mu8_DiEle12_CaloIdL_TrackIdL_DZ",
    "DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
    "TripleMu_12_10_5",
    "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
    "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
};

std::string triggerBranchName(std::size_t bit) {
  return std::string("HLT_") + kHltPaths.at(bit);
}

struct Options {
  std::string input;
  std::string output;
};

struct Lepton {
  double pt = 0.0;
  double eta = 0.0;
  double phi = 0.0;
  double mass = 0.0;
  int charge = 0;
  int pdgId = 0;
  double dxy = 0.0;
  double dz = 0.0;
  double dxyErr = 0.0;
  double dzErr = 0.0;
  double sip3d = 0.0;
  double ip3d = 0.0;
  double pfRelIso03 = 0.0;
  double miniRelIso = 0.0;
  int muMedium = 0;
  int muTight = 0;
  int muGlobal = 0;
  int muPF = 0;
  int elCutBased = 0;
  int elMvaWP90 = 0;
  int elMvaWP80 = 0;
  double elMvaBdt = -99.0;
};

struct P4 {
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double e = 0.0;
};

struct Candidate {
  double mll = 0.0;
  double ptll = 0.0;
  double etall = 0.0;
  double phill = 0.0;
  double yll = 0.0;
  int finalState = 0;
  std::array<Lepton, 2> pair{};
};

struct OutputLepton {
  // Kinematics and particle identity; pdgId follows the NanoAOD convention.
  Float_t pt = 0.0F;
  Float_t eta = 0.0F;
  Float_t phi = 0.0F;
  Float_t mass = 0.0F;
  Int_t charge = 0;
  Int_t pdgId = 0;
  // Dilepton candidates use zId=1 for both selected leptons.
  Int_t zId = 0;
  // Transverse/longitudinal impact parameters and their uncertainties.
  Float_t dxy = 0.0F;
  Float_t dz = 0.0F;
  Float_t dxyErr = 0.0F;
  Float_t dzErr = 0.0F;
  Float_t sip3d = 0.0F;
  Float_t ip3d = 0.0F;
  Float_t pfRelIso03 = 0.0F;
  Float_t miniRelIso = 0.0F;
  // Input NanoAOD muon/electron identification quantities.
  Int_t muMedium = 0;
  Int_t muTight = 0;
  Int_t muGlobal = 0;
  Int_t muPF = 0;
  Int_t elCutBased = 0;
  Int_t elMvaWP90 = 0;
  Int_t elMvaWP80 = 0;
  Float_t elMvaBdt = 0.0F;
};

// `dilepTree` has 63 scalar branches: 10 event/PV, 7 candidate/category,
// and 2 * 23 lepton branches.  finalState is 0=mumu and 1=ee.
// trigBits stores the trigger decisions in the bit order of kHltPaths.
struct OutputRow {
  // Event identifier and primary-vertex observables.
  Int_t run = 0;
  Int_t lumi = 0;
  Long64_t event = 0;
  Int_t nPV = 0;
  Float_t pvX = 0.0F;
  Float_t pvY = 0.0F;
  Float_t pvZ = 0.0F;
  Float_t pvChi2 = 0.0F;
  Float_t pvNdof = 0.0F;
  Float_t pvScore = 0.0F;
  // Dilepton kinematics.
  Float_t mll = 0.0F;
  Float_t ptll = 0.0F;
  Float_t etall = 0.0F;
  Float_t phill = 0.0F;
  Float_t yll = 0.0F;
  Int_t finalState = 0;
  Int_t trigBits = 0;
  std::array<OutputLepton, 2> leptons{};
};

[[noreturn]] void fail(const std::string &message) {
  throw std::runtime_error(message);
}

void printUsage(std::ostream &stream) {
  stream << "Usage: dilep_ntuplizer --input FILE.root --output FILE.root\n";
}

std::string takeValue(int &index, int argc, char **argv,
                      const std::string &option) {
  if (++index >= argc)
    fail("Missing value for " + option);
  return argv[index];
}

Options parseOptions(int argc, char **argv) {
  Options options;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (option == "--help" || option == "-h") {
      printUsage(std::cout);
      std::exit(0);
    }
    if (option == "--input") {
      options.input = takeValue(index, argc, argv, option);
    } else if (option == "--output") {
      options.output = takeValue(index, argc, argv, option);
    } else {
      fail("Unknown option: " + option);
    }
  }

  if (options.input.empty() || options.output.empty()) {
    printUsage(std::cerr);
    fail("--input and --output are required");
  }
  return options;
}

P4 p4(const Lepton &lepton) {
  const double cphi = std::cos(lepton.phi);
  const double sphi = std::sin(lepton.phi);
  const double coshEta = std::cosh(lepton.eta);
  const double sinhEta = std::sinh(lepton.eta);
  return {lepton.pt * cphi, lepton.pt * sphi, lepton.pt * sinhEta,
          std::hypot(lepton.pt * coshEta, lepton.mass)};
}

P4 addP4(const P4 &left, const P4 &right) {
  return {left.px + right.px, left.py + right.py, left.pz + right.pz,
          left.e + right.e};
}

double invariantMass(const Lepton &first, const Lepton &second) {
  const P4 sum = addP4(p4(first), p4(second));
  return std::sqrt(std::max(sum.e * sum.e - sum.px * sum.px - sum.py * sum.py -
                                sum.pz * sum.pz,
                            0.0));
}

double deltaR2(const Lepton &first, const Lepton &second) {
  double dphi = std::fmod(first.phi - second.phi + kPi, kTwoPi);
  if (dphi < 0.0)
    dphi += kTwoPi;
  dphi -= kPi;
  const double deta = first.eta - second.eta;
  return deta * deta + dphi * dphi;
}

bool isBetterPair(double distance, double leadingPt, double subleadingPt,
                  std::size_t firstIndex, std::size_t secondIndex,
                  double bestDistance, double bestLeadingPt,
                  double bestSubleadingPt, std::size_t bestFirstIndex,
                  std::size_t bestSecondIndex) {
  // The primary ranking is closeness to mZ.  The remaining keys make equal
  // mass-distance cases reproducible without introducing a mass requirement.
  constexpr double kTieTolerance = 1e-12;
  if (distance < bestDistance - kTieTolerance)
    return true;
  if (distance > bestDistance + kTieTolerance)
    return false;
  if (leadingPt > bestLeadingPt + kTieTolerance)
    return true;
  if (leadingPt < bestLeadingPt - kTieTolerance)
    return false;
  if (subleadingPt > bestSubleadingPt + kTieTolerance)
    return true;
  if (subleadingPt < bestSubleadingPt - kTieTolerance)
    return false;
  if (firstIndex != bestFirstIndex)
    return firstIndex < bestFirstIndex;
  return secondIndex < bestSecondIndex;
}

// Retain one OSSF pair with deltaR > 0.02, leading/subleading pT >= 25/15
// GeV, and mass closest to the nominal Z mass.  There is no mll requirement.
Candidate findBestPair(const std::vector<Lepton> &leptons, bool &found) {
  found = false;
  Candidate best;
  double bestDistance = std::numeric_limits<double>::infinity();
  double bestLeadingPt = -std::numeric_limits<double>::infinity();
  double bestSubleadingPt = -std::numeric_limits<double>::infinity();
  std::size_t bestFirstIndex = std::numeric_limits<std::size_t>::max();
  std::size_t bestSecondIndex = std::numeric_limits<std::size_t>::max();

  for (std::size_t first = 0; first + 1 < leptons.size(); ++first) {
    for (std::size_t second = first + 1; second < leptons.size(); ++second) {
      const Lepton &left = leptons[first];
      const Lepton &right = leptons[second];
      if (std::abs(left.pdgId) != std::abs(right.pdgId) ||
          left.charge == right.charge)
        continue;
      if (deltaR2(left, right) < 4e-4)
        continue;

      const double leadingPt = std::max(left.pt, right.pt);
      const double subleadingPt = std::min(left.pt, right.pt);
      if (leadingPt < 25.0 || subleadingPt < 15.0)
        continue;

      const double mll = invariantMass(left, right);
      const double distance = std::abs(mll - kMZ);
      if (!isBetterPair(distance, leadingPt, subleadingPt, first, second,
                        bestDistance, bestLeadingPt, bestSubleadingPt,
                        bestFirstIndex, bestSecondIndex)) {
        continue;
      }

      std::array<Lepton, 2> pair = {left, right};
      if (pair[0].pt < pair[1].pt)
        std::swap(pair[0], pair[1]);
      const P4 total = addP4(p4(pair[0]), p4(pair[1]));
      const double ptll = std::hypot(total.px, total.py);
      best.mll = mll;
      best.ptll = ptll;
      best.etall = std::asinh(total.pz / (ptll + 1e-30));
      best.phill = std::atan2(total.py, total.px);
      best.yll =
          0.5 * std::log((total.e + total.pz) / (total.e - total.pz + 1e-30));
      best.finalState = std::abs(pair[0].pdgId) == 13 ? 0 : 1;
      best.pair = pair;
      bestDistance = distance;
      bestLeadingPt = leadingPt;
      bestSubleadingPt = subleadingPt;
      bestFirstIndex = first;
      bestSecondIndex = second;
      found = true;
    }
  }
  return best;
}

template <typename T>
std::unique_ptr<TTreeReaderValue<T>>
optionalValue(TTreeReader &reader, TTree &tree, const char *name) {
  if (tree.GetBranch(name) == nullptr)
    return nullptr;
  return std::make_unique<TTreeReaderValue<T>>(reader, name);
}

template <typename T>
std::unique_ptr<TTreeReaderArray<T>>
optionalArray(TTreeReader &reader, TTree &tree, const char *name) {
  if (tree.GetBranch(name) == nullptr)
    return nullptr;
  return std::make_unique<TTreeReaderArray<T>>(reader, name);
}

void requireBranch(TTree &tree, const char *name) {
  if (tree.GetBranch(name) == nullptr)
    fail("Input Events tree is missing required branch '" + std::string(name) +
         "'");
}

template <typename T>
std::unique_ptr<TTreeReaderValue<T>>
requiredValuePtr(TTreeReader &reader, TTree &tree, const char *name) {
  requireBranch(tree, name);
  return std::make_unique<TTreeReaderValue<T>>(reader, name);
}

template <typename T>
T valueOr(const std::unique_ptr<TTreeReaderValue<T>> &value, T fallback) {
  return value ? **value : fallback;
}

template <typename T>
T arrayOr(const std::unique_ptr<TTreeReaderArray<T>> &value, std::size_t index,
          T fallback) {
  return value && index < value->GetSize() ? (*value)[index] : fallback;
}

void requireCollectionSize(const char *name, std::size_t count,
                           std::size_t available) {
  if (count > available) {
    fail("Input Events tree has inconsistent " + std::string(name) +
         " collection size");
  }
}

void addOutputBranches(TTree &tree, OutputRow &row) {
  tree.Branch("run", &row.run, "run/I");
  tree.Branch("lumi", &row.lumi, "lumi/I");
  tree.Branch("event", &row.event, "event/L");
  tree.Branch("nPV", &row.nPV, "nPV/I");
  tree.Branch("pvX", &row.pvX, "pvX/F");
  tree.Branch("pvY", &row.pvY, "pvY/F");
  tree.Branch("pvZ", &row.pvZ, "pvZ/F");
  tree.Branch("pvChi2", &row.pvChi2, "pvChi2/F");
  tree.Branch("pvNdof", &row.pvNdof, "pvNdof/F");
  tree.Branch("pvScore", &row.pvScore, "pvScore/F");
  tree.Branch("mll", &row.mll, "mll/F");
  tree.Branch("ptll", &row.ptll, "ptll/F");
  tree.Branch("etall", &row.etall, "etall/F");
  tree.Branch("phill", &row.phill, "phill/F");
  tree.Branch("yll", &row.yll, "yll/F");
  tree.Branch("finalState", &row.finalState, "finalState/I");
  tree.Branch("trigBits", &row.trigBits, "trigBits/I");

  for (std::size_t index = 0; index < row.leptons.size(); ++index) {
    const std::string prefix = "l" + std::to_string(index + 1);
    OutputLepton &lepton = row.leptons[index];
    tree.Branch((prefix + "pt").c_str(), &lepton.pt, (prefix + "pt/F").c_str());
    tree.Branch((prefix + "eta").c_str(), &lepton.eta,
                (prefix + "eta/F").c_str());
    tree.Branch((prefix + "phi").c_str(), &lepton.phi,
                (prefix + "phi/F").c_str());
    tree.Branch((prefix + "mass").c_str(), &lepton.mass,
                (prefix + "mass/F").c_str());
    tree.Branch((prefix + "charge").c_str(), &lepton.charge,
                (prefix + "charge/I").c_str());
    tree.Branch((prefix + "pdgId").c_str(), &lepton.pdgId,
                (prefix + "pdgId/I").c_str());
    tree.Branch((prefix + "zId").c_str(), &lepton.zId,
                (prefix + "zId/I").c_str());
    tree.Branch((prefix + "dxy").c_str(), &lepton.dxy,
                (prefix + "dxy/F").c_str());
    tree.Branch((prefix + "dz").c_str(), &lepton.dz, (prefix + "dz/F").c_str());
    tree.Branch((prefix + "dxyErr").c_str(), &lepton.dxyErr,
                (prefix + "dxyErr/F").c_str());
    tree.Branch((prefix + "dzErr").c_str(), &lepton.dzErr,
                (prefix + "dzErr/F").c_str());
    tree.Branch((prefix + "sip3d").c_str(), &lepton.sip3d,
                (prefix + "sip3d/F").c_str());
    tree.Branch((prefix + "ip3d").c_str(), &lepton.ip3d,
                (prefix + "ip3d/F").c_str());
    tree.Branch((prefix + "pfRelIso03").c_str(), &lepton.pfRelIso03,
                (prefix + "pfRelIso03/F").c_str());
    tree.Branch((prefix + "miniRelIso").c_str(), &lepton.miniRelIso,
                (prefix + "miniRelIso/F").c_str());
    tree.Branch((prefix + "muMedium").c_str(), &lepton.muMedium,
                (prefix + "muMedium/I").c_str());
    tree.Branch((prefix + "muTight").c_str(), &lepton.muTight,
                (prefix + "muTight/I").c_str());
    tree.Branch((prefix + "muGlobal").c_str(), &lepton.muGlobal,
                (prefix + "muGlobal/I").c_str());
    tree.Branch((prefix + "muPF").c_str(), &lepton.muPF,
                (prefix + "muPF/I").c_str());
    tree.Branch((prefix + "elCutBased").c_str(), &lepton.elCutBased,
                (prefix + "elCutBased/I").c_str());
    tree.Branch((prefix + "elMvaWP90").c_str(), &lepton.elMvaWP90,
                (prefix + "elMvaWP90/I").c_str());
    tree.Branch((prefix + "elMvaWP80").c_str(), &lepton.elMvaWP80,
                (prefix + "elMvaWP80/I").c_str());
    tree.Branch((prefix + "elMvaBdt").c_str(), &lepton.elMvaBdt,
                (prefix + "elMvaBdt/F").c_str());
  }
}

void fillOutputLepton(OutputLepton &output, const Lepton &input) {
  output.pt = static_cast<Float_t>(input.pt);
  output.eta = static_cast<Float_t>(input.eta);
  output.phi = static_cast<Float_t>(input.phi);
  output.mass = static_cast<Float_t>(input.mass);
  output.charge = input.charge;
  output.pdgId = input.pdgId;
  output.zId = 1;
  output.dxy = static_cast<Float_t>(input.dxy);
  output.dz = static_cast<Float_t>(input.dz);
  output.dxyErr = static_cast<Float_t>(input.dxyErr);
  output.dzErr = static_cast<Float_t>(input.dzErr);
  output.sip3d = static_cast<Float_t>(input.sip3d);
  output.ip3d = static_cast<Float_t>(input.ip3d);
  output.pfRelIso03 = static_cast<Float_t>(input.pfRelIso03);
  output.miniRelIso = static_cast<Float_t>(input.miniRelIso);
  output.muMedium = input.muMedium;
  output.muTight = input.muTight;
  output.muGlobal = input.muGlobal;
  output.muPF = input.muPF;
  output.elCutBased = input.elCutBased;
  output.elMvaWP90 = input.elMvaWP90;
  output.elMvaWP80 = input.elMvaWP80;
  output.elMvaBdt = static_cast<Float_t>(input.elMvaBdt);
}

int run(const Options &options) {
  std::unique_ptr<TFile> input(TFile::Open(options.input.c_str(), "READ"));
  if (!input || input->IsZombie())
    fail("Cannot open input ROOT file: " + options.input);
  TTree *events = dynamic_cast<TTree *>(input->Get("Events"));
  if (events == nullptr)
    fail("Input ROOT file is missing the Events tree: " + options.input);
  if (events->GetBranch("nMuon") == nullptr) {
    fail("Input Events tree lacks reconstructed nMuon");
  }
  std::unique_ptr<TFile> output(
      TFile::Open(options.output.c_str(), "RECREATE"));
  if (!output || output->IsZombie())
    fail("Cannot create output ROOT file: " + options.output);

  OutputRow row;
  TTree dilepTree("dilepTree", "Best dilepton control-region candidates "
                               "selected from reconstructed NanoAOD");
  addOutputBranches(dilepTree, row);

  TTreeReader reader(events);
  auto runNumber = requiredValuePtr<UInt_t>(reader, *events, "run");
  auto lumiBlock = requiredValuePtr<UInt_t>(reader, *events, "luminosityBlock");
  auto eventNumber = requiredValuePtr<ULong64_t>(reader, *events, "event");

  auto pvNpv = optionalValue<Int_t>(reader, *events, "PV_npvs");
  auto pvNpvGood = optionalValue<Int_t>(reader, *events, "PV_npvsGood");
  auto pvX = optionalValue<Float_t>(reader, *events, "PV_x");
  auto pvY = optionalValue<Float_t>(reader, *events, "PV_y");
  auto pvZ = optionalValue<Float_t>(reader, *events, "PV_z");
  auto pvChi2 = optionalValue<Float_t>(reader, *events, "PV_chi2");
  auto pvNdof = optionalValue<Float_t>(reader, *events, "PV_ndof");
  auto pvScore = optionalValue<Float_t>(reader, *events, "PV_score");

  std::array<std::unique_ptr<TTreeReaderValue<Bool_t>>, kHltPaths.size()>
      hltBits;
  for (std::size_t bit = 0; bit < hltBits.size(); ++bit) {
    hltBits[bit] =
        optionalValue<Bool_t>(reader, *events, triggerBranchName(bit).c_str());
  }

  auto nMuon = optionalValue<UInt_t>(reader, *events, "nMuon");
  auto muonPt = optionalArray<Float_t>(reader, *events, "Muon_pt");
  auto muonEta = optionalArray<Float_t>(reader, *events, "Muon_eta");
  auto muonPhi = optionalArray<Float_t>(reader, *events, "Muon_phi");
  auto muonMass = optionalArray<Float_t>(reader, *events, "Muon_mass");
  auto muonCharge = optionalArray<Int_t>(reader, *events, "Muon_charge");
  auto muonDxy = optionalArray<Float_t>(reader, *events, "Muon_dxy");
  auto muonDz = optionalArray<Float_t>(reader, *events, "Muon_dz");
  auto muonDxyErr = optionalArray<Float_t>(reader, *events, "Muon_dxyErr");
  auto muonDzErr = optionalArray<Float_t>(reader, *events, "Muon_dzErr");
  auto muonSip3d = optionalArray<Float_t>(reader, *events, "Muon_sip3d");
  auto muonIp3d = optionalArray<Float_t>(reader, *events, "Muon_ip3d");
  auto muonIso = optionalArray<Float_t>(reader, *events, "Muon_pfRelIso03_all");
  auto muonMiniIso =
      optionalArray<Float_t>(reader, *events, "Muon_miniPFRelIso_all");
  auto muonMedium = optionalArray<Bool_t>(reader, *events, "Muon_mediumId");
  auto muonTight = optionalArray<Bool_t>(reader, *events, "Muon_tightId");
  auto muonPf = optionalArray<Bool_t>(reader, *events, "Muon_isPFcand");
  auto muonGlobal = optionalArray<Bool_t>(reader, *events, "Muon_isGlobal");

  auto nElectron = optionalValue<UInt_t>(reader, *events, "nElectron");
  auto electronPt = optionalArray<Float_t>(reader, *events, "Electron_pt");
  auto electronEta = optionalArray<Float_t>(reader, *events, "Electron_eta");
  auto electronPhi = optionalArray<Float_t>(reader, *events, "Electron_phi");
  auto electronMass = optionalArray<Float_t>(reader, *events, "Electron_mass");
  auto electronCharge =
      optionalArray<Int_t>(reader, *events, "Electron_charge");
  auto electronDxy = optionalArray<Float_t>(reader, *events, "Electron_dxy");
  auto electronDz = optionalArray<Float_t>(reader, *events, "Electron_dz");
  auto electronDxyErr =
      optionalArray<Float_t>(reader, *events, "Electron_dxyErr");
  auto electronDzErr =
      optionalArray<Float_t>(reader, *events, "Electron_dzErr");
  auto electronSip3d =
      optionalArray<Float_t>(reader, *events, "Electron_sip3d");
  auto electronIp3d = optionalArray<Float_t>(reader, *events, "Electron_ip3d");
  auto electronIso =
      optionalArray<Float_t>(reader, *events, "Electron_pfRelIso03_all");
  auto electronMiniIso =
      optionalArray<Float_t>(reader, *events, "Electron_miniPFRelIso_all");
  auto electronCutBased =
      optionalArray<Int_t>(reader, *events, "Electron_cutBased");
  auto electronMvaFall17Wp90 =
      optionalArray<Bool_t>(reader, *events, "Electron_mvaFall17V2Iso_WP90");
  auto electronMvaFall17Wp80 =
      optionalArray<Bool_t>(reader, *events, "Electron_mvaFall17V2Iso_WP80");
  auto electronMvaFall17 =
      optionalArray<Float_t>(reader, *events, "Electron_mvaFall17V2Iso");
  auto electronMvaIsoWp90 =
      optionalArray<Bool_t>(reader, *events, "Electron_mvaIso_WP90");
  auto electronMvaIsoWp80 =
      optionalArray<Bool_t>(reader, *events, "Electron_mvaIso_WP80");
  auto electronMvaIso =
      optionalArray<Float_t>(reader, *events, "Electron_mvaIso");

  if (!nMuon || !muonPt || !muonEta || !muonPhi || !muonMass || !muonCharge ||
      !muonDxy || !muonDz || !muonDxyErr || !muonDzErr || !muonSip3d ||
      !muonIso || !muonMedium || !muonTight || !muonPf || !muonGlobal ||
      !nElectron || !electronPt || !electronEta || !electronPhi ||
      !electronMass || !electronCharge || !electronDxy || !electronDz ||
      !electronDxyErr || !electronDzErr || !electronSip3d || !electronIso) {
    fail("Input Events tree is missing required reconstructed lepton branches");
  }

  // Loose reconstructed leptons: muons require pT >= 5 GeV, |eta| <= 2.4,
  // |dxy| <= 0.5 cm, |dz| <= 1 cm, sip3d <= 4, isolation <= 0.35, and
  // either the PF-candidate or global-muon flag.  Electrons require pT >= 7
  // GeV, |eta| <= 2.5, |dxy| <= 0.5 cm, |dz| <= 1 cm, sip3d <= 4, and
  // isolation <= 0.35.
  while (reader.Next()) {
    std::vector<Lepton> leptons;

    const std::size_t muonCount = **nMuon;
    requireCollectionSize("Muon", muonCount, muonPt->GetSize());
    requireCollectionSize("Muon", muonCount, muonEta->GetSize());
    requireCollectionSize("Muon", muonCount, muonPhi->GetSize());
    requireCollectionSize("Muon", muonCount, muonMass->GetSize());
    requireCollectionSize("Muon", muonCount, muonCharge->GetSize());
    requireCollectionSize("Muon", muonCount, muonDxy->GetSize());
    requireCollectionSize("Muon", muonCount, muonDz->GetSize());
    requireCollectionSize("Muon", muonCount, muonDxyErr->GetSize());
    requireCollectionSize("Muon", muonCount, muonDzErr->GetSize());
    requireCollectionSize("Muon", muonCount, muonSip3d->GetSize());
    requireCollectionSize("Muon", muonCount, muonIso->GetSize());
    requireCollectionSize("Muon", muonCount, muonMedium->GetSize());
    requireCollectionSize("Muon", muonCount, muonTight->GetSize());
    requireCollectionSize("Muon", muonCount, muonPf->GetSize());
    requireCollectionSize("Muon", muonCount, muonGlobal->GetSize());
    for (std::size_t index = 0; index < muonCount; ++index) {
      const double pt = (*muonPt)[index];
      const double eta = (*muonEta)[index];
      if (pt < 5.0 || std::abs(eta) > 2.4)
        continue;
      const double dxy = (*muonDxy)[index];
      const double dz = (*muonDz)[index];
      const double sip3d = (*muonSip3d)[index];
      const double iso = (*muonIso)[index];
      const bool isPf = (*muonPf)[index];
      const bool isGlobal = (*muonGlobal)[index];
      if (std::abs(dxy) > 0.5 || std::abs(dz) > 1.0 || sip3d > 4.0 ||
          iso > 0.35 || (!isPf && !isGlobal))
        continue;

      Lepton lepton;
      lepton.pt = pt;
      lepton.eta = eta;
      lepton.phi = (*muonPhi)[index];
      lepton.mass = arrayOr(muonMass, index, 0.10566F);
      lepton.charge = (*muonCharge)[index];
      lepton.pdgId = -13 * lepton.charge;
      lepton.dxy = dxy;
      lepton.dz = dz;
      lepton.dxyErr = (*muonDxyErr)[index];
      lepton.dzErr = (*muonDzErr)[index];
      lepton.sip3d = sip3d;
      lepton.ip3d = arrayOr(muonIp3d, index, -99.0F);
      lepton.pfRelIso03 = iso;
      lepton.miniRelIso = arrayOr(muonMiniIso, index, -99.0F);
      lepton.muMedium = (*muonMedium)[index];
      lepton.muTight = (*muonTight)[index];
      lepton.muGlobal = isGlobal;
      lepton.muPF = isPf;
      leptons.push_back(lepton);
    }

    const std::size_t electronCount = **nElectron;
    requireCollectionSize("Electron", electronCount, electronPt->GetSize());
    requireCollectionSize("Electron", electronCount, electronEta->GetSize());
    requireCollectionSize("Electron", electronCount, electronPhi->GetSize());
    requireCollectionSize("Electron", electronCount, electronMass->GetSize());
    requireCollectionSize("Electron", electronCount, electronCharge->GetSize());
    requireCollectionSize("Electron", electronCount, electronDxy->GetSize());
    requireCollectionSize("Electron", electronCount, electronDz->GetSize());
    requireCollectionSize("Electron", electronCount, electronDxyErr->GetSize());
    requireCollectionSize("Electron", electronCount, electronDzErr->GetSize());
    requireCollectionSize("Electron", electronCount, electronSip3d->GetSize());
    requireCollectionSize("Electron", electronCount, electronIso->GetSize());
    for (std::size_t index = 0; index < electronCount; ++index) {
      const double pt = (*electronPt)[index];
      const double eta = (*electronEta)[index];
      if (pt < 7.0 || std::abs(eta) > 2.5)
        continue;
      const double dxy = (*electronDxy)[index];
      const double dz = (*electronDz)[index];
      const double sip3d = (*electronSip3d)[index];
      const double iso = (*electronIso)[index];
      if (std::abs(dxy) > 0.5 || std::abs(dz) > 1.0 || sip3d > 4.0 ||
          iso > 0.35)
        continue;

      const bool wp90Fall17 =
          arrayOr(electronMvaFall17Wp90, index, static_cast<Bool_t>(false));
      const bool wp80Fall17 =
          arrayOr(electronMvaFall17Wp80, index, static_cast<Bool_t>(false));
      const bool wp90Iso =
          arrayOr(electronMvaIsoWp90, index, static_cast<Bool_t>(false));
      const bool wp80Iso =
          arrayOr(electronMvaIsoWp80, index, static_cast<Bool_t>(false));
      Lepton lepton;
      lepton.pt = pt;
      lepton.eta = eta;
      lepton.phi = (*electronPhi)[index];
      lepton.mass = arrayOr(electronMass, index, 0.000511F);
      lepton.charge = (*electronCharge)[index];
      lepton.pdgId = -11 * lepton.charge;
      lepton.dxy = dxy;
      lepton.dz = dz;
      lepton.dxyErr = (*electronDxyErr)[index];
      lepton.dzErr = (*electronDzErr)[index];
      lepton.sip3d = sip3d;
      lepton.ip3d = arrayOr(electronIp3d, index, -99.0F);
      lepton.pfRelIso03 = iso;
      lepton.miniRelIso = arrayOr(electronMiniIso, index, -99.0F);
      lepton.elCutBased = arrayOr(electronCutBased, index, 0);
      lepton.elMvaWP90 = wp90Fall17 || wp90Iso;
      lepton.elMvaWP80 = wp80Fall17 || wp80Iso;
      lepton.elMvaBdt = electronMvaFall17
                            ? arrayOr(electronMvaFall17, index, -99.0F)
                            : arrayOr(electronMvaIso, index, -99.0F);
      leptons.push_back(lepton);
    }

    if (leptons.size() < 2)
      continue;
    std::stable_sort(leptons.begin(), leptons.end(),
                     [](const Lepton &left, const Lepton &right) {
                       return left.pt > right.pt;
                     });

    bool found = false;
    const Candidate candidate = findBestPair(leptons, found);
    if (!found)
      continue;

    row.run = static_cast<Int_t>(**runNumber);
    row.lumi = static_cast<Int_t>(**lumiBlock);
    row.event = static_cast<Long64_t>(**eventNumber);
    const Int_t npv = valueOr(pvNpv, 0);
    row.nPV = npv != 0 ? npv : valueOr(pvNpvGood, 1);
    row.pvX = valueOr(pvX, 0.0F);
    row.pvY = valueOr(pvY, 0.0F);
    row.pvZ = valueOr(pvZ, 0.0F);
    row.pvChi2 = valueOr(pvChi2, 0.0F);
    row.pvNdof = valueOr(pvNdof, 0.0F);
    row.pvScore = valueOr(pvScore, 0.0F);
    row.mll = static_cast<Float_t>(candidate.mll);
    row.ptll = static_cast<Float_t>(candidate.ptll);
    row.etall = static_cast<Float_t>(candidate.etall);
    row.phill = static_cast<Float_t>(candidate.phill);
    row.yll = static_cast<Float_t>(candidate.yll);
    row.finalState = candidate.finalState;
    row.trigBits = 0;
    for (std::size_t bit = 0; bit < hltBits.size(); ++bit) {
      if (hltBits[bit] && **hltBits[bit])
        row.trigBits |= (1 << bit);
    }
    for (std::size_t index = 0; index < candidate.pair.size(); ++index) {
      fillOutputLepton(row.leptons[index], candidate.pair[index]);
    }
    dilepTree.Fill();
  }

  if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd &&
      reader.GetEntryStatus() != TTreeReader::kEntryValid) {
    fail("Error while reading Events tree from " + options.input);
  }

  output->cd();
  dilepTree.Write();
  output->Close();
  return 0;
}

} // namespace

int main(int argc, char **argv) {
  try {
    return run(parseOptions(argc, argv));
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
}
