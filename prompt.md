# Analysis Title: H->4L measurment with CMS Open Data

Perform a complete H→4ℓ analysis with key observable distributions, a VBF category, an angular-variable neural-network discriminant in the four-lepton rest frame, and a simultaneous multi-category fit for mH and μ; use DY+jets MC as the fake-background model.

Compare every meaningfully comparable result with JHEP 11 (2017) 047 / CMS-HIG-16-041 and the PDG/world average where available, state why any result is not comparable, and use the same mass-fit window as that publication.

Use the supplied Z→ℓℓ control data to validate and calibrate the lepton momentum/energy scale, then propagate an uncertainty from this calibration to the fitted Higgs-boson mass.

Document the control-region check, the correction or residual treatment you adopt, and the resulting mH scale systematic.

h4l data and MC are produced by ./ntuplizer/h4l_ntuplizer.cpp, while Z control data and MC are produced by ./ntuplizer/dilep_ntuplizer.cpp.

Data: /eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/cms_opendata_2017_full_production/fake_data/fake_data_10fb.root (10 fb^-1)

Z control data: /eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/cms_opendata_2017_full_production/fake_data/fake_data_dilep_10fb.root (10 fb^-1)

MC: /eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/cms_opendata_2017_full_production/h4l_mc_nominal

Z control MC: /eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/cms_opendata_2017_full_production/dilep_mc_nominal

#GluGluToHToZZ.root cross-section:6.024e-03 pb     fullname:GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_minloHJJ_JHUGenV7011_pythia8
#VBF_HToZZ.root   cross-section:4.8794e-04 pb     fullname:VBF_HToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8
#ZHToZZ.root      cross-section:9.8394e-05 pb     fullname:ZH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8
#WPHToZZ.root     cross-section:1.072352e-04 pb  fullname:WplusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8
#WMHToZZ.root     cross-section:6.706e-05 pb     fullname:WminusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8
#ZZTo4L.root      cross-section:1.325e+00 pb     fullname:ZZTo4L_TuneCP5_13TeV_powheg_pythia8
#DYJetsToLL.root  cross-section:5.396e+03 pb     fullname:DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
#TTBar.root       cross-section:5.270e+01 pb     fullname:TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8
#GGZZ2E2Mu.root   cross-section:3.185e-03 pb     fullname:GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8
#GGZZ4Mu.root     cross-section:1.575e-03 pb     fullname:GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8
#GGZZ4E.root      cross-section:1.619e-03 pb     fullname:GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8

Note:
- These cross sections are the effective sample cross sections for the originally generated MC entries, including the decay or final-state definition implied by each sample name, and should be used to compute the MC event weight so that the simulated sample is normalized to the target data luminosity.
- For each MC ROOT file, the Metadata tree records the number of originally generated entries; the total generated event count should be obtained by summing nEvents over all entries in the Metadata tree.
