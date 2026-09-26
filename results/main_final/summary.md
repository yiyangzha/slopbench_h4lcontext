# H -> ZZ* -> 4l, UL16 PFNano pseudo-data, 20.0 fb^-1 (results v5/r2)

## Benchmark results (set b; RESULT.json)

* mu = 0.903 +0.189/-0.169 (stat +0.180/-0.163, syst 0.051); expected (Asimov, pre-fit) +0.193/-0.173
* m_H = 124.964 +0.396/-0.388 (stat +0.392/-0.385, syst 0.051) GeV; expected +-0.303 (pre-fit), +-0.319 (post-fit)
* local significance 8.46 (expected 10.08); discovery True
* saturated GoF (toys): p = 0.422 (toys: saturated, binned per channel in m4l (5 GeV) x D_bkg^kin template bins, toys from the best fit refitted with the full unbinned likelihood; q_saturated 190.3 on 735 bins)
* coverage of the 68 % interval of mu (toys): 0.647

## Paper-style results (set a)

* mu = 0.890 +0.200/-0.175 (stat +0.178/-0.161, syst 0.076)
* m_H = 124.960 +0.397/-0.390 (stat +0.394/-0.387, syst 0.052) GeV
* mu at m_H = 125.09 GeV (set b): 0.884 +0.200/-0.176 (stat +0.178/-0.162, syst 0.078)

## Calibration and efficiencies (RESULT.json: data-weighted averages over the Z -> ll calibration-sample data leptons)

| flavour | scale_shift | smear | sel_eff |
|---|---|---|---|
| muon | -0.01993 +- 0.00037 | 0.0116 +- 0.0005 | 0.9941 +- 0.0018 |
| electron | -0.01987 +- 0.00049 | 0.0113 +- 0.0021 | 1.0004 +- 0.0031 |

At the reference point pT = 45 GeV, |eta| = 1.2:

| flavour | scale shift | smear | efficiency SF |
|---|---|---|---|
| muon | -0.02030 +- 0.00037 | 0.0109 +- 0.0006 | 0.9941 +- 0.0002 |
| electron | -0.02017 +- 0.00050 | 0.0097 +- 0.0027 | 1.0005 +- 0.0005 |

eff_correction: 4mu 0.9871, 4e 1.0079, 2e2mu 0.9975, inclusive 0.9959

Reconstruction efficiency from the MC (user decision); cross-check, Z -> ll yield ratio data/MC against the T&P SF(full)^2: muon 0.9822 vs 0.9875 (remaining 0.9973 per lepton); electron 0.9686 vs 1.0001 (remaining 0.9841 per lepton)

## Mass fits (Table 6)

| fit | observed m_H [GeV] | expected uncertainty |
|---|---|---|
| 1D_refit | 124.557 +0.358/-0.332 (stat +0.354/-0.330, syst 0.045) | +-0.340 |
| 1D_norefit | 124.809 +0.485/-0.463 (stat +0.482/-0.461, syst 0.049) | +-0.407 |
| 2D_refit | 124.941 +0.391/-0.381 (stat +0.387/-0.380, syst 0.048) | +-0.310 |
| 2D_norefit | 125.028 +0.412/-0.407 (stat +0.409/-0.405, syst 0.042) | +-0.359 |
| 3D_refit | 124.964 +0.396/-0.388 (stat +0.392/-0.385, syst 0.051) | +-0.303 |
| 3D_norefit | 125.037 +0.409/-0.441 (stat +0.406/-0.439, syst 0.045) | +-0.349 |
| 2D_Dkin_refit | 124.562 +0.371/-0.342 (stat +0.368/-0.339, syst 0.046) | +-0.330 |
| 2D_Dkin_norefit | 124.802 +0.475/-0.489 (stat +0.472/-0.486, syst 0.050) | +-0.394 |

## Per final state

| final state | mu | m_H [GeV] (separate fit) |
|---|---|---|
| 4mu | 1.035 +0.356/-0.296 (stat +0.333/-0.285, syst 0.102) | 124.598 +0.421/-0.397 (stat +0.418/-0.394, syst 0.050) |
| 4e | 0.627 +0.481/-0.378 (stat +0.456/-0.365, syst 0.123) | 129.168 +0.794/-1.011 (stat +0.790/-1.009, syst 0.076) |
| 2e2mu | 0.850 +0.269/-0.229 (stat +0.251/-0.218, syst 0.078) | 124.637 +0.751/-0.864 (stat +0.750/-0.863, syst 0.040) |

m_H compatibility of the final states: q = 11.87 (2 dof), p = 0.003

## Signal strengths per category, mode, (mu_F, mu_V) and STXS stage 0

Paper style (Eq. 10.1): L(m4l, D_bkg^kin) at m_H = 125.09 GeV, POIs in [0, 20]; STXS stage 0 equals the per-mode strengths here (every reconstructed signal event has |y_H| < 2.5).

* mu_category: mu_VBF-2jet 1.15 +0.75/-0.53 (stat +0.74/-0.53, syst 0.07); mu_VH-hadronic 0.00 +0.54/-0.00 (stat +0.54/-0.00, syst 0.00); mu_VH-leptonic 0.00 +1.78/-0.00 (stat +1.74/-0.00, syst 0.09); mu_ttH 0.00 +4.19/-0.00 (stat +4.18/-0.00, syst 0.00); mu_VH-MET 0.00 +6.14/-0.00 (stat +6.14/-0.00, syst 0.00); mu_VBF-1jet 0.16 +0.28/-0.15 (stat +0.27/-0.15, syst 0.02); mu_Untagged 1.08 +0.25/-0.22 (stat +0.22/-0.20, syst 0.10)
* mu_mode: mu_ggH 0.92 +0.25/-0.22 (stat +0.22/-0.20, syst 0.09); mu_VBF 0.92 +1.39/-0.92 (stat +1.38/-0.92, syst 0.05); mu_VHhad 0.00 +1.82/-0.00 (stat +1.81/-0.00, syst 0.00); mu_VHlep 0.00 +1.94/-0.00 (stat +1.91/-0.00, syst 0.00)
* mu_fv: mu_F 0.96 +0.26/-0.23 (stat +0.23/-0.21, syst 0.09); mu_V 0.35 +1.05/-0.35 (stat +1.04/-0.35, syst 0.09)
* stxs0: r_ggH 0.92 +0.25/-0.22 (stat +0.22/-0.20, syst 0.09); r_VBF 0.92 +1.39/-0.92 (stat +1.38/-0.92, syst 0.05); r_VHhad 0.00 +1.82/-0.00 (stat +1.81/-0.00, syst 0.00); r_VHlep 0.00 +1.94/-0.00 (stat +1.91/-0.00, syst 0.00)

## Fiducial cross sections

* integrated (m_H 125.09): sigma_fid = 2.502 fb (r = 0.783 +0.195/-0.174 (stat +0.183/-0.165, syst 0.058)); SM 3.194 fb
* integrated (m_H profiled): sigma_fid = 2.543 fb
* per final state: 4mu 0.958 fb; 4e 0.544 fb; 2e2mu 1.002 fb; compatibility p = 0.465
* Table 5 (A_fid relative to sigma_eff of the delivered LO samples): ggH A_fid 0.474, eps 0.667, f_nonfid 0.045; VBF A_fid 0.487, eps 0.672, f_nonfid 0.049; VH A_fid 0.374, eps 0.654, f_nonfid 0.178
* pt4l: bin 0 0.678 (SM 0.679) fb; bin 1 0.513 (SM 0.666) fb; bin 2 0.801 (SM 1.125) fb; bin 3 0.427 (SM 0.571) fb; bin 4 0.200 (SM 0.152) fb
* njets: bin 0 1.886 (SM 1.853) fb; bin 1 0.566 (SM 0.921) fb; bin 2 0.251 (SM 0.319) fb; bin 3 0.000 (SM 0.100) fb
* ptj1: bin 0 1.827 (SM 1.853) fb; bin 1 0.612 (SM 0.507) fb; bin 2 0.000 (SM 0.414) fb; bin 3 0.225 (SM 0.317) fb; bin 4 0.077 (SM 0.102) fb

## Width and Z -> 4l

* Gamma_H < 7.05 GeV (95 % CL), expected < 2.15 GeV
* m_Z(4l) = 91.225 +- 0.182 (stat) +- 0.035 (syst) GeV

## Validation

* toys r2_toys_mu (500 per point): mu0_mh125: mu median -0.060 pull -1.16/2.45, m_H median 125.178 pull -0.88/12.83, coverage mu 0.31, Z>3 0.002; mu0.5_mh125: mu median 0.499 pull -0.00/1.01, m_H median 125.000 pull +0.01/1.04, coverage mu 0.69, Z>3 0.980; mu1_mh125: mu median 0.994 pull -0.05/1.04, m_H median 124.987 pull -0.03/1.02, coverage mu 0.65, Z>3 1.000; mu2_mh125: mu median 2.004 pull +0.01/0.98, m_H median 124.986 pull -0.10/1.00, coverage mu 0.69, Z>3 1.000; mu3_mh125: mu median 2.976 pull -0.12/1.00, m_H median 124.994 pull -0.03/0.94, coverage mu 0.64, Z>3 1.000; slopes {'mu_at_mh125': 1.0074924282518052}
* toys r2_toys_mh (400 per point): mu1_mh121: mu median 0.990 pull -0.03/1.12, m_H median 120.976 pull -0.06/1.09, coverage mu 0.66, Z>3 0.997; mu1_mh123: mu median 1.006 pull -0.00/0.96, m_H median 123.016 pull -0.01/1.01, coverage mu 0.70, Z>3 1.000; mu1_mh127: mu median 1.012 pull +0.05/1.01, m_H median 126.962 pull -0.06/0.96, coverage mu 0.66, Z>3 1.000; mu1_mh129: mu median 1.000 pull -0.00/1.06, m_H median 128.979 pull -0.03/1.04, coverage mu 0.64, Z>3 1.000; slopes {'mH_at_mu1': 0.9975624598892879}
* paired injection: mu0_mh125 dmu 1.116 +- 0.018; mu1_mh125 dmu 0.999 +- 0.010 (the mu = 0 base has m_H undefined: m_H floats on the background, median mu -0.06)
* MC-event closure (truth mu 1, m_H 125): 1D_refit mu 1.0007 m_H 125.0052; 1D_norefit mu 0.9999 m_H 125.0021; 2Dmass_refit mu 1.0047 m_H 125.0131; 2Dmass_norefit mu 1.0026 m_H 125.0122; 3D_refit mu 1.0181 m_H 124.9820; 3D_norefit mu 1.0082 m_H 124.9377
* T&P efficiency closure (report point): k1 electron 0.99989 +- 0.00041 (injected 1.0); k1 muon 1.00003 +- 0.00022 (injected 1.0); k2 electron 0.94956 +- 0.00115 (injected 0.95); k2 muon 0.97031 +- 0.00028 (injected 0.97); k3 electron 0.89888 +- 0.00061 (injected 0.9); k3 muon 0.93038 +- 0.00028 (injected 0.93)

## Declared cuts (N-1)

| cut | threshold | role | Z with | Z without |
|---|---|---|---|---|
| mz1_window | 40 < m_Z1 < 120 GeV | sensitivity | 7.380 | 7.354 |
| mz2_window | 10 < m_Z2 < 120 GeV (signal region; smart cut with the same threshold) | sensitivity | 7.380 | 6.382 |
| iso_rel_fsr | < 0.35 (muons and electrons) | sensitivity | 7.380 | 5.139 |
| sip3d | < 4 | quality | 7.380 | 7.071 |
| min_pt_lep | > 5 GeV (muons), > 7 GeV (electrons) | quality | 7.380 | 6.437 |
| lead_sublead_20_10 | > 20 GeV and > 10 GeV | quality | 7.380 | 7.383 |
| os_pair_mass | m_ll > 4 GeV | quality | 7.380 | 7.390 |
| id_muon | PF muon (tracker high-pT ID above 200 GeV) | quality | 7.380 | 7.170 |
| id_electron | mvaFall17V2noIso WPL | quality | 7.380 | 4.313 |
| dxy_lep | < 0.5 cm | quality | 7.380 | null |
| dz_lep | < 1 cm | quality | 7.380 | null |
| eta_lep | < 2.4 (muons), < 2.5 (electrons) | quality | 7.380 | null |
| dr_lep | > 0.02 | quality | 7.380 | null |
| m4l_min | > 70 GeV | quality | 7.380 | 7.380 |
| trigger_or | fired | quality | 7.380 | null |
| ghost_cross_cleaning | removed | quality | 7.380 | null |

## Systematic coverage (groups; set (a))

* lepton_scale: yes (scale_mu and scale_e nuisances (per flavour, from the in-situ Z calibration and its closures) move the signal peak)
* lepton_resolution: yes (res_mu and res_e nuisances (per flavour) scale the signal width)
* lepton_efficiency: yes (eff_mu and eff_e nuisances from the tag-and-probe SF uncertainties (statistical and fit model) on the signal, qqZZ and ggZZ yields)
* luminosity: yes (set (a): 2.5 % on every MC yield (the benchmark set (b) has none: the pseudo-data luminosity is exact))
* branching_ratio: yes (set (a): 2 % on the signal (none in set (b)))
* reducible_bkg_Zjets: yes (zx_<final state> log-normal nuisances from the data-driven OS/SS combination (statistics, MC-closure systematic, method envelope))
* qqZZ_theory: yes (set (a): 0.045 on qqZZ (QCD scale and PDF))
* ggZZ_kfactor: yes (set (a): 10 % on ggZZ)
* signal_QCDscale: yes (set (a): per production mode {'ggH': 0.054, 'VBF': 0.004, 'VH': 0.017})
* signal_PDF_acceptance: yes (set (a): PDF + alpha_s per mode {'ggH': 0.032, 'VBF': 0.021, 'VH': 0.018}; the acceptance x efficiency m_H dependence through the morph nuisance)
* signal_shape: yes (the DCB peak position and width per flavour (scale and resolution nuisances) and the m_H morphing (scale vs shift); the DCB tails fixed from the MC fit)
* background_shape: partial (Z+X shape from the combined OS/SS fit (no shape nuisance), qqZZ/ggZZ m4l and D_bkg^kin shapes from the MC; no jet-energy variations (not available in the inputs))

Figures: plots/paper_comparison.png; the fit and scan plots under production_v3/inference/<select>/<model>/<label>/plots/.
