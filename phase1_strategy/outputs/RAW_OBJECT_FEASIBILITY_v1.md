# Raw-object VBF/FSR feasibility record v1

## Scope and conclusion

This is a feasibility record for the binding reconstructed-jet VBF category
and FSR recovery. It is not a proxy category, a production selection, or a
claim that the supplied flat fake-data file can already be re-ntuplized.

**MC conclusion — demonstrated.** The local production area contains matching
raw NanoAOD for every supplied fit-MC process. A matched VBF NanoAOD pilot was
actually re-ntuplized with `uproot`, preserving event identity, cleaned jets,
DeepJet b-tag scores, `FsrPhoton` candidates, and generic photons.

**Data conclusion — unresolved.** The supplied `fake_data_10fb.root` and
`fake_data_dilep_10fb.root` are flat products. No raw-data NanoAOD location or
recoverable raw-data manifest is present under the supplied production root or
its immediately relevant owner area. The VBF/FSR programme is therefore ready
for MC development but must not be used in the final fake-data fit until the
data provenance is supplied and passes the same event-match pilot.

## Local provenance investigation

The following non-destructive checks were run on 2026-07-15.

| Check | Concrete result | Interpretation |
|---|---|---|
| Production-root directory inventory | Found `nanoaod_modified`, flat MC products, and `fake_data`; no raw-data NanoAOD directory. | Raw MC is locally available; raw data is not exposed in this production root. |
| Owner-area NanoAOD/H4l directory inventory | The only H4l NanoAOD directory is `h4l/.../nanoaod_modified`. | A broader owner-area search found no alternate data NanoAOD path. |
| `uproot` schema of supplied fake data | The H4l file has `h4lTree;1` and `Metadata;1`; the dilepton file has `dilepTree` and `Metadata`; both retain run/lumi/event but have zero `Jet_*`, `Photon_*`, or `FsrPhoton_*` branches. | Event identifiers alone do not reconstruct missing raw objects. |
| Mixer provenance | The available dilepton mixer summary names `/pool/condor/dir_933352/modified_dilep_mixer_manifest.json`; that path is absent locally, and the summary has no raw-data input manifest. | The available record cannot recover a data raw source. |

No proxy based on $p_\mathrm{T}^{4\ell}$, candidate multiplicity, or flat-tree
lepton information is allowed in place of the VBF category.

## MC raw-input coverage

One representative NanoAOD file was opened for each supplied fit-MC process.
Every one had an `Events` tree and all of the following: `run`,
`luminosityBlock`, `event`; `Jet_pt`, `Jet_eta`, `Jet_phi`, `Jet_mass`,
`Jet_jetId`, `Jet_rawFactor`, `Jet_btagDeepFlavB`; `FsrPhoton_pt`,
`FsrPhoton_eta`, `FsrPhoton_phi`, `FsrPhoton_relIso03`,
`FsrPhoton_muonIdx`; and `Photon_pt`, `Photon_eta`, `Photon_phi`,
`Photon_pfRelIso03_all`.

| Process | Representative raw input | Entries | Missing required branches |
|---|---|---:|---|
| DYJetsToLL | `nanoaod_modified/DYJetsToLL/130000/04935ACA-..._scale_smear.root` | 2,083,510 | none |
| GGZZ2E2Mu | `.../GGZZ2E2Mu/110000/1755CC95-..._scale_smear.root` | 2,000 | none |
| GGZZ4E | `.../GGZZ4E/230000/63E0F24B-..._scale_smear.root` | 3,000 | none |
| GGZZ4Mu | `.../GGZZ4Mu/110000/07AA19D1-..._scale_smear.root` | 2,920 | none |
| GluGluToHToZZ | `.../GluGluToHToZZ/2430000/0D7D1623-..._scale_smear.root` | 5,790 | none |
| TTBar | `.../TTBar/2820000/25C8DD83-..._scale_smear.root` | 880,778 | none |
| VBF_HToZZ | `.../VBF_HToZZ/230000/3602787B-..._scale_smear.root` | 8,000 | none |
| WMHToZZ | `.../WMHToZZ/100000/5FCCB7BB-..._scale_smear.root` | 9,228 | none |
| WPHToZZ | `.../WPHToZZ/100000/D6D7974E-..._scale_smear.root` | 24,679 | none |
| ZHToZZ | `.../ZHToZZ/260000/1203829B-..._scale_smear.root` | 512,586 | none |
| ZZTo4L | `.../ZZTo4L/110000/02054071-..._scale_smear.root` | 660,000 | none |

The abbreviated paths are under
`/eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/cms_opendata_2017_full_production/nanoaod_modified/`.

## Actual matched VBF pilot

The reproducible task is:

```bash
pixi run raw-object-pilot
```

It reads the VBF raw file
`nanoaod_modified/VBF_HToZZ/230000/3602787B-E4E1-784B-95EE-663C48DB1CE6_scale_smear.root`
and its UUID-matched flat output
`h4l_mc_modified/VBF_HToZZ/230000/3602787B-E4E1-784B-95EE-663C48DB1CE6_scale_smear_h4l.root`.
The persisted output is [raw_object_pilot_vbf_v1.root](raw_object_pilot_vbf_v1.root),
with the machine-readable summary
[raw_object_pilot_vbf_v1.json](raw_object_pilot_vbf_v1.json).

| Test | Result |
|---|---:|
| Raw `Events` entries | 8,000 |
| Existing flat candidates | 1,198 |
| Exact run/lumi/event matches | 1,198 / 1,198 (100%) |
| Matched candidates in $105<m_{4\ell}<140$ GeV | 1,134 |
| Matched candidates with at least two cleaned raw jets | 711 |
| In-window candidates with at least two cleaned raw jets | 671 |
| Matched candidates with `FsrPhoton` | 37 |
| Matched candidates with a muon-indexed FSR association | 37 |
| Required jet/b-tag/FSR/photon branches missing | none |

For this object-preservation test, jets satisfy $p_\mathrm{T}>30$ GeV,
$|\eta|<4.7$, `Jet_jetId>=2`, and $\Delta R(\mathrm{jet},\ell)>0.4$ to the
four preserved candidate leptons. This is an object baseline only; it is not
an estimate of VBF purity or an optimized VBF tag. Each pilot row contains the
selected-jet four-vector, `Jet_rawFactor`, `Jet_btagDeepFlavB`, leading-dijet
mass, all `FsrPhoton` fields, and generic photon fields. Thus the required
objects are actually carried forward, rather than merely observed in a schema.

## Binding production route and gates

1. Re-run the raw-object producer for each nominal and systematic MC input,
   preserving run/lumi/event, selected jets, raw/JEC/JER-relevant fields,
   b tags, FSR photons, generic photons, and all four calibrated leptons.
2. Re-run the pilot after a raw-data path is supplied. It must establish a
   one-to-one event match to `fake_data_10fb.root`; otherwise the final data
   VBF/FSR category is blocked rather than approximated.
3. Form candidate masses only after FSR association/recovery and calibrated
   lepton four-vectors. Muon FSR uses `FsrPhoton_muonIdx`; electron FSR is
   tested from retained generic photons and electron association/geometry.
4. Use the HIG-16-041 reference-like two-jet discriminant as the final tag
   after its implementation is independently checked. The raw-object pilot is
   deliberately not substituted for that tag.
5. On matched MC, report before/after FSR/calibration pairing agreement and
   distributions of all five angular-NN inputs. Any pairing or angle change is
   quantified by final state and propagated to NN templates before fitting.
6. Before unblinding, publish a per-process/final-state/category expected-yield
   table, VBF signal efficiency, purity, background composition, JEC/JER and
   b-tag migration, and the fraction of VBF-bin Fisher information that remains
   after nuisance profiling. A category with no non-prior-dominated information
   requires a formal strategy revision, not a proxy.

## Status

- **Resolved for MC feasibility:** actual raw input, branch coverage, exact
  event matching, and a persisted jet/FSR re-ntupling pilot.
- **Unresolved for supplied fake data:** raw data location/manifest and
  data-to-raw event matching. This is a hard prerequisite for the final
  data VBF/FSR fit and is tracked in `COMMITMENTS.md` as open, not downscoped.
