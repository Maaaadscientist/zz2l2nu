# Data files

This directory contains auxiliary files with various corrections. They are organized into several subdirectories:


## JERC

* `Summer16_07Aug2017*_V11_*_AK4PFchs.txt` <br />
  Files defining JEC and its total uncertainty for 2016. Downloaded from the JERC repository by running `./download_jec.py 2016`. Recommended [here](https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC?rev=166#Jet_Energy_Corrections_in_Run2).
* `Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt` <br />
  Jet p<sub>T</sub> resolution in `Summer16` simulation. File copied from [here](https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Summer16_25nsV1_MC/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt).
* `Summer16_25nsV1_MC_SF_AK4PFchs.txt` <br />
  Data-to-simulation scale factors for jet p<sub>T</sub> resolution for 2016. File copied from [here](https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Summer16_25nsV1_MC/Summer16_25nsV1_MC_SF_AK4PFchs.txt).

* `Fall17_17Nov2017*_V33_*_AK4PFchs.txt` <br />
  Files defining JEC and its total uncertainty for 2017. Downloaded from the JERC repository by running `./download_jec.py 2017`. Recommended [here](https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC?rev=170#Jet_Energy_Corrections_in_Run2).
* `Fall17_V3_MC_PtResolution_AK4PFchs.txt` <br />
  Jet p<sub>T</sub> resolution in `Fall17` simulation. File copied from [here](https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Fall17_V3_MC/Fall17_V3_MC_PtResolution_AK4PFchs.txt).
* `Fall17_V3_MC_SF_AK4PFchs.txt` <br />
  Data-to-simulation scale factors for jet p<sub>T</sub> resolution for 2017. File copied from [here](https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Fall17_V3_MC/Fall17_V3_MC_SF_AK4PFchs.txt).

## LeptonSF

* `2016Muon_SF_ID.root` and `2016Muon_SF_ISO.root`  <br />
  2016 Muon SFs copied from Muon POG [here](https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2016LegacyRereco?rev=10).
* `2017Muon_SF_ID.root` and `2017Muon_SF_ISO.root`  <br />
  2017 Muon SFs copied from Muon POG [here](https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2017?rev=30).
* `2018Muon_SF_ID.root` and `2018Muon_SF_ISO.root`  <br />
  2018 Muon SFs copied from Muon POG [here](https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2018?rev=8).
* `201*Electron_SF_IDISO.root` and `201*Electron_SF_RECO.root`  <br />
  Electron scale factors for all years copied from EGamma POG [here](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations?rev=15#Electron_Scale_Factors).
  There are several sets of electron scale factors. 
  For 2016 RECO SFs, we use `legacy`, `ET > 20GeV`.
  For ID&ISO SFs, we use `Fall17v2` `cutBasedElectronID-Fall17-94X-V2-tight` [here](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations?rev=15#Fall17v2)

## b-tagging efficiencies

* `DeepJet_2016LegacySF_V1.csv`  <br />
  2016 b-tagging SFs copied from BTag POG [here](https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation2016Legacy/DeepJet_2016LegacySF_WP_V1.csv).
* `DeepFlavour_94XSF_V4_B_F.csv`  <br />
  2017 b-tagging SFs copied from BTag POG [here](https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation94X/DeepFlavour_94XSF_WP_V3_B_F.csv).
* `DeepJet_102XSF_V1.csv`  <br />
  2018 b-tagging SFs copied from BTag POG [here](https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation102X/DeepJet_102XSF_WP_V1.csv).

* `btagging_efficiencies_deep2016.root` <br />
  2016 b-tagging efficiencies and mistag-rates computed and stored in 2-d histograms for zz &rarr; 2&ell;2&nu; process.
* `btagging_efficiencies_deep2017.root` <br />
  2017 b-tagging efficiencies and mistag-rates computed and stored in 2-d histograms for zz &rarr; 2&ell;2&nu; process.

