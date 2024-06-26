# Data files

This directory contains auxiliary files with various corrections. They are organized into several subdirectories:


## BTag

* `DeepJet_2016LegacySF_V1.csv`  <br />
  2016 b-tagging SFs copied from BTag POG [here](https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation2016Legacy/DeepJet_2016LegacySF_WP_V1.csv).
* `DeepFlavour_94XSF_V4_B_F.csv`  <br />
  2017 b-tagging SFs copied from BTag POG [here](https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation94X/DeepFlavour_94XSF_WP_V3_B_F.csv).
* `DeepJet_102XSF_V1.csv`  <br />
  2018 b-tagging SFs copied from BTag POG [here](https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation102X/DeepJet_102XSF_WP_V1.csv).

* `btagging_efficiencies_deep2016.root` <br />
  2016 b-tagging efficiencies and mistag-rates computed and stored in 2-d histograms for ZZ &rarr; 2&ell;2&nu; process.
* `btagging_efficiencies_deep2017.root` <br />
  2017 b-tagging efficiencies and mistag-rates computed and stored in 2-d histograms for ZZ &rarr; 2&ell;2&nu; process.
* `btagging_efficiencies_deep2018.root` <br />
  2018 b-tagging efficiencies and mistag-rates computed and stored in 2-d histograms for ZZ &rarr; 2&ell;2&nu; process.


## InstrMetReweighting

* `lineshape_mass_201*.root` <br />
  Mass lineshapes computed from dilepton data with pTmiss < 125 GeV. Used in the CR to give a mass to the photon. Obtained from `compute_mass_lineshape.py`.
* `weight_nvtx_201*.root` <br />
  2D weights in number of vertices vs boson pT (in bins of pT threshold), computed from dilepton and photon data with pTmiss < 125 GeV. Used to reweight the data in the photon CR. Obtained from `compute_instrMET_weights.py`.
* `weight_pt_201*.root` <br />
  Weights in number of vertices vs boson pT (in bins of pT threshold), computed from dilepton and photon data with pTmiss < 125 GeV on top of nvtx weights. Used to reweight the data in the photon CR. Obtained from `compute_instrMET_weights.py`.
* `meanWeights_201*.root` <br />
  Mean weights per bin in mT, to be used with the statistical analysis. Computed from the code in `WeightsAndDatadriven/InstrMET/meanWeights`.


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

* `Autumn18*_V19_*_AK4PFchs.txt` <br />
  Files defining JEC and its total uncertainty for 2018. Downloaded from the JERC repository by running `./download_jec.py 2018`. Recommended [here](https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC?rev=174#Jet_Energy_Corrections_in_Run2).
* `Autumn18_V7b_MC_PtResolution_AK4PFchs.txt` <br />
  Jet p<sub>T</sub> resolution in `Autumn18` simulation. File copied from [here](https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Autumn18_V7b_MC/Autumn18_V7b_MC_PtResolution_AK4PFchs.txt).
* `Autumn18_V7b_MC_SF_AK4PFchs.txt` <br />
  Data-to-simulation scale factors for jet p<sub>T</sub> resolution for 2018. File copied from [here](https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Autumn18_V7b_MC/Autumn18_V7b_MC_SF_AK4PFchs.txt).


## LeptonSF

* `201X/Muon_efficiencies_id_iso_201X.root`
  SFs for both muon ID and Iso for customized lepton WPs for year 201X, approval at [here](https://indico.cern.ch/event/943782/contributions/3968378/)
  For ID, we use `Muon_mediumPromptId`, for Iso, we use `Muon_pfRelIso03_all < 0.15`

* `201X/Electron_efficiencies_id_iso_201X.root`
  SFs for both ID and Iso for customized lepton WPs for year 201X, approval at [here](https://indico.cern.ch/event/879930/#4-approval-of-h-zz-id-and-hlt)
  For ID, we use `Electron_mvaFall17V2noIso_WP90`, for Iso, we use `Electron_pfRelIso03_all < 0.1`

* `201X/Electron_efficiencies_tracking_201X.root`
  Electron reconstruction scale factors for all years copied from EGamma POG [here](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale).
  For 2016 RECO SFs, we use `legacy`, `ET > 20GeV`.


## Lumi

* `lumi.yaml` <br />
  Per-run recorded integrated luminosities for Run 2. Computed with ‘golden’ certification files and ‘physics’ calibration.


## PhotonSF

* `2016_Photon_SF_ID.root`  <br />
  2016 Photon ID SF (tight WP) copied from EGamma POG [here](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2?rev=116#94X_series_Fall17V2_Scale_fa_AN2).
  These SF are for ID `Fall17v2` based on EGamma POG recommendation ID for 2016 [here](https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2?rev=52#Cut_Based_Photon_ID_for_Run_2).

* `2017_Photon_SF_ID.root`  <br />
  2016 Photon ID SF (tight WP) copied from EGamma POG [here](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2?rev=116#94X_series_Fall17V2_IDs_Scale_fa).
  These SF are for ID `Fall17v2` based on EGamma POG recommendation ID for 2017 [here](https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2?rev=52#Cut_Based_Photon_ID_for_Run_2).

* `2017_Photon_SF_ID.root`  <br />
  2016 Photon ID SF (tight WP) copied from EGamma POG [here](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2?rev=116#102X_series_Fall17V2_IDs_Sca_AN1).
  These SF are for ID `Fall17v2` based on EGamma POG recommendation ID for 2018 [here](https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2?rev=52#Cut_Based_Photon_ID_for_Run_2).


## PileupID

* `scalefactorsPUID_81Xtraining.root` <br />
  Scale factors for pileup ID with 81X training. Copied from the [pileup ID wiki page](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID?rev=61).
* `pileup_eff_nanoaodv6.xgb` <br />
  XGBoost model that parameterizes the efficiency of the pileup ID in simulation. It was trained with NanoAODv6. Details are provided in [this talk](https://indico.cern.ch/event/934117/#9-parametrization-of-puid-effc).


## VBF Discriminants' g and c constants

* `gConstant_VBF_g2.root` <br />
   `g` constants of `a2` anomalous coupling stored in `TSpline`. They are used in `DjjVBF` discriminant calculation for the model assuming `a2` anomalous coupling in HVV scattering amplitude. Copied form [here](https://github.com/usarica/CMS3NtupleTools/blob/aaa69bc1d2d9794c4653bd6a5b3f038ead76fae3/AnalysisTree/data/RecoMEConstants/gConstant_VBF_g2.root).
* `gConstant_VBF_g4.root` <br />
   `g` constants `a3` anomalous coupling stored in `TSpline`. They are used in `DjjVBF` discriminant calculation for the model assuming `a3` anomalous coupling in HVV scattering amplitude. Copied from [here](https://github.com/usarica/CMS3NtupleTools/blob/aaa69bc1d2d9794c4653bd6a5b3f038ead76fae3/AnalysisTree/data/RecoMEConstants/gConstant_VBF_g4.root).
* `gConstant_VBF_L1.root` <br />
   `g` constants `L1` anomalous coupling stored in `TSpline`. They are used in `DjjVBF` discriminant calculation for the model assuming `L1` anomalous coupling in HVV scattering amplitude. Copied from [here](https://github.com/usarica/CMS3NtupleTools/blob/aaa69bc1d2d9794c4653bd6a5b3f038ead76fae3/AnalysisTree/data/RecoMEConstants/gConstant_VBF_L1.root).
* `SmoothKDConstant_m4l_DjjVBF_13TeV.root` <br />
   `c` constants stored in TSpline. They are used in DjjVBF discriminants calculation for all target models. Copied from [here](https://github.com/usarica/CMS3NtupleTools/blob/aaa69bc1d2d9794c4653bd6a5b3f038ead76fae3/AnalysisTree/data/RecoMEConstants/SmoothKDConstant_m4l_DjjVBF_13TeV.root).
