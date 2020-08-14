#include <FileInPath.h>
#include <PhotonTrees.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include <TVector2.h>

#include <Utils.h>


namespace po = boost::program_options;


int const PhotonTrees::maxSize_;


PhotonTrees::PhotonTrees(Options const &options, Dataset &dataset)
    : EventTrees{options, dataset},
      storeMoreVariables_{options.Exists("more-vars")},
      srcEvent_{dataset.Reader(), "event"},
      photonPrescales_{dataset, options},
      photonWeight_{dataset, options, &photonBuilder_},
      gJetsWeight_{dataset, &photonBuilder_},
      srcNumPVGood_{dataset.Reader(), "PV_npvsGood"} {

  photonBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});
  jetBuilder_.EnableCleaning({&photonBuilder_});
  ptMissBuilder_.PullCalibration({&photonBuilder_});

  weightCollector_.Add(&photonWeight_);
  weightCollector_.Add(&gJetsWeight_);

  CreateWeightBranches();

  AddBranch("jet_cat", &jetCat_);
  AddBranch("photon_pt", &photonPt_);
  AddBranch("photon_eta", &photonEta_);
  AddBranch("photon_phi", &photonPhi_);
  AddBranch("photon_mass", &photonMass_);
  AddBranch("ptmiss", &missPt_);
  AddBranch("ptmiss_phi", &missPhi_);
  AddBranch("mT", &mT_);
  AddBranch("num_pv_good", &numPVGood_);
  AddBranch("trigger_weight", &triggerWeight_);
  AddBranch("photon_reweighting", &photonReweighting_);
  AddBranch("photon_nvtx_reweighting", &photonNvtxReweighting_);
  AddBranch("mean_weight", &meanWeight_);

  if (storeMoreVariables_) {
    AddBranch("event", &event_);

    AddBranch("jet_size", &jetSize_);
    AddBranch("jet_pt", jetPt_, "jet_pt[jet_size]/F");
    AddBranch("jet_eta", jetEta_, "jet_eta[jet_size]/F");
    AddBranch("jet_phi", jetPhi_, "jet_phi[jet_size]/F");
    AddBranch("jet_mass", jetMass_, "jet_mass[jet_size]/F");
  }

  auto const WGSettingsNode = dataset.Info().Parameters()["wgamma_lnugamma"];
  if (WGSettingsNode and not WGSettingsNode.IsNull()) {
    labelWGamma_ = WGSettingsNode.as<std::string>();
    genPhotonBuilder_.emplace(dataset);
  }
  auto const ZGSettingsNode = dataset.Info().Parameters()["zgamma_2nugamma"];
  if (ZGSettingsNode and not ZGSettingsNode.IsNull()) {
    labelZGamma_ = ZGSettingsNode.as<std::string>();
    genPhotonBuilder_.emplace(dataset);
  }

  auto const &isQCDNode = dataset.Info().Parameters()["mc_qcd"];
  isQCD_ = (isQCDNode and not isQCDNode.IsNull() and isQCDNode.as<bool>());

  // FIXME temporary. These will be replaced by a new class, much more practical. For now, still use old functions from Utils.
  v_jetCat_ = {"_eq0jets","_geq1jets","_vbf"};
  auto const base_path = std::string(std::getenv("HZZ2L2NU_BASE")) + "/";
  nvtxReweightingFile_ = FileInPath::Resolve(Options::NodeAs<std::string>(
    options.GetConfig(), {"photon_reweighting", "nvtx_reweighting"}));
  ptReweightingFile_ = FileInPath::Resolve(Options::NodeAs<std::string>(
    options.GetConfig(), {"photon_reweighting", "pt_reweighting"}));
  massLineshapeFile_ = FileInPath::Resolve(Options::NodeAs<std::string>(
    options.GetConfig(), {"photon_reweighting", "mass_lineshape"}));
  meanWeightsFile_ = FileInPath::Resolve(Options::NodeAs<std::string>(
    options.GetConfig(), {"photon_reweighting", "mean_weights"}));
  weight_NVtx_exist_ = utils::file_exist(nvtxReweightingFile_);
  weight_Pt_exist_ = utils::file_exist(ptReweightingFile_);
  weight_Mass_exist_ = utils::file_exist(massLineshapeFile_);
  meanWeight_exist_ = utils::file_exist(meanWeightsFile_);
  utils::loadInstrMETWeights(weight_NVtx_exist_, weight_Pt_exist_, weight_Mass_exist_, nvtxReweightingFile_, ptReweightingFile_, massLineshapeFile_, nVtxWeight_map_, ptWeight_map_, lineshapeMassWeight_map_, v_jetCat_);
  utils::loadMeanWeights(meanWeight_exist_, meanWeightsFile_, meanWeight_map_, v_jetCat_);
}


po::options_description PhotonTrees::OptionsDescription() {
  auto optionsDescription = AnalysisCommon::OptionsDescription();
  optionsDescription.add_options()
    ("more-vars", "Store additional variables");
  return optionsDescription;
}


bool PhotonTrees::ProcessEvent() {
  if (not meKinFilter_() or not metFilters_())
    return false;

  auto const photon = CheckPhotons();
  if (photon == nullptr)
    return false;

  // Avoid double counting for ZGamma overlap between 2 samples
  if (labelZGamma_ != "" and genPhotonBuilder_) {
    double genPhotonPt = genPhotonBuilder_->P4Gamma().Pt();
    if (labelZGamma_ == "inclusive" and genPhotonPt >= 130)
      return false;
    if (labelZGamma_ == "pt130" and genPhotonPt < 130)
      return false;
  }

  // Avoid double counting for WGamma overlap between 2 samples
  if (labelWGamma_ != "" and genPhotonBuilder_) {
    double genPhotonPt = genPhotonBuilder_->P4Gamma().Pt();
    if (labelWGamma_ == "inclusive" and genPhotonPt >= 130)
      return false;
    if (labelWGamma_ == "pt130" and (genPhotonPt < 130 or genPhotonPt >= 500))
      return false;
    if (labelWGamma_ == "pt500" and genPhotonPt < 500)
      return false;
  }

  // Resolve G+jet/QCD mixing (avoid double counting of photons):
  // QCD samples allow prompt photons of pT > 10, for gamma+jets it's 25
  if (isQCD_) {
    auto const &photons = photonBuilder_.Get();
    for (int i = 0 ; i < int(photons.size()) ; i++) {
      if (photons[i].flavour == Photon::Origin::PromptPhoton 
            and photons[i].genP4.Pt() > 25.) {
        return false;  // Remove all QCD events with a photon of pT > 25
      }
    }
  }

  if (photon->p4.Pt() < minPtLL_)
    return false;


  auto const &p4Miss = ptMissBuilder_.Get().p4;
  missPt_ = p4Miss.Pt();
  missPhi_ = p4Miss.Phi();

  if (std::abs(TVector2::Phi_mpi_pi(photon->p4.Phi() - p4Miss.Phi())) 
        < minDphiLLPtMiss_)
    return false;


  auto const &jets = jetBuilder_.Get();

  for (auto const &jet : jets) {
    if (bTagger_(jet))
      return false;

    if (std::abs(TVector2::Phi_mpi_pi(jet.p4.Phi() - p4Miss.Phi())) 
          < minDphiJetsPtMiss_)
      return false;
  }

  if (DPhiLeptonsJetsSystemPtMiss(false) < minDphiLeptonsJetsPtMiss_)
    return false;

  if (jets.size() == 0)
    jetCat_ = int(JetCat::kEq0J);
  else if (utils::PassVbfCuts(jets, photon->p4))
    jetCat_ = int(JetCat::kVbf);
  else
    jetCat_ = int(JetCat::kGEq1J);

  // FIXME temporary. These will be replaced by a new class, much more practical. For now, still use old functions from Utils.
  // Reweighting
  photonReweighting_ = 1.;
  photonNvtxReweighting_ = 1.;
  // In nvtx
  if (weight_NVtx_exist_) {
    std::map<std::pair<double,double>, std::pair<double,double> >::iterator itlow;
    itlow = nVtxWeight_map_["_ll"].upper_bound(std::make_pair(*srcNumPVGood_,photon->p4.Pt())); //look at which bin in the map currentEvt.rho corresponds
    if (itlow == nVtxWeight_map_["_ll"].begin()) 
      throw std::out_of_range("You are trying to access your NVtx reweighting map outside of bin boundaries");
    itlow--;
    photonReweighting_ *= itlow->second.first;
    photonNvtxReweighting_ *= itlow->second.first;
  }
  // In pT
  if (weight_Pt_exist_) {
    std::map<double, std::pair<double,double> >::iterator itlow;
    itlow = ptWeight_map_["_ll"+v_jetCat_[jetCat_]].upper_bound(photon->p4.Pt()); //look at which bin in the map currentEvt.pT corresponds
    if (itlow == ptWeight_map_["_ll"+v_jetCat_[jetCat_]].begin()) 
      throw std::out_of_range("You are trying to access your Pt reweighting map outside of bin boundaries");
    itlow--;
    photonReweighting_ *= itlow->second.first; //don't apply for first element of the map which is the normal one without reweighting.
  }

  // Give mass to photon
  double photonMass = 0;
  if (weight_Mass_exist_) {
    photonMass = lineshapeMassWeight_map_["_ll"]->GetRandom();
  }
  TLorentzVector photonWithMass;
  photonWithMass.SetPtEtaPhiM(photon->p4.Pt(), photon->p4.Eta(), 
    photon->p4.Phi(), photonMass);

  photonPt_ = photonWithMass.Pt();
  photonEta_ = photonWithMass.Eta();
  photonPhi_ = photonWithMass.Phi();
  photonMass_ = photonWithMass.M();


  double const eT =
      std::sqrt(std::pow(photonWithMass.Pt(), 2)
          + std::pow(photonWithMass.M(), 2))
      + std::sqrt(std::pow(p4Miss.Pt(), 2) + std::pow(kNominalMZ_, 2));
  mT_ = std::sqrt(
      std::pow(eT, 2) - std::pow((photonWithMass + p4Miss).Pt(), 2));

  triggerWeight_ = photonPrescales_.GetPhotonPrescale(photonWithMass.Pt());
  if (triggerWeight_ == 0)
    return false;

  // FIXME temporary. These will be replaced by a new class, much more practical. For now, still use old functions from Utils.
  // Get mean weights
  if (meanWeight_exist_) {
    std::map<double, double>::iterator itlow;
    itlow = meanWeight_map_[v_jetCat_[jetCat_]].upper_bound(mT_); //look at which bin in the map mt corresponds
    if (itlow == meanWeight_map_[v_jetCat_[jetCat_]].begin()) 
      throw std::out_of_range("You are trying to access your mean weight map outside of bin boundaries");
    itlow--;
    meanWeight_ = itlow->second;
  }

  if (meanWeight_ == 0)
    return false;

  numPVGood_ = *srcNumPVGood_;


  if (storeMoreVariables_)
    FillMoreVariables(jets);

  FillTree();
  return true;

}

Photon const * PhotonTrees::CheckPhotons() const {
  auto const &looseElectrons = electronBuilder_.GetLoose();
  auto const &looseMuons = muonBuilder_.GetLoose();
  auto const &photons = photonBuilder_.Get();

  if (looseElectrons.size() + looseMuons.size() > 0)
    return nullptr;
  if (photons.size() != 1)
    return nullptr;
  if (isotrkBuilder_.Get().size() > 0)
    return nullptr;

  return &photons[0];
}

void PhotonTrees::FillMoreVariables(std::vector<Jet> const &jets) {

  event_ = *srcEvent_;

  jetSize_ = std::min<int>(jets.size(), maxSize_);

  for (int i = 0; i < jetSize_; ++i) {
    auto const &p4 = jets[i].p4;
    jetPt_[i] = p4.Pt();
    jetEta_[i] = p4.Eta();
    jetPhi_[i] = p4.Phi();
    jetMass_[i] = p4.M();
  }

}

// Still missing:
// - Computation of weights, uncertainties on the weights etc...
