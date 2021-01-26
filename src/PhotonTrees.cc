#include <FileInPath.h>
#include <PhotonTrees.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <TVector2.h>

#include <HZZException.h>
#include <Utils.h>


namespace po = boost::program_options;


int const PhotonTrees::maxSize_;


PhotonTrees::PhotonTrees(Options const &options, Dataset &dataset)
    : EventTrees{options, dataset},
      storeMoreVariables_{options.Exists("more-vars")},
      srcRun_{dataset.Reader(), "run"},
      srcLumi_{dataset.Reader(), "luminosityBlock"},
      srcEvent_{dataset.Reader(), "event"},
      photonBuilder_{dataset},
      photonPrescales_{dataset, options},
      photonWeight_{dataset, options, &photonBuilder_},
      gJetsWeight_{dataset, &photonBuilder_},
      photonFilter_{dataset, options},
      srcNumPVGood_{dataset.Reader(), "PV_npvsGood"} {

  photonBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});
  jetBuilder_.EnableCleaning({&photonBuilder_});
  ptMissBuilder_.PullCalibration({&photonBuilder_});

  weightCollector_.Add(&photonWeight_);
  weightCollector_.Add(&gJetsWeight_);

  CreateWeightBranches();

  AddBranch("jet_cat", &jetCat_);
  AddBranch("jet_size", &jetSize_);
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
    AddBranch("run", &run_);
    AddBranch("lumi", &lumi_);
    AddBranch("event", &event_);

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
  v_jetCat_ = {"_eq0jets","_eq1jets","_geq2jets"};
  applyNvtxWeights_ = Options::NodeAs<bool>(
    options.GetConfig(), {"photon_reweighting", "apply_nvtx_reweighting"});
  applyPtWeights_ = Options::NodeAs<bool>(
    options.GetConfig(), {"photon_reweighting", "apply_pt_reweighting"});
  applyMassLineshape_ = Options::NodeAs<bool>(
    options.GetConfig(), {"photon_reweighting", "apply_mass_lineshape"});
  applyMeanWeights_ = Options::NodeAs<bool>(
    options.GetConfig(), {"photon_reweighting", "apply_mean_weights"});
  utils::loadInstrMETWeights(applyNvtxWeights_, applyPtWeights_, applyMassLineshape_, nVtxWeight_map_, ptWeight_map_, lineshapeMassWeight_map_, v_jetCat_, options);
  utils::loadMeanWeights(applyMeanWeights_, meanWeight_map_, v_jetCat_, options);
}


po::options_description PhotonTrees::OptionsDescription() {
  auto optionsDescription = AnalysisCommon::OptionsDescription();
  optionsDescription.add_options()
    ("more-vars", "Store additional variables");
  return optionsDescription;
}


bool PhotonTrees::ProcessEvent() {
  if (not ApplyCommonFilters())
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

  if (DPhiPtMiss({&jetBuilder_, &photonBuilder_}) < minDphiLeptonsJetsPtMiss_)
    return false;

  if (jets.size() == 0)
    jetCat_ = int(JetCat::kEq0J);
  else if (jets.size() == 1)
    jetCat_ = int(JetCat::kEq1J);
  else
    jetCat_ = int(JetCat::kGEq2J);

  // Only consider photons in the barrel except for Njet >= 2
  if (jets.size() < 2 && !photon->isEB)
    return false;

  // Veto event if the photon fails the following identification
  if (!photon->passElecVeto)
    return false;

  if (photon->sieie < 0.001)
    return false;

  // FIXME temporary. These will be replaced by a new class, much more practical. For now, still use old functions from Utils.
  // Reweighting
  photonReweighting_ = 1.;
  photonNvtxReweighting_ = 1.;
  // In nvtx
  if (applyNvtxWeights_) {
    std::map<double, std::pair<double,double> >::iterator itlow;
    itlow = nVtxWeight_map_["_ll"].upper_bound(*srcNumPVGood_);
    if (itlow == nVtxWeight_map_["_ll"].begin()) 
      throw HZZException{
        "You are trying to access your NVtx reweighting map outside of bin "
        "boundaries."
      };
    itlow--;
    photonReweighting_ *= itlow->second.first;
    photonNvtxReweighting_ *= itlow->second.first;
  }
  // In pT
  if (applyPtWeights_) {
    std::map<double, std::pair<double,double> >::iterator itlow;
    itlow = ptWeight_map_["_ll"+v_jetCat_[jetCat_]].upper_bound(photon->p4.Pt()); //look at which bin in the map currentEvt.pT corresponds
    if (itlow == ptWeight_map_["_ll" + v_jetCat_[jetCat_]].begin())
      throw HZZException{
        "You are trying to access your Pt reweighting map outside of bin "
        "boundaries."
      };
    itlow--;
    photonReweighting_ *= itlow->second.first; //don't apply for first element of the map which is the normal one without reweighting.
  }

  // Give mass to photon
  double photonMass = 0;
  if (applyMassLineshape_) {
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

  // Apply the event veto for photons failing an OR of the addtional IDs of photon PF ID, sigmaIPhiIPhi > 0.001,
  //  MIPTotalEnergy < 4.9, |seedtime| < 2ns, and seedtime < 1ns for 2018
  if (!photonFilter_())
    return false;

  // FIXME temporary. These will be replaced by a new class, much more practical. For now, still use old functions from Utils.
  // Get mean weights
  if (applyMeanWeights_) {
    std::map<double, double>::iterator itlow;
    itlow = meanWeight_map_[v_jetCat_[jetCat_]].upper_bound(mT_); //look at which bin in the map mt corresponds
    if (itlow == meanWeight_map_[v_jetCat_[jetCat_]].begin())
      throw HZZException{
        "You are trying to access your mean weight map outside of bin "
        "boundaries."
      };
    itlow--;
    meanWeight_ = itlow->second;
  }

  if (applyMeanWeights_ and meanWeight_ == 0)
    return false;

  numPVGood_ = *srcNumPVGood_;


  if (storeMoreVariables_)
    FillMoreVariables(jets);

  FillTree();
  return true;

}


Photon const *PhotonTrees::CheckPhotons() const {
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
  run_ = *srcRun_;
  lumi_ = *srcLumi_;
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
