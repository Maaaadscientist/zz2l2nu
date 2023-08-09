#include <fstream>
#include <FileInPath.h>
#include <ElectronTrees.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <TVector2.h>

#include <HZZException.h>
#include <Utils.h>



namespace po = boost::program_options;


int const ElectronTrees::maxSize_;


ElectronTrees::ElectronTrees(Options const &options, Dataset &dataset)
    : EventTrees{options, dataset},
      storeMoreVariables_{options.Exists("more-vars")},
      srcRun_{dataset.Reader(), "run"},
      srcLumi_{dataset.Reader(), "luminosityBlock"},
      srcEvent_{dataset.Reader(), "event"},
      photonBuilder_{dataset},
      // photonWeight_{dataset, options, &photonBuilder_},
      // photonFilter_{dataset, options},
      srcNumPVGood_{dataset.Reader(), "PV_npvsGood"} {

  if (isSim_) {
    numGenPart_.reset(new TTreeReaderValue<UInt_t>(dataset.Reader(), "nGenPart"));
    genPartPdgId_.reset(new TTreeReaderArray<Int_t>(dataset.Reader(), "GenPart_pdgId"));
    genPartPt_.reset(new TTreeReaderArray<Float_t>(dataset.Reader(), "GenPart_pt"));
    genPartStatus_.reset(new TTreeReaderArray<Int_t>(dataset.Reader(), "GenPart_status"));
    genPartStatusFlags_.reset(new TTreeReaderArray<Int_t>(dataset.Reader(), "GenPart_statusFlags"));
  }

  photonBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});

  // weightCollector_.Add(&photonWeight_);

  CreateWeightBranches();

  AddBranch("jet_cat", &jetCat_);
  AddBranch("jet_size", &jetSize_);
  AddBranch("electron_pt", &electronPt_);
  AddBranch("electron_eta", &electronEta_);
  AddBranch("electron_phi", &electronPhi_);
  AddBranch("electron_M", &electronM_);
  AddBranch("electron_eta_sc", &electronEtaSc_);
  // AddBranch("dijet_M", &dijetM_);
  AddBranch("ptmiss", &missPt_);
  AddBranch("ptmiss_phi", &missPhi_);
  AddBranch("electron_MET_deltaPhi", &electronMetDeltaPhi_);
  AddBranch("electron_MET_Mt", &electronMetMt_);
  AddBranch("num_pv_good", &numPVGood_);


  if (storeMoreVariables_) {
    AddBranch("run", &run_);
    AddBranch("lumi", &lumi_);
    AddBranch("event", &event_);

    AddBranch("jet_pt", jetPt_, "jet_pt[jet_size]/F");
    AddBranch("jet_eta", jetEta_, "jet_eta[jet_size]/F");
    AddBranch("jet_phi", jetPhi_, "jet_phi[jet_size]/F");
    AddBranch("jet_mass", jetMass_, "jet_mass[jet_size]/F");
  }

  auto const &isQCDNode = dataset.Info().Parameters()["mc_qcd"];
  isQCD_ = (isQCDNode and not isQCDNode.IsNull() and isQCDNode.as<bool>());

}


po::options_description ElectronTrees::OptionsDescription() {
  auto optionsDescription = AnalysisCommon::OptionsDescription();
  optionsDescription.add_options()
    ("more-vars", "Store additional variables");
  return optionsDescription;
}


bool ElectronTrees::ProcessEvent() {

  // Resolve G+jet/QCD mixing (avoid double counting of photons):
  // QCD samples allow prompt photons of pT > 10, for gamma+jets it's 25
  if (isSim_ && isQCD_) {
    for (unsigned i = 0; i < genPartPdgId_->GetSize(); ++i) {
      // Particle is in the final state
      if (not (genPartStatus_->At(i) == 1))
        continue;

      // is a photon
      if (not (genPartPdgId_->At(i) == 22))
        continue;

      // isPrompt
      if (not ((genPartStatusFlags_->At(i) & 1)))
        continue;

      // pT > 25 GeV
      if (not (genPartPt_->At(i) > 25.))
        continue;

      // std::cout << "Prompt photon with pT > 25 GeV found!" << std::endl;
      return false;
    }
  }

  if (not ApplyCommonFilters()) {
    //if(sel) std::cout << "not pass common filters" <<std::endl;
    return false;
  }

  auto const electron = CheckElectron();
  if (electron == nullptr) {
    return false;
  }

  if (electron->p4.Pt() < minPtLL_) {
    return false;
  }

  auto const &p4Miss = ptMissBuilder_.Get().p4;
  missPt_ = p4Miss.Pt();
  missPhi_ = p4Miss.Phi();

  if (std::abs(TVector2::Phi_mpi_pi(electron->p4.Phi() - p4Miss.Phi()))
        < minDphiLLPtMiss_)
  {
    //if(sel) std::cout<<"photon, met dphi = "<<std::abs(TVector2::Phi_mpi_pi(photon->p4.Phi() - p4Miss.Phi())) <<" < min dphi"<<std::endl;
    return false;
  }

  auto const &jets = jetBuilder_.Get();
  jetSize_ = jets.size();

  if (jets.size() == 0)
    jetCat_ = int(JetCat::kEq0J);
  else if (jets.size() == 1)
    jetCat_ = int(JetCat::kEq1J);
  else
    jetCat_ = int(JetCat::kGEq2J);

  // Only consider electrons in the barrel except for Njet >= 2
  // if (jets.size() < 2 && !(std::abs(electron->etaSc) < 1.4442)) {
  //   return false;
  // }

  // if (jets.size() < 2) {
  //   return false;
  // }

  for (auto const &jet : jets) {
    if (bTagger_(jet)) {
      //if (sel) std::cout<<"b-tag jet found"<<std::endl;
      return false;
    }

    if (std::abs(TVector2::Phi_mpi_pi(jet.p4.Phi() - p4Miss.Phi())) 
          < minDphiJetsPtMiss_) {
      return false;
      }
  }

  if (DPhiPtMiss({&jetBuilder_, &electronBuilder_}) < minDphiLeptonsJetsPtMiss_) {
    return false;

  }

  // if (jets[0].p4.Eta() * jets[1].p4.Eta() >= 0)
  //   return false;

  // auto dijetP4 = jets[0].p4 + jets[1].p4;
  // dijetM_ = dijetP4.M();
  // if (dijetM_ <= 400.0)
  //   return false;

  // if (std::abs(jets[0].p4.Eta() - jets[1].p4.Eta()) <= 2.4)
  //   return false;

  electronPt_ = electron->p4.Pt();
  electronEta_ = electron->p4.Eta();
  electronPhi_ = electron->p4.Phi();
  electronM_ = electron->p4.M();

  electronEtaSc_ = electron->etaSc;

  numPVGood_ = *srcNumPVGood_;

  electronMetDeltaPhi_ = TVector2::Phi_mpi_pi(electron->p4.Phi() - p4Miss.Phi());
  electronMetMt_ = std::sqrt(2 * electron->p4.Pt() * p4Miss.Pt() * (1 - std::cos(electronMetDeltaPhi_)));

  if (storeMoreVariables_)
    FillMoreVariables(jets);

  FillTree();
  // std::cout << eventInfo <<std::endl;
  return true;

}

Electron const *ElectronTrees::CheckElectron() const {
  auto const &looseElectrons = electronBuilder_.GetLoose();
  auto const &tightElectrons = electronBuilder_.GetLoose();

  auto const &looseMuons = muonBuilder_.GetLoose();
  auto const &photons = photonBuilder_.Get();

  if (looseElectrons.size() != tightElectrons.size()) {
    return nullptr;
  }
  if (looseMuons.size() > 0) {
    return nullptr;
  }
  if (photons.size() > 0) {
    return nullptr;
  }
  if (isotrkBuilder_.Get().size() > 0) {
    return nullptr;
  }

  if (tightElectrons.size() != 1) {
    return nullptr;
  }

  return &tightElectrons[0];
}

void ElectronTrees::FillMoreVariables(std::vector<Jet> const &jets) {
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
