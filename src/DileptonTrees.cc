#include <DileptonTrees.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include <TVector2.h>

#include <LeptonsEfficiencySF.h>
#include <Utils.h>


namespace po = boost::program_options;


DileptonTrees::DileptonTrees(Options const &options, Dataset &dataset)
    : dataset_{dataset},
      intLumi_{options.GetConfig()["luminosity"].as<double>()},
      storeMoreVariables_{options.Exists("more-vars")},
      randomGenerator_{options.GetAs<unsigned>("seed")},
      bTagger_{options},
      electronBuilder_{dataset_, options},
      muonBuilder_{dataset_, options, randomGenerator_},
      jetBuilder_{dataset_, options, randomGenerator_},
      ptMissBuilder_{dataset_},
      meKinFilter_{dataset_}, metFilters_{dataset_},
      bTagWeight_{options, bTagger_},
      outputFile_{options.GetAs<std::string>("output").c_str(), "recreate"},
      p4LL_{nullptr}, p4Miss_{nullptr},
      weight_{1.} {
        
  bool const isSim = dataset_.Info().IsSimulation();

  if (isSim) {
    auto const &node = dataset_.Info().Parameters()["zz_2l2nu"];

    if (node and not node.IsNull() and node.as<bool>())
      genZZBuilder_.reset(new GenZZBuilder(dataset));
  }

  if (isSim) {
    genJetBuilder_.reset(new GenJetBuilder(dataset_, options));
    jetBuilder_.SetGenJetBuilder(genJetBuilder_.get());
  }
  
  jetBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});
  
  // Type 1 corrections for ptmiss
  ptMissBuilder_.PullCalibration(
      {&muonBuilder_, &electronBuilder_, &jetBuilder_});

  if (isSim) {
    genWeight_.reset(new GenWeight{dataset_});
    ewCorrectionWeight_.reset(new EWCorrectionWeight{dataset_, options});
    pileUpWeight_.reset(new PileUpWeight{dataset_, options});
    kFactorCorrection_.reset(new KFactorCorrection{dataset_, options});
  }


  tree_ = new TTree("Vars", "Variables in dilepton selection");
  tree_->SetDirectory(&outputFile_);
  tree_->Branch("leptonCat", &leptonCat_);
  tree_->Branch("jetCat", &jetCat_);
  tree_->Branch("p4LL", &p4LL_);
  tree_->Branch("p4Miss", &p4Miss_);
  tree_->Branch("mT", &mT_);

  if (storeMoreVariables_) {
    if (genZZBuilder_)
      tree_->Branch("genMZZ", &genMZZ_);

    tree_->Branch("lepton_charge", leptonCharge_, "lepton_charge[2]/I");
    tree_->Branch("lepton_pt", leptonPt_, "lepton_pt[2]/F");
    tree_->Branch("lepton_eta", leptonEta_, "lepton_eta_[2]/F");
    tree_->Branch("lepton_phi", leptonPhi_, "lepton_phi[2]/F");
    tree_->Branch("lepton_mass", leptonMass_, "lepton_mass[2]/F");

    tree_->Branch("jet_size", &jetSize_);
    tree_->Branch("jet_pt", jetPt_, "jet_pt[jet_size]/F");
    tree_->Branch("jet_eta", jetEta_, "jet_eta[jet_size]/F");
    tree_->Branch("jet_phi", jetPhi_, "jet_phi[jet_size]/F");
    tree_->Branch("jet_mass", jetMass_, "jet_mass[jet_size]/F");
  }

  tree_->Branch("weight", &weight_);
}


po::options_description DileptonTrees::OptionsDescription() {
  po::options_description optionsDescription{"Analysis-specific options"};
  optionsDescription.add_options()
    ("syst", po::value<std::string>()->default_value(""),
     "Requested systematic variation")
    ("output,o", po::value<std::string>()->default_value("output.root"),
     "Name for output file with histograms")
    ("seed", po::value<unsigned>()->default_value(0),
     "Seed for random number generator; 0 means a unique seed")
    ("more-vars", "Store additional variables");
  return optionsDescription;
}


void DileptonTrees::PostProcessing() {
  outputFile_.Write();
  outputFile_.Close();
}


bool DileptonTrees::ProcessEvent() {
  if (not meKinFilter_() or not metFilters_())
    return false;

  auto const leptonResult = CheckLeptons();
  if (not leptonResult)
    return false;
  auto const &[leptonCat, l1, l2] = leptonResult.value();

  leptonCat_ = int(leptonCat);
  *p4LL_ = l1->p4 + l2->p4;

  if (std::abs(p4LL_->M() - kNominalMZ_) > 15.)
    return false;

  if (p4LL_->Pt() < 55.)
    return false;


  *p4Miss_ = ptMissBuilder_.Get().p4;

  if (p4Miss_->Pt() < 80.)
    return false;

  if (std::abs(TVector2::Phi_mpi_pi(p4LL_->Phi() - p4Miss_->Phi())) < 0.5)
    return false;


  auto const &jets = jetBuilder_.Get();

  for (auto const &jet : jets) {
    if (bTagger_(jet))
      return false;

    if (std::abs(TVector2::Phi_mpi_pi(jet.p4.Phi() - p4Miss_->Phi())) < 0.5)
      return false;
  }

  if (jets.size() == 0)
    jetCat_ = int(JetCat::kEq0J);
  else if (utils::PassVbfCuts(jets, *p4LL_))
    jetCat_ = int(JetCat::kVbf);
  else
    jetCat_ = int(JetCat::kGEq1J);


  double const eT =
      std::sqrt(std::pow(p4LL_->Pt(), 2) + std::pow(p4LL_->M(), 2))
      + std::sqrt(std::pow(p4Miss_->Pt(), 2) + std::pow(kNominalMZ_, 2));
  mT_ = std::sqrt(std::pow(eT, 2) - std::pow((*p4LL_ + *p4Miss_).Pt(), 2));
  

  if (storeMoreVariables_)
    FillMoreVariables({*l1, *l2}, jets);

  if (dataset_.Info().IsSimulation())
    weight_ = SimWeight(leptonCat);

  tree_->Fill();
  return true;
}


std::optional<std::tuple<DileptonTrees::LeptonCat, Lepton const *, Lepton const *>>
DileptonTrees::CheckLeptons() const {
  auto const &tightElectrons = electronBuilder_.GetTight();
  auto const &looseElectrons = electronBuilder_.GetLoose();

  auto const &tightMuons = muonBuilder_.GetTight();
  auto const &looseMuons = muonBuilder_.GetLoose();

  if (looseElectrons.size() + looseMuons.size() != 2)
    return {};

  LeptonCat leptonCat;
  Lepton const *l1, *l2;

  if (tightElectrons.size() == 2) {
    leptonCat = LeptonCat::kEE;
    l1 = &tightElectrons[0];
    l2 = &tightElectrons[1];
  } else if (tightMuons.size() == 2) {
    leptonCat = LeptonCat::kMuMu;
    l1 = &tightMuons[0];
    l2 = &tightMuons[1];
  } else if (tightElectrons.size() == 1 and tightMuons.size() == 1) {
    leptonCat = LeptonCat::kEMu;
    l1 = &tightElectrons[0];
    l2 = &tightMuons[0];
    if (l1->p4.Pt() < l2->p4.Pt())
      std::swap(l1, l2);
  } else
    return {};

  return std::make_tuple(leptonCat, l1, l2);
}


void DileptonTrees::FillMoreVariables(
    std::array<Lepton, 2> const &leptons, std::vector<Jet> const &jets) {
  
  if (genZZBuilder_)
    genMZZ_ = genZZBuilder_->P4ZZ().M();

  for (int i = 0; i < 2; ++i) {
    leptonCharge_[i] = leptons[i].charge;
    auto const &p4 = leptons[i].p4;
    leptonPt_[i] = p4.Pt();
    leptonEta_[i] = p4.Eta();
    leptonPhi_[i] = p4.Phi();
    leptonMass_[i] = p4.M();
  }

  jetSize_ = std::min<int>(jets.size(), maxSize_);

  for (int i = 0; i < jetSize_; ++i) {
    auto const &p4 = jets[i].p4;
    jetPt_[i] = p4.Pt();
    jetEta_[i] = p4.Eta();
    jetPhi_[i] = p4.Phi();
    jetMass_[i] = p4.M();
  }
}


double DileptonTrees::LeptonSFWeight(LeptonCat leptonCat) const {
  auto const &electrons = electronBuilder_.GetTight();
  auto const &muons = muonBuilder_.GetTight();

  switch (leptonCat) {
    case LeptonCat::kEE:
      return trigAndIDsfs::diElectronEventSFs(
        utils::CutVersion::CutSet::Moriond17Cut,
        electrons[0].p4.Pt(), electrons[0].etaSc,
        electrons[1].p4.Pt(), electrons[1].etaSc);

    case LeptonCat::kMuMu:
      return trigAndIDsfs::diMuonEventSFs(
        utils::CutVersion::CutSet::Moriond17Cut,
        muons[0].uncorrP4.Pt(), muons[0].uncorrP4.Eta(),
        muons[1].uncorrP4.Pt(), muons[1].uncorrP4.Eta());

    case LeptonCat::kEMu:
      return trigAndIDsfs::EMuEventSFs(
        utils::CutVersion::CutSet::Moriond17Cut,
        muons[0].uncorrP4.Pt(), muons[0].uncorrP4.Eta(),
        electrons[0].p4.Pt(), electrons[0].etaSc);
  }

  // This should never happen
  throw std::logic_error("Unhandled lepton category.");
}


double DileptonTrees::SimWeight(LeptonCat leptonCat) const {
  double weight = 1.;
  weight *= (*genWeight_)() * intLumi_;
  weight *= (*ewCorrectionWeight_)() * (*kFactorCorrection_)();
  weight *= (*pileUpWeight_)();
  weight *= LeptonSFWeight(leptonCat);
  weight *= bTagWeight_(jetBuilder_.Get());

  return weight;
}

