#include <DileptonTrees.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include <TVector2.h>

#include <Utils.h>


namespace po = boost::program_options;


int const DileptonTrees::maxSize_;


DileptonTrees::DileptonTrees(Options const &options, Dataset &dataset)
    : AnalysisCommon{options, dataset},
      dataset_{dataset},
      storeMoreVariables_{options.Exists("more-vars")},
      srcEvent_{dataset_.Reader(), "event"},
      outputFile_{options.GetAs<std::string>("output").c_str(), "recreate"},
      p4LL_{nullptr}, p4Miss_{nullptr},
      weight_{1.} {
        
  if (dataset_.Info().IsSimulation()) {
    auto const &node = dataset_.Info().Parameters()["zz_2l2nu"];

    if (node and not node.IsNull() and node.as<bool>())
      genZZBuilder_.reset(new GenZZBuilder(dataset));
  }


  tree_ = new TTree("Vars", "Variables in dilepton selection");
  tree_->SetDirectory(&outputFile_);
  tree_->Branch("leptonCat", &leptonCat_);
  tree_->Branch("jetCat", &jetCat_);
  tree_->Branch("p4LL", &p4LL_);
  tree_->Branch("p4Miss", &p4Miss_);
  tree_->Branch("mT", &mT_);

  if (storeMoreVariables_) {
    tree_->Branch("event", &event_);

    if (genZZBuilder_)
      tree_->Branch("genMZZ", &genMZZ_);

    tree_->Branch("lepton_charge", leptonCharge_, "lepton_charge[2]/I");
    tree_->Branch("lepton_pt", leptonPt_, "lepton_pt[2]/F");
    tree_->Branch("lepton_eta", leptonEta_, "lepton_eta[2]/F");
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
    weight_ = SimWeight();

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
  
  event_ = *srcEvent_;

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

double DileptonTrees::SimWeight() const {
  double weight = 1.;
  weight *= (*genWeight_)() * intLumi_;
  weight *= (*ewCorrectionWeight_)() * (*kFactorCorrection_)();
  weight *= (*pileUpWeight_)();
  weight *= leptonWeight_();
  weight *= bTagWeight_();

  return weight;
}

