#include <DileptonTrees.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include <TLorentzVector.h>
#include <TVector2.h>

#include <Utils.h>


namespace po = boost::program_options;


int const DileptonTrees::maxSize_;


DileptonTrees::DileptonTrees(Options const &options, Dataset &dataset)
    : EventTrees{options, dataset},
      storeMoreVariables_{options.Exists("more-vars")},
      srcEvent_{dataset.Reader(), "event"} {

  if (isSim_) {
    auto const &node = dataset.Info().Parameters()["zz_2l2nu"];

    if (node and not node.IsNull() and node.as<bool>())
      genZZBuilder_.emplace(dataset);
  }

  CreateWeightBranches();

  AddBranch("lepton_cat", &leptonCat_);
  AddBranch("jet_cat", &jetCat_);
  AddBranch("ll_pt", &llPt_);
  AddBranch("ll_eta", &llEta_);
  AddBranch("ll_phi", &llPhi_);
  AddBranch("ll_mass", &llMass_);
  AddBranch("ptmiss", &missPt_);
  AddBranch("ptmiss_phi", &missPhi_);
  AddBranch("mT", &mT_);

  if (storeMoreVariables_) {
    AddBranch("event", &event_);

    if (genZZBuilder_)
      AddBranch("gen_mzz", &genMZZ_);

    AddBranch("lepton_charge", leptonCharge_, "lepton_charge[2]/I");
    AddBranch("lepton_pt", leptonPt_, "lepton_pt[2]/F");
    AddBranch("lepton_eta", leptonEta_, "lepton_eta[2]/F");
    AddBranch("lepton_phi", leptonPhi_, "lepton_phi[2]/F");
    AddBranch("lepton_mass", leptonMass_, "lepton_mass[2]/F");

    AddBranch("jet_size", &jetSize_);
    AddBranch("jet_pt", jetPt_, "jet_pt[jet_size]/F");
    AddBranch("jet_eta", jetEta_, "jet_eta[jet_size]/F");
    AddBranch("jet_phi", jetPhi_, "jet_phi[jet_size]/F");
    AddBranch("jet_mass", jetMass_, "jet_mass[jet_size]/F");
  }
}


po::options_description DileptonTrees::OptionsDescription() {
  auto optionsDescription = AnalysisCommon::OptionsDescription();
  optionsDescription.add_options()
    ("more-vars", "Store additional variables");
  return optionsDescription;
}


bool DileptonTrees::ProcessEvent() {
  if (not meKinFilter_() or not metFilters_())
    return false;

  auto const leptonResult = CheckLeptons();
  if (not leptonResult)
    return false;

  if (isotrkBuilder_.Get().size() > 0)
    return false;

  auto const &[leptonCat, l1, l2] = leptonResult.value();

  leptonCat_ = int(leptonCat);
  TLorentzVector const p4LL = l1->p4 + l2->p4;
  llPt_ = p4LL.Pt();
  llEta_ = p4LL.Eta();
  llPhi_ = p4LL.Phi();
  llMass_ = p4LL.M();

  if (std::abs(p4LL.M() - kNominalMZ_) > zMassWindow_)
    return false;

  if (p4LL.Pt() < minPtLL_)
    return false;


  auto const &p4Miss = ptMissBuilder_.Get().p4;
  missPt_ = p4Miss.Pt();
  missPhi_ = p4Miss.Phi();

  if (p4Miss.Pt() < 80.)
    return false;

  if (std::abs(
        TVector2::Phi_mpi_pi(p4LL.Phi() - p4Miss.Phi())) < minDphiLLPtMiss_)
    return false;


  auto const &jets = jetBuilder_.Get();

  for (auto const &jet : jets) {
    if (bTagger_(jet))
      return false;

    if (std::abs(TVector2::Phi_mpi_pi(
            jet.p4.Phi() - p4Miss.Phi())) < minDphiJetsPtMiss_)
      return false;
  }

  if (DPhiPtMiss({&jetBuilder_, &muonBuilder_, &electronBuilder_})
      < minDphiLeptonsJetsPtMiss_)
    return false;

  if (jets.size() == 0)
    jetCat_ = int(JetCat::kEq0J);
  else if (utils::PassVbfCuts(jets, p4LL))
    jetCat_ = int(JetCat::kVbf);
  else
    jetCat_ = int(JetCat::kGEq1J);


  double const eT =
      std::sqrt(std::pow(p4LL.Pt(), 2) + std::pow(p4LL.M(), 2))
      + std::sqrt(std::pow(p4Miss.Pt(), 2) + std::pow(kNominalMZ_, 2));
  mT_ = std::sqrt(std::pow(eT, 2) - std::pow((p4LL + p4Miss).Pt(), 2));


  if (storeMoreVariables_)
    FillMoreVariables({*l1, *l2}, jets);

  FillTree();
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
