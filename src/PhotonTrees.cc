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
      photonBuilder_{dataset, options},
      photonPrescales_{dataset, options},
      photonWeight_{dataset, options, &photonBuilder_},
      p4Photon_{nullptr}, p4Miss_{nullptr} {

  photonBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});
  jetBuilder_.EnableCleaning({&photonBuilder_});
  ptMissBuilder_.PullCalibration({&photonBuilder_});

  weightCollector_.Add(&photonWeight_);

  CreateWeightBranches();

  AddBranch("jetCat", &jetCat_);
  AddBranch("p4Photon", &p4Photon_);
  AddBranch("p4Miss", &p4Miss_);
  AddBranch("mT", &mT_);
  AddBranch("triggerWeight", &triggerWeight_);

  if (storeMoreVariables_) {
    AddBranch("event", &event_);

    AddBranch("jet_size", &jetSize_);
    AddBranch("jet_pt", jetPt_, "jet_pt[jet_size]/F");
    AddBranch("jet_eta", jetEta_, "jet_eta[jet_size]/F");
    AddBranch("jet_phi", jetPhi_, "jet_phi[jet_size]/F");
    AddBranch("jet_mass", jetMass_, "jet_mass[jet_size]/F");
  }

}


po::options_description PhotonTrees::OptionsDescription() {
  po::options_description optionsDescription{"Analysis-specific options"};
  optionsDescription.add_options()
    ("syst", po::value<std::string>()->default_value(""),
     "Requested systematic variation")
    ("output,o", po::value<std::string>()->default_value("output.root"),
     "Name for output file with histograms")
    ("more-vars", "Store additional variables");
  return optionsDescription;
}


bool PhotonTrees::ProcessEvent() {
  if (not meKinFilter_() or not metFilters_())
    return false;

  auto const photonResult = CheckPhotons();
  if (not photonResult)
    return false;
  auto const photon = photonResult.value();

  *p4Photon_ = photon->p4;

  if (p4Photon_->Pt() < 55.)
    return false;


  *p4Miss_ = ptMissBuilder_.Get().p4;

  if (std::abs(TVector2::Phi_mpi_pi(p4Photon_->Phi() - p4Miss_->Phi())) < 0.5)
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
  else if (utils::PassVbfCuts(jets, *p4Photon_))
    jetCat_ = int(JetCat::kVbf);
  else
    jetCat_ = int(JetCat::kGEq1J);


  // FIXME eT definition should include mass from the mass lineshape, computed separately.
  double const eT =
      std::sqrt(std::pow(p4Photon_->Pt(), 2) + std::pow(p4Photon_->M(), 2))
      + std::sqrt(std::pow(p4Miss_->Pt(), 2) + std::pow(kNominalMZ_, 2));
  mT_ = std::sqrt(std::pow(eT, 2) - std::pow((*p4Photon_ + *p4Miss_).Pt(), 2));

  triggerWeight_ = photonPrescales_.GetWeight(p4Photon_->Pt());
  if (triggerWeight_ == 0)
    return false;
  

  if (storeMoreVariables_)
    FillMoreVariables(jets);

  FillTree();
  return true;

}

std::optional<Photon const *> PhotonTrees::CheckPhotons() const {
  auto const &looseElectrons = electronBuilder_.GetLoose();
  auto const &looseMuons = muonBuilder_.GetLoose();
  auto const &photons = photonBuilder_.Get();

  if (looseElectrons.size() + looseMuons.size() > 0 or photons.size() != 1) {
    return {};
  }
  else {
    return &photons[0];
  }
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
// - prescales: see how we handle them
// - >= 1 PV ?
// - MC QCD
// - ZGamma overlap
// - Add the number of vertices to the list of branches
// - Application of weights
// - Give mass to photon
// - Computation of weights, uncertainties etc...
