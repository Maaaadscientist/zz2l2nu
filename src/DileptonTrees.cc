#include <DileptonTrees.h>

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
        
  if (dataset_.Info().IsSimulation()) {
    genJetBuilder_.reset(new GenJetBuilder(dataset_, options));
    jetBuilder_.SetGenJetBuilder(genJetBuilder_.get());
  }
  
  jetBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});
  
  // Type 1 corrections for ptmiss
  ptMissBuilder_.PullCalibration(
      {&muonBuilder_, &electronBuilder_, &jetBuilder_});

  if (dataset_.Info().IsSimulation()) {
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
     "Seed for random number generator; 0 means a unique seed");
  return optionsDescription;
}


void DileptonTrees::PostProcessing() {
  outputFile_.Write();
  outputFile_.Close();
}


bool DileptonTrees::ProcessEvent() {
  if (not meKinFilter_() or not metFilters_())
    return false;


  auto const &tightElectrons = electronBuilder_.GetTight();
  auto const &looseElectrons = electronBuilder_.GetLoose();

  auto const &tightMuons = muonBuilder_.GetTight();
  auto const &looseMuons = muonBuilder_.GetLoose();

  if (looseElectrons.size() + looseMuons.size() != 2)
    return false;

  LeptonCat leptonCat;

  if (tightElectrons.size() == 2) {
    leptonCat = LeptonCat::kEE;
    *p4LL_ = tightElectrons[0].p4 + tightElectrons[1].p4;
  } else if (tightMuons.size() == 2) {
    leptonCat = LeptonCat::kMuMu;
    *p4LL_ = tightMuons[0].p4 + tightMuons[1].p4;
  } else if (tightElectrons.size() == 1 and tightMuons.size() == 1) {
    leptonCat = LeptonCat::kEMu;
    *p4LL_ = tightElectrons[0].p4 + tightMuons[0].p4;
  } else
    return false;

  leptonCat_ = int(leptonCat);

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
  

  if (dataset_.Info().IsSimulation())
    weight_ = SimWeight(leptonCat);

  tree_->Fill();
  return true;
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

