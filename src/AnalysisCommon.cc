#include <AnalysisCommon.h>

#include <initializer_list>

namespace po = boost::program_options;


AnalysisCommon::AnalysisCommon(Options const &options, Dataset &dataset)
    : intLumi_{options.GetConfig()["luminosity"].as<double>()},
      isSim_{dataset.Info().IsSimulation()},
      tabulatedRngEngine_{dataset},
      bTagger_{options}, pileUpIdFilter_{options},
      electronBuilder_{dataset, options},
      muonBuilder_{dataset, options, tabulatedRngEngine_},
      jetBuilder_{dataset, options, tabulatedRngEngine_},
      ptMissBuilder_{dataset, options},
      meKinFilter_{dataset}, metFilters_{options, dataset},
      leptonWeight_{dataset, options, &electronBuilder_, &muonBuilder_},
      bTagWeight_{dataset, options, &bTagger_, &jetBuilder_} {

  YAML::Node const &selectionCutsNode = options.GetConfig()["selection_cuts"];
  zMassWindow_ = selectionCutsNode["z_mass_window"].as<double>();
  minPtLL_ = selectionCutsNode["min_pt_ll"].as<double>();
  minDphiLLPtMiss_ = selectionCutsNode["min_dphi_ll_ptmiss"].as<double>();
  minDphiJetsPtMiss_ = selectionCutsNode["min_dphi_jets_ptmiss"].as<double>();
  minDphiLeptonsJetsPtMiss_ =
    selectionCutsNode["min_dphi_leptonsjets_ptmiss"].as<double>();

  if (isSim_) {
    genJetBuilder_.emplace(dataset, options);
    jetBuilder_.SetGenJetBuilder(&genJetBuilder_.value());
  }

  jetBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});
  jetBuilder_.SetPileUpIdFilter(&pileUpIdFilter_);

  // Type 1 corrections for ptmiss
  ptMissBuilder_.PullCalibration(
      {&muonBuilder_, &electronBuilder_, &jetBuilder_});

  if (isSim_) {
    genWeight_.emplace(dataset, options);
    kFactorCorrection_.emplace(dataset, options);
    ewCorrectionWeight_.emplace(dataset, options);
    pileUpWeight_.emplace(dataset, options);
    l1tPrefiringWeight_.emplace(dataset, options);

    weightCollector_.Add(&genWeight_.value());
    weightCollector_.Add(&kFactorCorrection_.value());
    weightCollector_.Add(&ewCorrectionWeight_.value());
    weightCollector_.Add(&pileUpWeight_.value());
    weightCollector_.Add(&l1tPrefiringWeight_.value());
    weightCollector_.Add(&leptonWeight_);
    weightCollector_.Add(&bTagWeight_);
  }
}


po::options_description AnalysisCommon::OptionsDescription() {
  po::options_description optionsDescription{"Analysis-specific options"};
  optionsDescription.add_options()
    ("syst", po::value<std::string>()->default_value(""),
     "Requested systematic variation")
    ("output,o", po::value<std::string>()->default_value("output.root"),
     "Name for output file with histograms");
  return optionsDescription;
}


double AnalysisCommon::DPhiLeptonsJetsSystemPtMiss() {

  auto const &electrons = electronBuilder_.GetTight();
  auto const &muons = muonBuilder_.GetTight();
  auto const &jets = jetBuilder_.Get();
  TLorentzVector const &p4Miss = ptMissBuilder_.Get().p4;
  TLorentzVector p4LeptonsJets(0, 0, 0, 0);

  for (auto const &electron : electrons)
    p4LeptonsJets += electron.p4;

  for (auto const &muon : muons)
    p4LeptonsJets += muon.p4;

  for (auto const &jet : jets)
    p4LeptonsJets += jet.p4;

  return std::abs(TVector2::Phi_mpi_pi(p4LeptonsJets.Phi() - p4Miss.Phi()));
}

