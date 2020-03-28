#include <AnalysisCommon.h>

#include <initializer_list>


AnalysisCommon::AnalysisCommon(Options const &options, Dataset &dataset)
    : intLumi_{options.GetConfig()["luminosity"].as<double>()},
      isSim_{dataset.Info().IsSimulation()},
      tabulatedRngEngine_{dataset},
      bTagger_{options},
      electronBuilder_{dataset, options},
      muonBuilder_{dataset, options, tabulatedRngEngine_},
      jetBuilder_{dataset, options, tabulatedRngEngine_},
      ptMissBuilder_{dataset, options},
      meKinFilter_{dataset}, metFilters_{dataset},
      leptonWeight_{dataset, options, &electronBuilder_, &muonBuilder_},
      bTagWeight_{dataset, options, &bTagger_, &jetBuilder_} {

  YAML::Node const &selectionCutsNode = options.GetConfig()["selection_cuts"];
  zMassWindow_ = selectionCutsNode["z_mass_window"].as<double>();
  minPtLL_ = selectionCutsNode["min_pt_ll"].as<double>();
  minDphiLLPtMiss_ = selectionCutsNode["min_dphi_ll_ptmiss"].as<double>();
  minDphiJetPtMiss_ = selectionCutsNode["min_dphi_jet_ptmiss"].as<double>();

  if (isSim_) {
    genJetBuilder_.emplace(dataset, options);
    jetBuilder_.SetGenJetBuilder(&genJetBuilder_.value());
  }

  jetBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});

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

