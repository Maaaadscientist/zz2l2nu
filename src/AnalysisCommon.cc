#include <AnalysisCommon.h>


AnalysisCommon::AnalysisCommon(Options const &options, Dataset &dataset)
    : intLumi_{options.GetConfig()["luminosity"].as<double>()},
      randomGenerator_{options.GetAs<unsigned>("seed")},
      bTagger_{options},
      electronBuilder_{dataset, options},
      muonBuilder_{dataset, options, randomGenerator_},
      jetBuilder_{dataset, options, randomGenerator_},
      ptMissBuilder_{dataset},
      meKinFilter_{dataset}, metFilters_{dataset},
      bTagWeight_{options, bTagger_} {

  bool const isSim = dataset.Info().IsSimulation();

  if (isSim) {
    genJetBuilder_.reset(new GenJetBuilder(dataset, options));
    jetBuilder_.SetGenJetBuilder(genJetBuilder_.get());
  }

  jetBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});

  // Type 1 corrections for ptmiss
  ptMissBuilder_.PullCalibration(
      {&muonBuilder_, &electronBuilder_, &jetBuilder_});

  if (isSim) {
    genWeight_.reset(new GenWeight{dataset});
    ewCorrectionWeight_.reset(new EWCorrectionWeight{dataset, options});
    pileUpWeight_.reset(new PileUpWeight{dataset, options});
    kFactorCorrection_.reset(new KFactorCorrection{dataset, options});
  }
}
