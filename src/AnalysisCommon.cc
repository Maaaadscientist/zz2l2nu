#include <AnalysisCommon.h>

#include <initializer_list>


AnalysisCommon::AnalysisCommon(Options const &options, Dataset &dataset)
    : intLumi_{options.GetConfig()["luminosity"].as<double>()},
      tabulatedRngEngine_{dataset},
      bTagger_{options},
      electronBuilder_{dataset, options},
      muonBuilder_{dataset, options, tabulatedRngEngine_},
      jetBuilder_{dataset, options, tabulatedRngEngine_},
      ptMissBuilder_{dataset, options},
      meKinFilter_{dataset}, metFilters_{dataset},
      leptonWeight_{dataset, options, &electronBuilder_, &muonBuilder_},
      bTagWeight_{dataset, options, &bTagger_, &jetBuilder_} {

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
    genWeight_.reset(new GenWeight{dataset, options});
    ewCorrectionWeight_.reset(new EWCorrectionWeight{dataset, options});
    pileUpWeight_.reset(new PileUpWeight{dataset, options});
    l1tPrefiringWeight_.reset(new L1TPrefiringWeight{dataset, options});
    kFactorCorrection_.reset(new KFactorCorrection{dataset, options});

    weightCollector_.Add(genWeight_.get());
    weightCollector_.Add(ewCorrectionWeight_.get());
    weightCollector_.Add(pileUpWeight_.get());
    weightCollector_.Add(l1tPrefiringWeight_.get());
    weightCollector_.Add(kFactorCorrection_.get());
    weightCollector_.Add(&leptonWeight_);
    weightCollector_.Add(&bTagWeight_);
  }
}

