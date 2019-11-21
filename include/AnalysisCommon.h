#ifndef HZZ2L2NU_INCLUDE_ANALYSISCOMMON_H_
#define HZZ2L2NU_INCLUDE_ANALYSISCOMMON_H_

#include <memory>

#include <BTagger.h>
#include <BTagWeight.h>
#include <Dataset.h>
#include <ElectronBuilder.h>
#include <EWCorrectionWeight.h>
#include <GenJetBuilder.h>
#include <GenWeight.h>
#include <JetBuilder.h>
#include <KFactorCorrection.h>
#include <LeptonWeight.h>
#include <MeKinFilter.h>
#include <MetFilters.h>
#include <MuonBuilder.h>
#include <Options.h>
#include <PileUpWeight.h>
#include <PtMissBuilder.h>
#include <TabulatedRandomGenerator.h>


/**
 * \brief Implements common aspects of analyses
 *
 * This class creates a number of blocks used in multiple analyses, thus
 * avoiding code duplication. Provided blocks include
 * - builders for leptons, jets, and missing pt,
 * - common event filters,
 * - common reweighting objects,
 * and some others. They are provided to the inheriting classes as protected
 * data members.
 */
class AnalysisCommon {
 public:
  AnalysisCommon(Options const &options, Dataset &dataset);

 protected:
  /// Integrated luminosity, 1/pb
  double intLumi_;

  /// Common random number generator engine
  TabulatedRngEngine tabulatedRngEngine_;

  BTagger bTagger_;

  ElectronBuilder electronBuilder_;
  MuonBuilder muonBuilder_;
  std::unique_ptr<GenJetBuilder> genJetBuilder_;
  JetBuilder jetBuilder_;
  PtMissBuilder ptMissBuilder_;

  LeptonWeight leptonWeight_;
  MeKinFilter meKinFilter_;
  MetFilters metFilters_;

  std::unique_ptr<GenWeight> genWeight_;
  std::unique_ptr<EWCorrectionWeight> ewCorrectionWeight_;
  std::unique_ptr<PileUpWeight> pileUpWeight_;
  std::unique_ptr<KFactorCorrection> kFactorCorrection_;
  BTagWeight bTagWeight_;
};

#endif  // HZZ2L2NU_INCLUDE_ANALYSISCOMMON_H_

