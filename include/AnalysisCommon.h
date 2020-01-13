#ifndef HZZ2L2NU_INCLUDE_ANALYSISCOMMON_H_
#define HZZ2L2NU_INCLUDE_ANALYSISCOMMON_H_

#include <optional>

#include <BTagger.h>
#include <BTagWeight.h>
#include <Dataset.h>
#include <ElectronBuilder.h>
#include <EWCorrectionWeight.h>
#include <GenJetBuilder.h>
#include <GenWeight.h>
#include <JetBuilder.h>
#include <KFactorCorrection.h>
#include <L1TPrefiringWeight.h>
#include <LeptonWeight.h>
#include <MeKinFilter.h>
#include <MetFilters.h>
#include <MuonBuilder.h>
#include <Options.h>
#include <PhotonPrescales.h>
#include <PileUpWeight.h>
#include <PtMissBuilder.h>
#include <TabulatedRandomGenerator.h>
#include <WeightCollector.h>


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
  std::optional<GenJetBuilder> genJetBuilder_;
  JetBuilder jetBuilder_;
  PtMissBuilder ptMissBuilder_;

  PhotonPrescales photonPrescales_;

  LeptonWeight leptonWeight_;
  MeKinFilter meKinFilter_;
  MetFilters metFilters_;

  std::optional<GenWeight> genWeight_;
  std::optional<KFactorCorrection> kFactorCorrection_;
  std::optional<EWCorrectionWeight> ewCorrectionWeight_;
  std::optional<PileUpWeight> pileUpWeight_;
  std::optional<L1TPrefiringWeight> l1tPrefiringWeight_;
  BTagWeight bTagWeight_;
  WeightCollector weightCollector_;
};

#endif  // HZZ2L2NU_INCLUDE_ANALYSISCOMMON_H_

