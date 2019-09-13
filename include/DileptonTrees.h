#ifndef HZZ2L2NU_INCLUDE_DILEPTONTREES_H_
#define HZZ2L2NU_INCLUDE_DILEPTONTREES_H_

#include <memory>

#include <boost/program_options.hpp>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TTree.h>

#include <BTagger.h>
#include <BTagWeight.h>
#include <Dataset.h>
#include <ElectronBuilder.h>
#include <EWCorrectionWeight.h>
#include <GenJetBuilder.h>
#include <GenWeight.h>
#include <GenZZBuilder.h>
#include <JetBuilder.h>
#include <KFactorCorrection.h>
#include <MeKinFilter.h>
#include <MetFilters.h>
#include <MuonBuilder.h>
#include <Options.h>
#include <PileUpWeight.h>
#include <PtMissBuilder.h>


/**
 * \brief Implements an analysis in the dilepton channel
 *
 * This class applies a slightly looser version of the standard analysis
 * selection in the dilepton channel: the leptons are not required to be of the
 * same flavour, and the cut on ptmiss is loosened. For each selected event, a
 * few observables, together with the event weight, are stored in a ROOT tree.
 */
class DileptonTrees {
 public:
  DileptonTrees(Options const &options, Dataset &dataset);

  /// Constructs descriptions for command line options
  static boost::program_options::options_description OptionsDescription();

  /// Writes the output file
  void PostProcessing();

  /// Performs the event selection and fills the output tree
  bool ProcessEvent();

 private:
  enum class LeptonCat : int {
    kEE = 0,
    kMuMu = 1,
    kEMu = 2
  };

  enum class JetCat : int {
    kEq0J,
    kGEq1J,
    kVbf
  };

  /// Computes weight to account for lepton scale factors
  double LeptonSFWeight(LeptonCat leptonCat) const;

  /// Computes full weight to be applied to simulation
  double SimWeight(LeptonCat leptonCat) const;

  double const kNominalMZ_ = 91.1876;

  Dataset &dataset_;
  double intLumi_;

  TRandom3 randomGenerator_;
  BTagger bTagger_;

  /**
   * \brief An object to reconstruct generator-level ZZ system
   *
   * Only created for datasets for ZZ production with decays to 2l2nu.
   */
  std::unique_ptr<GenZZBuilder> genZZBuilder_;

  ElectronBuilder electronBuilder_;
  MuonBuilder muonBuilder_;
  std::unique_ptr<GenJetBuilder> genJetBuilder_;
  JetBuilder jetBuilder_;
  PtMissBuilder ptMissBuilder_;

  MeKinFilter meKinFilter_;
  MetFilters metFilters_;

  std::unique_ptr<GenWeight> genWeight_;
  std::unique_ptr<EWCorrectionWeight> ewCorrectionWeight_;
  std::unique_ptr<PileUpWeight> pileUpWeight_;
  std::unique_ptr<KFactorCorrection> kFactorCorrection_;
  BTagWeight bTagWeight_;

  TFile outputFile_;
  TTree *tree_;
  Int_t leptonCat_, jetCat_;
  TLorentzVector *p4LL_, *p4Miss_;
  Float_t mT_;
  Float_t genMZZ_;
  Float_t weight_;
};

#endif  // HZZ2L2NU_INCLUDE_DILEPTONTREES_H_

