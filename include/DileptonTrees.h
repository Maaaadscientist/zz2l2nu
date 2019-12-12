#ifndef HZZ2L2NU_INCLUDE_DILEPTONTREES_H_
#define HZZ2L2NU_INCLUDE_DILEPTONTREES_H_

#include <array>
#include <memory>
#include <optional>
#include <tuple>

#include <boost/program_options.hpp>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TTreeReaderValue.h>

#include <AnalysisCommon.h>
#include <Dataset.h>
#include <GenZZBuilder.h>
#include <Options.h>


/**
 * \brief Implements an analysis in the dilepton channel
 *
 * This class applies a slightly looser version of the standard analysis
 * selection in the dilepton channel: the leptons are not required to be of the
 * same flavour, and the cut on ptmiss is loosened. For each selected event, a
 * few observables, together with the event weight, are stored in a ROOT tree.
 * In addition, momenta and other properties of jets and leptons can be stored
 * if flag --more-vars is provided.
 */
class DileptonTrees : public AnalysisCommon {
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

  /**
   * \brief Performs selection on leptons
   *
   * \return If the current event does not pass the selection, the returned
   * optional is empty. Otherwise it contains a tuple of the determined lepton
   * category and pointers to the two leptons, which are ordered by pt.
   */
  std::optional<std::tuple<LeptonCat, Lepton const *, Lepton const *>>
  CheckLeptons() const;

  /// Fills additional variables, mostly lepton and jet momenta
  void FillMoreVariables(std::array<Lepton, 2> const &leptons,
      std::vector<Jet> const &jets);
  
  static double constexpr kNominalMZ_ = 91.1876;

  Dataset &dataset_;

  /// Indicates that additional variables should be stored
  bool storeMoreVariables_;

  /**
   * \brief An object to reconstruct generator-level ZZ system
   *
   * Only created for datasets for ZZ production with decays to 2l2nu.
   */
  std::unique_ptr<GenZZBuilder> genZZBuilder_;

  TTreeReaderValue<ULong64_t> srcEvent_;

  TFile outputFile_;
  TTree *tree_;

  Int_t leptonCat_, jetCat_;
  TLorentzVector *p4LL_, *p4Miss_;
  Float_t mT_;

  ULong64_t event_;
  Float_t genMZZ_;
  Int_t leptonCharge_[2];
  Float_t leptonPt_[2], leptonEta_[2], leptonPhi_[2], leptonMass_[2];
  static int const maxSize_ = 32;
  Int_t jetSize_;
  Float_t jetPt_[maxSize_], jetEta_[maxSize_], jetPhi_[maxSize_],
          jetMass_[maxSize_];

  Float_t weight_;
};

#endif  // HZZ2L2NU_INCLUDE_DILEPTONTREES_H_

