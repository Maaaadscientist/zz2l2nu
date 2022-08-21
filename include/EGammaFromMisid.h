#ifndef HZZ2L2NU_INCLUDE_EGAMMAFROMMISID_H_
#define HZZ2L2NU_INCLUDE_EGAMMAFROMMISID_H_

#include <array>
#include <optional>
#include <tuple>
#include <vector>

#include <boost/program_options.hpp>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReaderValue.h>

#include <EventTrees.h>
#include <Dataset.h>
#include <GenZZBuilder.h>
#include <Options.h>
#include <TriggerFilter.h>


/**
 * \brief Implements an analysis in the dilepton channel
 *
 * This class applies a slightly looser version of the standard analysis
 * selection in the dilepton channel: the leptons are not required to be of the
 * same flavour, and the cut on ptmiss is loosened. For each selected event, a
 * few observables, together with the event weight, are stored in a ROOT tree.
 * In addition, momenta and other properties of jets and leptons can be stored
 * if flag --more-vars is provided.
 *
 * Event weights are saved with the help of the base class EventTrees.
 */
class EGammaFromMisid final : public EventTrees {
 public:
  EGammaFromMisid(Options const &options, Dataset &dataset);

  /// Constructs descriptions for command line options
  static boost::program_options::options_description OptionsDescription();

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
    kEq1J,
    kGEq2J
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

  /// Indicates that additional variables should be stored
  bool storeMoreVariables_;

  /// Specifies the cut on ptmiss. Default is 80.
  double ptMissCut_;

  /**
   * \brief An object to reconstruct generator-level ZZ system
   *
   * Only created for datasets for ZZ production with decays to 2l2nu.
   */
  std::optional<GenZZBuilder> genZZBuilder_;

  TriggerFilter triggerFilter_;

  TTreeReaderValue<ULong64_t> srcEvent_;

  Int_t leptonCat_, jetCat_, numPVGood_;
  Float_t llPt_, llEta_, llPhi_, llMass_;
  Float_t missPt_, missPhi_;
  Float_t mT_;

  TTreeReaderValue<int> srcNumPVGood_;

  ULong64_t event_;
  Float_t genMZZ_;
  Int_t leptonCharge_[2];
  Float_t leptonPt_[2], leptonEta_[2], leptonPhi_[2], leptonMass_[2];
  static int const maxSize_ = 32;
  Int_t jetSize_;
  Float_t jetPt_[maxSize_], jetEta_[maxSize_], jetPhi_[maxSize_],
          jetMass_[maxSize_];
};

#endif  // HZZ2L2NU_INCLUDE_EGAMMAFROMMISID_H_
