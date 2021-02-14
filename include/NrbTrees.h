#ifndef HZZ2L2NU_INCLUDE_NRBTREES_H_
#define HZZ2L2NU_INCLUDE_NRBTREES_H_

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
 * \brief Implements an analysis for non-resonant background 
 *
 * This class applies a looser version with control region of the standard 
 * analysis selection in the dilepton channel: 
 * the leptons are not required to be of the same flavour, 
 * and the cuts on ptmiss, btag, ll mass and ll pt are loosened . 
 * For each selected event, a few observables, 
 * together with the event weight, are stored in a ROOT tree.
 * In addition, momenta and other properties of jets and leptons can be stored
 * if flag --more-vars is provided.
 *
 * Event weights are saved with the help of the base class EventTrees.
 */
class NrbTrees final : public EventTrees {
 public:
  NrbTrees(Options const &options, Dataset &dataset);

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
   * category and pointers to the two leptons, which are ordered by pt in the 
	 * flavor-same channel, with electron as l1 and muon as l2 in flavor-opposite
	 * channel.
   */
  std::optional<std::tuple<LeptonCat, Lepton const *, Lepton const *>>
  CheckLeptons() const;

  /// Fills additional variables, mostly lepton and jet momenta
  void FillMoreVariables(std::array<Lepton, 2> const &leptons,
      std::vector<Jet> const &jets);

  static double constexpr kNominalMZ_ = 91.1876;

  /// Indicates that additional variables should be stored
  bool storeMoreVariables_;

  /// Specifies the cut on ptmiss. Default is 50.
  double ptMissCut_;

  /**
   * \brief An object to reconstruct generator-level ZZ system
   *
   * Only created for datasets for ZZ production with decays to 2l2nu.
   */
  std::optional<GenZZBuilder> genZZBuilder_;

  /// Buffer to save the nominal event weight, transfer factors for emu.
  Float_t weight_;

  /**
   * \brief Buffers to save alternative event weights
   *
   * These are full weights, i.e. they are not relative with respect to the
   * nominal one.
   */
  std::vector<Float_t> leptonWeights_;
  std::vector<Float_t> triggerWeights_;

  TriggerFilter triggerFilter_;

  LeptonWeight leptonEff_;
  TriggerWeight triggerEff_;

  TTreeReaderValue<ULong64_t> srcEvent_;

  Int_t leptonCat_, jetCat_, numPVGood_;
  Float_t llPt_, llEta_, llPhi_, llMass_;
  Float_t missPt_, missPhi_;
  Float_t mT_;
  Float_t l1Pt_, l2Pt_, l1Eta_, l2Eta_;
  Float_t eeLeptonWeight_[9], mumuLeptonWeight_[9];
  Float_t eeTriggerWeight_[7], mumuTriggerWeight_[7];
  Bool_t btagLoose_, btagMedium_, btagTight_;
  Bool_t btagLooseLowPt_, btagMediumLowPt_, btagTightLowPt_;

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

#endif  // HZZ2L2NU_INCLUDE_NRBTREES_H_
