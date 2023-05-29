#ifndef HZZ2L2NU_INCLUDE_ZGAMMATREES_H_
#define HZZ2L2NU_INCLUDE_ZGAMMATREES_H_

#include <optional>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <Dataset.h>
#include <EventNumberFilter.h>
#include <EventTrees.h>
#include <GJetsWeight.h>
#include <GenPhotonBuilder.h>
#include <Options.h>
#include <PhotonBuilder.h>
#include <PhotonPrescales.h>
#include <PhotonWeight.h>
#include <TriggerFilter.h>


class ZGammaTrees final : public EventTrees {
 public:
  ZGammaTrees(Options const &options, Dataset &dataset);

  /// Constructs descriptions for command line options
  static boost::program_options::options_description OptionsDescription();

  /// Performs the event selection and fills the output tree
  bool ProcessEvent();

 private:
  enum class LeptonCat : int {
    kEE = 0,
    kMuMu = 1,
    // kEMu = 2
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

  Photon const *CheckPhotons() const;

  void FillMoreVariables(std::vector<Jet> const &jets);

  static double constexpr kNominalMZ_ = 91.1876;

  /// Indicates that additional variables should be stored
  bool storeMoreVariables_;

  std::optional<GenPhotonBuilder> genPhotonBuilder_;

  TTreeReaderValue<UInt_t> srcRun_;
  TTreeReaderValue<UInt_t> srcLumi_;
  TTreeReaderValue<ULong64_t> srcEvent_;

  mutable std::unique_ptr<TTreeReaderValue<UInt_t>> numGenPart_;
  mutable std::unique_ptr<TTreeReaderArray<Int_t>> genPartPdgId_;
  mutable std::unique_ptr<TTreeReaderArray<Float_t>> genPartPt_, genPartEta_, genPartPhi_;
  mutable std::unique_ptr<TTreeReaderArray<Int_t>> genPartStatus_;
  mutable std::unique_ptr<TTreeReaderArray<Int_t>> genPartStatusFlags_;

  PhotonBuilder photonBuilder_;

  PhotonPrescales photonPrescales_;
  std::vector<PhotonTrigger> photonTriggers_;

  PhotonWeight photonWeight_;

  TriggerFilter triggerFilter_;
  GJetsWeight gJetsWeight_;

  //EventNumberFilter photonFilter_;

  std::string labelWGamma_ = "";
  std::string labelZGamma_ = "";
  std::string datasetName_ = "";

  Int_t leptonCat_, jetCat_, analysisCat_, numPVGood_;
  Float_t llPt_, llEta_, llPhi_, llMass_;
  Float_t photonPt_, photonEta_, photonPhi_, photonMass_;
  Float_t photonR9_, photonSieie_;
  Float_t missPt_, missPhi_;
  Float_t mT_, triggerWeight_;
  // Float_t photonReweighting_, photonNvtxReweighting_, photonEtaReweighting_;
  // Float_t meanWeight_;
  Float_t l1Pt_, l2Pt_;
  Bool_t isOverlapped_;

  TTreeReaderValue<int> srcNumPVGood_;

  UInt_t run_, lumi_;
  ULong64_t event_;
  static int const maxSize_ = 32;
  Int_t jetSize_;
  Float_t jetPt_[maxSize_], jetEta_[maxSize_], jetPhi_[maxSize_],
          jetMass_[maxSize_];

  bool isQCD_;

  // FIXME temporary. These will be replaced by a new class, much more practical. For now, still use old functions from Utils.
  std::vector<std::string> v_jetCat_, v_analysisCat_;
  bool applyNvtxWeights_, applyEtaWeights_, applyPtWeights_, applyMassLineshape_;
  bool applyMeanWeights_;
  std::map<TString, std::map<double, std::pair<double, double>>> nVtxWeight_map_;
  std::map<TString, std::map<double, std::pair<double, double>>> etaWeight_map_;
  std::map<TString, std::map<double, std::pair<double, double>>> ptWeight_map_;
  std::map<TString, std::map<double, double>> meanWeight_map_;
  std::map<TString, TH1 *> lineshapeMassWeight_map_;
};

#endif  // HZZ2L2NU_INCLUDE_ZGAMMATREES_H_
