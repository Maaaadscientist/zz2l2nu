#ifndef HZZ2L2NU_INCLUDE_ELECTRONTREES_H_
#define HZZ2L2NU_INCLUDE_ELECTRONTREES_H_

#include <optional>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReaderValue.h>

#include <Dataset.h>
// #include <EventNumberFilter.h>
#include <EventTrees.h>
// #include <GJetsWeight.h>
// #include <GenPhotonBuilder.h>
#include <Options.h>
#include <PhotonBuilder.h>
// #include <PhotonPrescales.h>
// #include <PhotonWeight.h>


class ElectronTrees final : public EventTrees {
 public:
  ElectronTrees(Options const &options, Dataset &dataset);

  /// Constructs descriptions for command line options
  static boost::program_options::options_description OptionsDescription();

  /// Performs the event selection and fills the output tree
  bool ProcessEvent();
  TTreeReaderValue<UInt_t> srcRun_;
  TTreeReaderValue<UInt_t> srcLumi_;
  TTreeReaderValue<ULong64_t> srcEvent_;

 private:
  enum class JetCat : int {
    kEq0J,
    kEq1J,
    kGEq2J
  };

  Electron const *CheckElectron() const;

  void FillMoreVariables(std::vector<Jet> const &jets);

  // static double constexpr kNominalMZ_ = 91.1876;

  /// Indicates that additional variables should be stored
  bool storeMoreVariables_;

  // std::optional<GenPhotonBuilder> genPhotonBuilder_;


  PhotonBuilder photonBuilder_;

  // PhotonPrescales photonPrescales_;
  // std::vector<PhotonTrigger> photonTriggers_;

  // PhotonWeight photonWeight_;

  // GJetsWeight gJetsWeight_;

  // EventNumberFilter photonFilter_;

  // std::string labelWGamma_ = "";
  // std::string labelZGamma_ = "";

  Int_t jetSize_, jetCat_, numPVGood_;
  Float_t electronPt_, electronEta_, electronPhi_;
  Float_t electronM_;
  // Float_t dijetM_;
  Float_t missPt_, missPhi_;
  Float_t electronMetDeltaPhi_;
  Float_t electronMetMt_;
  // Float_t mT_, triggerWeight_, photonReweighting_, photonNvtxReweighting_;
  // Float_t photonEtaReweighting_;
  // Float_t meanWeight_;

  TTreeReaderValue<int> srcNumPVGood_;

  UInt_t run_, lumi_;
  ULong64_t event_;
  static int const maxSize_ = 32;
  Float_t jetPt_[maxSize_], jetEta_[maxSize_], jetPhi_[maxSize_],
          jetMass_[maxSize_];

  // bool isQCD_;
  // std::vector<std::string> lines_;

};

#endif  // HZZ2L2NU_INCLUDE_ELECTRONTREES_H_
