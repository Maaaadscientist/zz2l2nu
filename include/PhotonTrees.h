#ifndef HZZ2L2NU_INCLUDE_PHOTONTREES_H_
#define HZZ2L2NU_INCLUDE_PHOTONTREES_H_

#include <optional>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TTreeReaderValue.h>

#include <Dataset.h>
#include <EventTrees.h>
#include <GJetsWeight.h>
#include <Options.h>
#include <PhotonBuilder.h>
#include <PhotonPrescales.h>
#include <PhotonWeight.h>


class PhotonTrees final : public EventTrees {
 public:
  PhotonTrees(Options const &options, Dataset &dataset);

  /// Constructs descriptions for command line options
  static boost::program_options::options_description OptionsDescription();

  /// Performs the event selection and fills the output tree
  bool ProcessEvent();

 private:
  enum class JetCat : int {
    kEq0J,
    kGEq1J,
    kVbf
  };

  Photon const *CheckPhotons() const;

  void FillMoreVariables(std::vector<Jet> const &jets);

  static double constexpr kNominalMZ_ = 91.1876;

  /// Indicates that additional variables should be stored
  bool storeMoreVariables_;

  TTreeReaderValue<ULong64_t> srcEvent_;

  PhotonBuilder photonBuilder_;

  PhotonPrescales photonPrescales_;
  std::vector<PhotonTrigger> photonTriggers_;

  PhotonWeight photonWeight_;

  GJetsWeight gJetsWeight_;

  std::string labelZGamma_ = "";

  Int_t jetCat_, numPVGood_;
  TLorentzVector *p4Photon_, *p4Miss_;
  Float_t mT_, triggerWeight_, photonReweighting_;

  TTreeReaderValue<int> srcNumPVGood_;

  ULong64_t event_;
  static int const maxSize_ = 32;
  Int_t jetSize_;
  Float_t jetPt_[maxSize_], jetEta_[maxSize_], jetPhi_[maxSize_],
          jetMass_[maxSize_];

  bool isQCD_;

  // FIXME temporary. These will be replaced by a new class, much more practical. For now, still use old functions from Utils.
  std::vector<std::string> v_jetCat_;
  bool weight_NVtx_exist_, weight_Pt_exist_, weight_Mass_exist_;
  std::map<TString, std::map<std::pair<double, double>, std::pair<double, double>>> nVtxWeight_map_;
  std::map<TString, std::map<double, std::pair<double, double>>> ptWeight_map_;
  std::map<TString, TH1 *> lineshapeMassWeight_map_;
};

#endif  // HZZ2L2NU_INCLUDE_PHOTONTREES_H_