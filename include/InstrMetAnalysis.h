#ifndef HZZ2L2NU_INCLUDE_INSTRMETANALYSIS_H_
#define HZZ2L2NU_INCLUDE_INSTRMETANALYSIS_H_

#include <map>
#include <list>
#include <memory>
#include <utility>
#include <vector>

#include <boost/program_options.hpp>
#include <TH1.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <AnalysisCommon.h>
#include <Dataset.h>
#include <Options.h>
#include <PhotonBuilder.h>
#include <PhotonPrescales.h>
#include <SmartSelectionMonitor_hzz.h>


/// Constructs weights for data-driven modelling of Z+jets background
class InstrMetAnalysis : public AnalysisCommon {
 public:
  InstrMetAnalysis(Options const &options, Dataset &dataset);

  /// Constructs descriptions for command line options
  static boost::program_options::options_description OptionsDescription();

  /// Saves histograms to the output ROOT file
  void PostProcessing();

  /**
   * \brief Processes the current event in the dataset
   *
   * \return True if the event has passed the full selection, false otherwise.
   */
  bool ProcessEvent();

 private:
  enum {ee, mumu, ll, lepCat_size};
  enum {eq0jets, geq1jets, vbf, jetCat_size};

  Dataset &dataset_;
  bool isMC_;
  std::string outputFile_;
  std::string syst_;

  PhotonBuilder photonBuilder_;

  PhotonPrescales photonPrescales_;
  std::vector<PhotonTrigger> photonTriggers_;

  mutable SmartSelectionMonitor_hzz mon_;
  bool divideFinalHistoByBinWidth_;
  std::vector<std::string> v_jetCat_, tagsR_;
  unsigned tagsR_size_;

  TString fileName_;
  bool isMC_QCD_, isMC_GJet_HT_, isMC_LO_ZNuNuGJets_,
      isMC_NLO_ZGTo2NuG_inclusive_, isMC_NLO_ZGTo2NuG_Pt130_;

  bool weight_NVtx_exist_, weight_Pt_exist_, weight_Mass_exist_;
  std::map<TString, std::map<std::pair<double, double>, std::pair<double, double>>> nVtxWeight_map_;
  std::map<TString, std::map<double, std::pair<double, double>>> ptWeight_map_;
  std::map<TString, TH1 *> lineshapeMassWeight_map_;

  TTreeReaderValue<UInt_t> run_ = {dataset_.Reader(), "run"};
  TTreeReaderValue<UInt_t> luminosityBlock_ = {dataset_.Reader(), "luminosityBlock"};
  TTreeReaderValue<ULong64_t> eventNumber_ = {dataset_.Reader(), "event"};
  TTreeReaderValue<Float_t> rho_ = {dataset_.Reader(), "fixedGridRhoFastjetAll"};
  TTreeReaderValue<Int_t> numPVGood_ = {dataset_.Reader(), "PV_npvsGood"};
  TTreeReaderArray<Float_t> muonPt_ = {dataset_.Reader(), "Muon_pt"};
  TTreeReaderArray<Float_t> electronPt_ = {dataset_.Reader(), "Electron_pt"};
  std::unique_ptr<TTreeReaderArray<Float_t>> genPartPt_, genPartEta_,
      genPartPhi_, genPartMass_;
  std::unique_ptr<TTreeReaderArray<Int_t>> genPartPdgId_, genPartMotherIndex_;
  std::unique_ptr<TTreeReaderArray<Int_t>> photonGenPartIndex_;
  std::unique_ptr<TTreeReaderArray<UChar_t>> photonGenPartFlavor_;
};

#endif  // HZZ2L2NU_INCLUDE_INSTRMETANALYSIS_H_

