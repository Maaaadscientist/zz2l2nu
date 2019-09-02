#ifndef HZZ2L2NU_INCLUDE_NRBANALYSIS_H_
#define HZZ2L2NU_INCLUDE_NRBANALYSIS_H_

#include <map>
#include <memory>
#include <utility>
#include <vector>

#include <TH1.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <BTagger.h>
#include <BTagWeight.h>
#include <Dataset.h>
#include <ElectronBuilder.h>
#include <EWCorrectionWeight.h>
#include <GenWeight.h>
#include <GenJetBuilder.h>
#include <JetBuilder.h>
#include <MeKinFilter.h>
#include <MetFilters.h>
#include <MuonBuilder.h>
#include <Options.h>
#include <PhotonBuilder.h>
#include <PileUpWeight.h>
#include <PtMissBuilder.h>
#include <SmartSelectionMonitor_hzz.h>


/// Constructs inputs for data-driven estimation of non-resonant background
class NrbAnalysis {
 public:
  NrbAnalysis(Options const &options, Dataset &dataset);

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

  void InitializeHistograms();
  
  Dataset &dataset_;
  bool isMC_;
  double intLumi_;
  std::string outputFile_;
  bool keepAllControlPlots_;
  bool isPhotonDatadriven_;
  std::string syst_;

  TRandom3 randomGenerator_;
  BTagger bTagger_;

  ElectronBuilder electronBuilder_;
  MuonBuilder muonBuilder_;
  PhotonBuilder photonBuilder_;
  std::unique_ptr<GenJetBuilder> genJetBuilder_;
  JetBuilder jetBuilder_;
  PtMissBuilder ptMissBuilder_;

  MeKinFilter meKinFilter_;
  MetFilters metFilters_;

  std::unique_ptr<GenWeight> genWeight_;
  std::unique_ptr<EWCorrectionWeight> ewCorrectionWeight_;
  std::unique_ptr<PileUpWeight> pileUpWeight_;
  BTagWeight bTagWeight_;

  mutable SmartSelectionMonitor_hzz mon_;
  bool divideFinalHistoByBinWidth_;
  std::vector<std::string> v_jetCat_, tagsR_;
  unsigned tagsR_size_;

  std::vector<double> optim_Cuts1_met_;

  TString fileName_;

  TTreeReaderValue<UInt_t> run_ = {dataset_.Reader(), "run"};
  TTreeReaderValue<Float_t> rho_ = {dataset_.Reader(), "fixedGridRhoFastjetAll"};
  TTreeReaderValue<Int_t> numPVGood_ = {dataset_.Reader(), "PV_npvsGood"};
  TTreeReaderArray<Float_t> muonPt_ = {dataset_.Reader(), "Muon_pt"};
  TTreeReaderArray<Float_t> electronPt_ = {dataset_.Reader(), "Electron_pt"};
  std::unique_ptr<TTreeReaderArray<int>> genPartPdgId_, genPartMotherIndex_;
};

#endif  // HZZ2L2NU_INCLUDE_NRBANALYSIS_H_

