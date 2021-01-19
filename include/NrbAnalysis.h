#ifndef HZZ2L2NU_INCLUDE_NRBANALYSIS_H_
#define HZZ2L2NU_INCLUDE_NRBANALYSIS_H_

#include <map>
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
#include <RunSampler.h>
#include <SmartSelectionMonitor_hzz.h>
#include <TriggerFilter.h>


/// Constructs inputs for data-driven estimation of non-resonant background
class NrbAnalysis : public AnalysisCommon {
 public:
  NrbAnalysis(Options const &options, Dataset &dataset);

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
  enum {eq0jets, eq1jets, geq2jets, jetCat_size};

  void InitializeHistograms();

  Dataset &dataset_;
  std::string outputFile_;
  bool keepAllControlPlots_;
  std::string syst_;

  RunSampler runSampler_;
  TriggerFilter triggerFilter_;

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
