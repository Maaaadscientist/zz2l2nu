#ifndef HZZ2L2NU_INCLUDE_MAINANALYSIS_H_
#define HZZ2L2NU_INCLUDE_MAINANALYSIS_H_

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

#include <BTagger.h>
#include <BTagWeight.h>
#include <Dataset.h>
#include <ElectronBuilder.h>
#include <EWCorrectionWeight.h>
#include <GenWeight.h>
#include <GenJetBuilder.h>
#include <JetBuilder.h>
#include <KFactorCorrection.h>
#include <MeKinFilter.h>
#include <MelaWeight.h>
#include <MetFilters.h>
#include <MuonBuilder.h>
#include <Options.h>
#include <PhotonBuilder.h>
#include <PhysicsObjects.h>
#include <PileUpWeight.h>
#include <PtMissBuilder.h>
#include <SmartSelectionMonitor_hzz.h>


/**
 * \brief Implements main analysis
 *
 * By default, performs pure MC-based analysis. If flag \c --dd-photon is given,
 * constructs a data-driven prediction for the Drell-Yan background. The results
 * are stored in the output ROOT file in the form of histograms.
 */
class MainAnalysis {
 public:
  MainAnalysis(Options const &options, Dataset &dataset);

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

  /**
   * \brief Fills histograms with jets passing b-tagging selection
   */
  void FillBTagEfficiency(
      std::vector<Jet> const &jets, double weight) const;

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
  std::unique_ptr<KFactorCorrection> kfactorCorrection_;
  MelaWeight melaWeight_;
  BTagWeight bTagWeight_;

  mutable SmartSelectionMonitor_hzz mon_;
  bool divideFinalHistoByBinWidth_;
  std::vector<std::string> v_jetCat_, tagsR_;
  unsigned tagsR_size_;
  std::vector<TH1 *> h_mT_;
  std::vector<int> h_mT_size_;
  TH1 *h_Vtx_, *h_pT_;

  bool isMC_NLO_ZGTo2NuG_inclusive_, isMC_NLO_ZGTo2NuG_Pt130_;

  std::vector<std::vector<std::vector<TH1 *>>> pdfReplicas_;
  std::unique_ptr<std::vector<std::vector<std::vector<std::vector<std::vector<std::pair<double, double>>>>>>> mT_InstrMET_map_;
  std::unique_ptr<std::vector<std::vector<std::vector<std::vector<std::pair<double, double>>>>>> photon_reweighting_;
  std::map<TString, std::map<double, std::pair<double, double>>> nVtxWeight_map_, ptWeight_map_;
  std::map<TString, TH1 *> lineshapeMassWeight_map_;

  TTreeReaderValue<UInt_t> run_ = {dataset_.Reader(), "run"};
  TTreeReaderValue<Float_t> rho_ = {dataset_.Reader(), "fixedGridRhoFastjetAll"};
  TTreeReaderValue<Int_t> numPVGood_ = {dataset_.Reader(), "PV_npvsGood"};
  TTreeReaderArray<Float_t> muonPt_ = {dataset_.Reader(), "Muon_pt"};
  TTreeReaderArray<Float_t> electronPt_ = {dataset_.Reader(), "Electron_pt"};
};

#endif  // HZZ2L2NU_INCLUDE_MAINANALYSIS_H_
