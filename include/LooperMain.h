#ifndef LooperMain_h
#define LooperMain_h

#include <Utils.h>

#include <cstdlib>
#include <iostream>
#include <filesystem>
#include <fstream>

#include <TFile.h>
#include <TString.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TRandom3.h> 

#include <BTagger.h>
#include <Dataset.h>
#include <Logger.h>
#include <Options.h>
#include <SmartSelectionMonitor_hzz.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

namespace fs = std::filesystem;
using std::vector;
using namespace std;

class LooperMain {
public :
   Options const &options_;
   Dataset dataset_;
   int maxEvents_;
   TString outputFile_;
   int isMC_;
   int isPhotonDatadriven_;
   std::string syst_;
   bool keepAllControlPlots_;
   
   
   // Readers to access the data
   TTreeReaderValue<UInt_t> EvtNum = {dataset_.Reader(), "EvtNum"};
   TTreeReaderValue<UInt_t> EvtRunNum = {dataset_.Reader(), "EvtRunNum"};
   TTreeReaderValue<Int_t> EvtLumiNum = {dataset_.Reader(), "EvtLumiNum"};
   TTreeReaderValue<Int_t> EvtVtxCnt = {dataset_.Reader(), "EvtVtxCnt"};
   TTreeReaderValue<Int_t> EvtPuCnt = {dataset_.Reader(), "EvtPuCnt"};
   TTreeReaderValue<Float_t> EvtFastJetRho = {dataset_.Reader(), "EvtFastJetRho"};
   TTreeReaderValue<ULong64_t> TrigMET = {dataset_.Reader(), "TrigMET"};
   TTreeReaderValue<ULong64_t> TrigHltPhot = {dataset_.Reader(), "TrigHltPhot"};
   TTreeReaderValue<ULong64_t> TrigHltMu = {dataset_.Reader(), "TrigHltMu"};
   TTreeReaderValue<ULong64_t> TrigHltDiMu = {dataset_.Reader(), "TrigHltDiMu"};
   TTreeReaderValue<ULong64_t> TrigHltEl = {dataset_.Reader(), "TrigHltEl"};
   TTreeReaderValue<ULong64_t> TrigHltDiEl = {dataset_.Reader(), "TrigHltDiEl"};
   TTreeReaderValue<ULong64_t> TrigHltElMu = {dataset_.Reader(), "TrigHltElMu"};
   TTreeReaderArray<unsigned int> TrigHltPhot_prescale = {dataset_.Reader(), "TrigHltPhot_prescale"};
   TTreeReaderArray<unsigned int> TrigHltMu_prescale = {dataset_.Reader(), "TrigHltMu_prescale"};
   TTreeReaderArray<unsigned int> TrigHltDiMu_prescale = {dataset_.Reader(), "TrigHltDiMu_prescale"};
   TTreeReaderArray<unsigned int> TrigHltEl_prescale = {dataset_.Reader(), "TrigHltEl_prescale"};
   TTreeReaderArray<unsigned int> TrigHltDiEl_prescale = {dataset_.Reader(), "TrigHltDiEl_prescale"};
   TTreeReaderArray<unsigned int> TrigHltElMu_prescale = {dataset_.Reader(), "TrigHltElMu_prescale"};
   TTreeReaderArray<float> METPtType1XY = {dataset_.Reader(), "METPtType1XY"};
   TTreeReaderArray<float> METPhiType1XY = {dataset_.Reader(), "METPhiType1XY"};
   TTreeReaderArray<float> METsigx2 = {dataset_.Reader(), "METsigx2"};
   TTreeReaderArray<float> METsigxy = {dataset_.Reader(), "METsigxy"};
   TTreeReaderArray<float> METsigy2 = {dataset_.Reader(), "METsigy2"};
   TTreeReaderArray<float> METsig = {dataset_.Reader(), "METsig"};
   TTreeReaderArray<float> GLepBarePt = {dataset_.Reader(), "GLepBarePt"};
   TTreeReaderArray<float> GLepBareEta = {dataset_.Reader(), "GLepBareEta"};
   TTreeReaderArray<float> GLepBarePhi = {dataset_.Reader(), "GLepBarePhi"};
   TTreeReaderArray<float> GLepBareE = {dataset_.Reader(), "GLepBareE"};
   TTreeReaderArray<int> GLepBareId = {dataset_.Reader(), "GLepBareId"};
   TTreeReaderArray<int> GLepBareMomId = {dataset_.Reader(), "GLepBareMomId"};
   TTreeReaderArray<float> GPhotPt = {dataset_.Reader(), "GPhotPt"};
   TTreeReaderArray<float> GPhotPrompt = {dataset_.Reader(), "GPhotPrompt"};
   TTreeReaderArray<float> MuPt = {dataset_.Reader(), "MuPt"};
   TTreeReaderArray<float> MuEta = {dataset_.Reader(), "MuEta"};
   TTreeReaderArray<float> MuPhi = {dataset_.Reader(), "MuPhi"};
   TTreeReaderArray<float> MuE = {dataset_.Reader(), "MuE"};
   TTreeReaderArray<unsigned int> MuId = {dataset_.Reader(), "MuId"};
   TTreeReaderArray<unsigned int> MuIdTight = {dataset_.Reader(), "MuIdTight"};
   TTreeReaderArray<float> MuCh = {dataset_.Reader(), "MuCh"};
   TTreeReaderArray<float> MuPfIso = {dataset_.Reader(), "MuPfIso"};
   TTreeReaderArray<unsigned int> MuHltMatch = {dataset_.Reader(), "MuHltMatch"};
   TTreeReaderArray<float> ElPt = {dataset_.Reader(), "ElPt"};
   TTreeReaderArray<float> ElEta = {dataset_.Reader(), "ElEta"};
   TTreeReaderArray<float> ElEtaSc = {dataset_.Reader(), "ElEtaSc"};
   TTreeReaderArray<float> ElPhi = {dataset_.Reader(), "ElPhi"};
   TTreeReaderArray<float> ElE = {dataset_.Reader(), "ElE"};
   TTreeReaderArray<unsigned int> ElId = {dataset_.Reader(), "ElId"};
   TTreeReaderArray<float> ElCh = {dataset_.Reader(), "ElCh"};
   TTreeReaderArray<float> ElPfIsoRho = {dataset_.Reader(), "ElPfIsoRho"};
   TTreeReaderArray<float> PhotPt = {dataset_.Reader(), "PhotPt"};
   TTreeReaderArray<float> PhotEta = {dataset_.Reader(), "PhotEta"};
   TTreeReaderArray<float> PhotPhi = {dataset_.Reader(), "PhotPhi"};
   TTreeReaderArray<float> PhotScEta = {dataset_.Reader(), "PhotScEta"};
   TTreeReaderArray<float> PhotPfIsoChHad = {dataset_.Reader(), "PhotPfIsoChHad"};
   TTreeReaderArray<float> PhotPfIsoNeutralHad = {dataset_.Reader(), "PhotPfIsoNeutralHad"};
   TTreeReaderArray<float> PhotPfIsoPhot = {dataset_.Reader(), "PhotPfIsoPhot"};
   TTreeReaderArray<float> PhotSigmaIetaIeta = {dataset_.Reader(), "PhotSigmaIetaIeta"};
   TTreeReaderArray<float> PhotSigmaIphiIphi = {dataset_.Reader(), "PhotSigmaIphiIphi"};
   TTreeReaderArray<float> PhotHoE = {dataset_.Reader(), "PhotHoE"};
   TTreeReaderArray<float> PhotR9 = {dataset_.Reader(), "PhotR9"};
   TTreeReaderValue<vector<bool>> PhotHasPixelSeed = {dataset_.Reader(), "PhotHasPixelSeed"};

   LooperMain(Options const &options);
   virtual void     Loop();
   virtual void     Loop_InstrMET();
   virtual void     Loop_NRB();

  /**
   * \brief Fills histograms with jets passing b-tagging selection
   */
  void FillBTagEfficiency(std::vector<Jet> const &jets, double weight,
                          SmartSelectionMonitor_hzz &mon, BTagger const &bTagger) const;

private :
   TRandom3 randomGenerator_;
};

#if defined(HZZ2l2nuLooper_cxx) || defined(InstrMETLooper_cxx)
LooperMain::LooperMain(Options const &options)
    : options_(options),
      dataset_{DatasetInfo{options.GetAs<std::string>("catalog"), options},
               options_.GetAs<int>("skip-files"),
               options_.GetAs<int>("max-files")},
      randomGenerator_(options.GetAs<unsigned>("seed")) {

  outputFile_ = options_.GetAs<std::string>("output");
  maxEvents_ = options_.GetAs<long long>("max-events");
  keepAllControlPlots_ = options_.Exists("all-control-plots");
  isPhotonDatadriven_ = options_.Exists("dd-photon");
  syst_ = options_.GetAs<std::string>("syst");

  isMC_ = dataset_.Info().IsSimulation();


  LOG_DEBUG << "Dataset definition file: " << dataset_.Info().DefinitionFile();
  LOG_DEBUG << "The output file is " << outputFile_;
  LOG_DEBUG << "Will run on a max of " << maxEvents_ << " events";

  if (syst_ == "")
    LOG_DEBUG << "Will not use systematic uncertainties";
  else
    LOG_DEBUG << "Will use the systematic " << syst_;

  if (isMC_)
    LOG_DEBUG << "This file is MC with a cross section of " <<
        dataset_.Info().CrossSection();
}
#endif // #if defined(HZZ2l2nuLooper_cxx) || defined(InstrMETLooper_cxx)

#endif
