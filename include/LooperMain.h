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
   double intLumi_;
   int maxEvents_;
   TString outputFile_;
   int isMC_;
   int isPhotonDatadriven_;
   std::string syst_;
   bool keepAllControlPlots_;
   
   
   // Readers to access the data
   TTreeReaderValue<ULong64_t> event = {dataset_.Reader(), "event"};
   TTreeReaderValue<UInt_t> run = {dataset_.Reader(), "run"};
   TTreeReaderValue<UInt_t> luminosityBlock = {dataset_.Reader(), "luminosityBlock"};
   TTreeReaderValue<Int_t> PV_npvsGood = {dataset_.Reader(), "PV_npvsGood"}; //FIXME or PV_npvs, need to check
   TTreeReaderValue<Float_t> fixedGridRhoFastjetAll = {dataset_.Reader(), "fixedGridRhoFastjetAll"};
   //TTreeReaderValue<Float_t> MET_covXX = {dataset_.Reader(), "MET_covXX"}; // Apparently not present in NanoAOD, unfortunately
   //TTreeReaderValue<Float_t> MET_covXY= {dataset_.Reader(), "MET_covXY"};
   //TTreeReaderValue<Float_t> MET_covYY= {dataset_.Reader(), "MET_covYY"};
   //TTreeReaderValue<Float_t> MET_significance = {dataset_.Reader(), "MET_significance"}; // Unfortunately, this is also absent.
   TTreeReaderArray<Float_t> Muon_pt = {dataset_.Reader(), "Muon_pt"};
   TTreeReaderArray<Float_t> Electron_pt = {dataset_.Reader(), "Electron_pt"};
   TTreeReaderArray<Float_t> Photon_pt = {dataset_.Reader(), "Photon_pt"};

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
      intLumi_{options.GetConfig()["luminosity"].as<double>()},
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
