//-*- c-basic-offset: 4; -*-
#include <iostream>
#include <thread>
#include <TString.h>
#include <TChain.h>

#include <ArgParser.h>
#include <LooperMain.h>

using namespace std;

int main(int argc, char **argv)
{

  //TString catalogInputFile = "/storage_mnt/storage/user/hbrun/myEOS/cms/store/user/hbrun/bonzais/Catalogs/Bonzai6octPruner/Bonzais-DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-all-ZZ2l2vPruner-MC_DMu.txt";
  TString catalogInputFile = "/user/npostiau/event_files/MC_ewk/Bonzais-catalog_test_ZZTo2L2Nu-ZZ2l2vPruner.txt";
  TString outputFile = "outputFile.root";
  int maxEvents = -1;
  int isMC = 1;
  double sampleXsection = -1;
  int skipFile = 0;
  int maxFile = 1;
  int isPhotonDatadriven = 0;
  int doInstrMETAnalysis = 0;
  int doTnPTree = 0;
  int doNRBAnalysis =0;
  bool keepAllControlPlots = true;
  TString syst = "";
  int seed = 0;

  //--- Parse the arguments -----------------------------------------------------
  if (argc > 1) {
    for (int i = 1; i < argc; ++i) {
      TString currentArg = argv[i];
      //--- possible options ---
      if (currentArg.BeginsWith("catalogInputFile=")) {
        getArg(currentArg, catalogInputFile);
      }
      else if (currentArg.BeginsWith("histosOutputFile=")) {
        getArg(currentArg, outputFile);
      }
      else if (currentArg.BeginsWith("skip-files=")) {
        getArg(currentArg, skipFile);
      }
      else if (currentArg.BeginsWith("max-files=")) {
        getArg(currentArg, maxFile);
      }
      else if (currentArg.BeginsWith("isMC=")) {
        getArg(currentArg, isMC);
      }
      else if (currentArg.BeginsWith("crossSection=")) {
        getArg(currentArg, sampleXsection);
      }
      else if (currentArg.BeginsWith("maxEvents=")) {
        getArg(currentArg, maxEvents);
      }
      else if (currentArg.BeginsWith("isPhotonDatadriven=")) { //launch the HZZanalyis but with photon data specific options. This is not useful for the Photon CR.
        getArg(currentArg, isPhotonDatadriven);
      }
      else if (currentArg.BeginsWith("doInstrMETAnalysis=")) { //launch the Photon CR analysis
        getArg(currentArg, doInstrMETAnalysis);
      }
      else if (currentArg.BeginsWith("doTnPTree=")) {
        getArg(currentArg, doTnPTree);
      }
      else if (currentArg.BeginsWith("doNRBAnalysis=")) {
        getArg(currentArg, doNRBAnalysis);
      }
      else if (currentArg.BeginsWith("syst=")) {
        getArg(currentArg, syst);
      }
      else if (currentArg.BeginsWith("keepAllControlPlots=")) {
        getArg(currentArg, keepAllControlPlots);
      }
      else if (currentArg.BeginsWith("seed")) {
        getArg(currentArg, seed);
      }
    }
  }

  cout << "The Input Catalog is " << catalogInputFile << endl;
  cout << "The output file is " << outputFile << endl;
  cout << "Will run on a max of " << maxEvents << " events" << endl;
  if(syst=="") cout << "Will not use systematic uncertainties" << endl;
  else cout << "Will use the systematic " << syst << endl;
  if (isMC) cout << "This file is MC with a cross section of " << sampleXsection <<  endl;
  LooperMain myHZZlooper(
    catalogInputFile, skipFile, maxFile, outputFile, maxEvents, isMC,
    sampleXsection, syst, keepAllControlPlots, isPhotonDatadriven, seed);
  if(doInstrMETAnalysis) myHZZlooper.Loop_InstrMET();
  else if(doTnPTree) myHZZlooper.Loop_TnP();
  else if(doNRBAnalysis) myHZZlooper.Loop_NRB();
  else myHZZlooper.Loop();
}
