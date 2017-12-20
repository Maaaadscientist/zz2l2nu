//-*- c-basic-offset: 4; -*-
#include <iostream>
#include <thread>
#include <TString.h>
#include <TChain.h>
#include "Tools/ArgParser.h"
#include "Loopers_Includes/LooperMain.h"

using namespace std;

int main(int argc, char **argv)
{

  TString catalogInputFile = "/storage_mnt/storage/user/hbrun/myEOS/cms/store/user/hbrun/bonzais/Catalogs/Bonzai6octPruner/Bonzais-DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-all-ZZ2l2vPruner-MC_DMu.txt";
  TString outputFile = "theOutputFile.root";
  int maxEvents = -1;
  int isMC = 0;
  double sampleXsection = -1;
  int skipFile = 0;
  int maxFile = 1;
  int doInstrMETAnalysis = 0;
  int doTnPTree = 0;

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
      else if (currentArg.BeginsWith("doInstrMETAnalysis=")) {
        getArg(currentArg, doInstrMETAnalysis);
      }
      else if (currentArg.BeginsWith("doTnPTree=")) {
        getArg(currentArg, doTnPTree);
      }
    }
  }

  cout << "The Input Catalog is " << catalogInputFile << endl;
  cout << "The output file is " << outputFile << endl;
  cout << "Will run on a max of " << maxEvents << " events" << endl;
  if (isMC) cout << "This file is MC with a cross section of " << sampleXsection <<  endl;
  LooperMain myHZZlooper(catalogInputFile, skipFile, maxFile, outputFile, maxEvents, isMC,  sampleXsection);
  if(doInstrMETAnalysis) myHZZlooper.Loop_InstrMET();
  else if(doTnPTree) myHZZlooper.Loop_TnP();
  else myHZZlooper.Loop();
}
