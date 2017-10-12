//-*- c-basic-offset: 4; -*-
#include <iostream>
#include <thread>
#include <TString.h>
#include <TChain.h>
#include "Includes/ArgParser.h"
#include "Includes/HZZ2l2nuLooper.h"
#include "Includes/SmartSelectionMonitor.h"

using namespace std;

int main(int argc, char **argv)
{
  TString bonzaiInputFile = "dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/hbrun/bonzais/testStart/CRAB_PrivateMC/crab_GluGluHToZZTo2L2Nu_M500_13TeV_powheg2_JHUgenV698_pythia8-ZZ2l2vPruner-DMu/170907_115538/0000/Bonzai-GluGluHToZZTo2L2Nu_M500_13TeV_powheg2_JHUgenV698_pythia8-ZZ2l2vPruner-DMu_1.root";

  TString outputFile = "theOutputFile.root";
  int maxEvents = -1;
  int isMC = 0;
  double sampleXsection = -1;


  //--- Parse the arguments -----------------------------------------------------
  if (argc > 1) {
    for (int i = 1; i < argc; ++i) {
        TString currentArg = argv[i];
        //--- possible options ---
        if (currentArg.BeginsWith("bonzaiInputFile=")) {
            getArg(currentArg, bonzaiInputFile);
        }
        else if (currentArg.BeginsWith("histosOutputFile=")) {
            getArg(currentArg, outputFile);
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
    }
  }

  cout << "The Input Direcotory is " << bonzaiInputFile << endl;
  cout << "The output file is " << outputFile << endl;
  cout << "Will run on a max of " << maxEvents << " events" << endl;
  if (isMC) cout << "This file is MC with a cross section of " << sampleXsection <<  endl;
  //HZZ2l2nuLooper myHZZlooper("mySkim.root"); //please add the parameter that you think are needed
  //HZZ2l2nuLooper myHZZlooper("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/hbrun/bonzais/testStart/CRAB_PrivateMC/crab_GluGluHToZZTo2L2Nu_M500_13TeV_powheg2_JHUgenV698_pythia8-ZZ2l2vPruner-DMu/170907_115538/0000/Bonzai-GluGluHToZZTo2L2Nu_M500_13TeV_powheg2_JHUgenV698_pythia8-ZZ2l2vPruner-DMu_1.root"); //please add the parameter that you think are needed
  HZZ2l2nuLooper myHZZlooper(bonzaiInputFile, outputFile, maxEvents, isMC,  sampleXsection);
  myHZZlooper.Loop(); //same here
}
