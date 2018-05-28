#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <iostream>
#include "TROOT.h"
#include <TGaxis.h>
#include <TLine.h>
#include <algorithm>
#include <TStyle.h>
#include <TKey.h>
#include "../samples.h"

void harvestInstrMET(TString suffix){
  std::cout<< "Harvesting the Instr.MET..."<< std::endl;
  TString base_path = std::string(getenv("CMSSW_BASE")) + "/src/shears/HZZ2l2nu/";
  TString fileDirectory= base_path+"OUTPUTS/"+suffix+"/MERGED";
  outputPrefixName = "outputPhotonDatadriven_"; //declared in samples.h

  TString files, weights;

  std::vector<MCentry> allSamples;
  TFile* dataFile = new TFile();
  takeHisto_InstrMET(allSamples, &dataFile, fileDirectory);

  //add data photon:
  files += " "+fileDirectory+"/"+outputPrefixName+"Data.root";
  weights += " 1";

  //Remove genuine MET from MC samples
  for (MCentry theEntry: allSamples){
    if(theEntry.InstrMETContribution == 0) continue;
    theEntry.sampleFile = new TFile(fileDirectory+"/"+outputPrefixName+theEntry.fileSuffix+".root");
    if(!theEntry.sampleFile->IsOpen()) continue;
    files += " "+fileDirectory+"/"+outputPrefixName+theEntry.fileSuffix+".root";

    TH1F *totEventInBaobab = (TH1F*) (theEntry.sampleFile)->Get("totEventInBaobab_tot");
    float norm = theEntry.InstrMETContribution*instLumi*theEntry.crossSection/totEventInBaobab->Integral();
    weights += " "+std::to_string(norm);

  }

  int result = system(base_path+"Tools/haddws/haddws "+files+" "+weights);
  //system(" mv result.root "+fileDirectory+"/outputHZZ_InstrMET.root");

  //After harvesting, set all negative bins histo to 0. Instr.MET doesn't predict negative values
  //TFile* file = TFile::Open(fileDirectory+"/outputHZZ_InstrMET.root", "UPDATE");
  TFile* file = TFile::Open("result.root");
  TFile *f_output = new TFile(fileDirectory+"/outputHZZ_InstrMET.root","RECREATE");
  TIter listPlots(file->GetListOfKeys());
  TKey *keyPlot;
  while ((keyPlot = (TKey*)listPlots())) {
    TH1F *histo = (TH1F*) file->Get(keyPlot->GetTitle());
    for(int i = 0; i <= histo->GetNbinsX(); i++){
      if(histo->GetBinContent(i) < 0) histo->SetBinContent(i, 0);
    }
    f_output->cd();
    histo->Write("");
  }
  delete keyPlot;

}


