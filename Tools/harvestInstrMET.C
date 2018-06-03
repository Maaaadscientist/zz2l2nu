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
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(kTRUE); //To ensure that all histograms are created with the sum of weights
  TH2::SetDefaultSumw2(kTRUE); //To ensure that all histograms are created with the sum of weights

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

  //After harvesting, set all negative bins histo to 0. Instr.MET doesn't predict negative values
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

  f_output->Close();
  //Final normalization of the Instr.MET
  std::vector<MCentry> allMCsamples;
  TFile* dileptonFile = new TFile();

  bool isDatadriven = true;
  outputPrefixName = "outputHZZ_";
  takeHisto_HZZanalysis(allMCsamples, &dileptonFile, fileDirectory, isDatadriven);

  for (MCentry &theEntry: allMCsamples){
    theEntry.sampleFile = new TFile(fileDirectory+"/"+outputPrefixName+theEntry.fileSuffix+".root");
  }

  TString theHistoName;
  double tot_evt = 0, Instr_evt = 0, all_MC_evt = 0;

  std::vector<TString> jetCat = {"_eq0jets","_geq1jets","_vbf", ""};
  std::vector<TString> lepCat = {"_ee","_mumu", "_ll"};
  std::map<TString, double> factor;

  for(unsigned int i = 0; i < jetCat.size(); i++){
    for(unsigned int j = 0; j < lepCat.size(); j++){

      theHistoName = "MET_InstrMET_reweighting"+jetCat[i]+lepCat[j];
      TH1F *HZ_data = (TH1F*) dileptonFile->Get(theHistoName);
      tot_evt = HZ_data->Integral();

      TH1F* MChistos[99]; //Only allow 99 MC processes
      int iteHisto=0;
      THStack *stackMCsamples = new THStack("stackMCsamples",theHistoName);
      for (MCentry theMCentry: allMCsamples){
        MChistos[iteHisto] = (TH1F*) (theMCentry.sampleFile)->Get(theHistoName);
        if (MChistos[iteHisto] == 0) continue;
        TH1F *totEventInBaobab = (TH1F*) (theMCentry.sampleFile)->Get("totEventInBaobab_tot");
        float norm = instLumi*theMCentry.crossSection/totEventInBaobab->Integral();
        if(theMCentry.crossSection != 0) MChistos[iteHisto]->Scale(norm);
        else{
          Instr_evt = MChistos[iteHisto]->Integral();
        }
        stackMCsamples->Add(MChistos[iteHisto]);
        delete totEventInBaobab;
        iteHisto++;
      }
      all_MC_evt = ((TH1*)(stackMCsamples->GetStack()->Last()))->Integral();
      factor[jetCat[i]+lepCat[j]] = (tot_evt - (all_MC_evt - Instr_evt)) / Instr_evt;
      std::cout << "factor = " << factor[jetCat[i]+lepCat[j]] << std::endl;

    }
  }


  TFile *f_final = new TFile(fileDirectory+"/final.root","RECREATE");

  for (MCentry &theEntry: allMCsamples){
    if(theEntry.fileSuffix != "InstrMET") continue;
    TIter listInstrMETPlots(theEntry.sampleFile->GetListOfKeys());
    TKey *keyInstrMETPlot;
    while ((keyInstrMETPlot = (TKey*)listInstrMETPlots())) {
      TString name = keyInstrMETPlot->GetTitle();
      TH1F *h_InstrMET = (TH1F*) theEntry.sampleFile->Get(name);
      for(unsigned int i = 0; i < jetCat.size(); i++){
        if(!name.Contains(jetCat[i])) continue;
        for(unsigned int j = 0; j < lepCat.size(); j++){
          if(!name.Contains(lepCat[j])) continue;
          h_InstrMET->Scale(factor[jetCat[i]+lepCat[j]]);
        }
      }
      f_final->cd();
      h_InstrMET->Write("");
    }
    delete keyPlot;
  }

system("mv "+fileDirectory+"/final.root "+fileDirectory+"/outputHZZ_InstrMET.root");

}


