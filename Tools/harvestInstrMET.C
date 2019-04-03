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
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <map>
#include <cstring>

using namespace std;

std::map<TString, double> norma_factor;

void constructNominalInstrMET(TString suffix, TString systType){

  std::cout<< "Harvesting the Instr.MET..."<< std::endl;
  TString base_path = std::string(getenv("HZZ2L2NU_BASE")) + "/";
  TString fileDirectory= base_path+"OUTPUTS/"+suffix+"/MERGED";
  outputPrefixName = "outputPhotonDatadriven_"; //declared in samples.h

  TString files, weights;

  std::vector<MCentry> allSamples;
  TFile* dataFile = new TFile();
  takeHisto_InstrMET(allSamples, &dataFile, fileDirectory);

  //add data photon:
  files += " "+fileDirectory+"/"+outputPrefixName+"Data"+systType+".root";
  weights += " 1";

  //Remove genuine MET from MC samples
  for (MCentry theEntry: allSamples){
    if(theEntry.InstrMETContribution == 0) continue;
    theEntry.sampleFile = new TFile(fileDirectory+"/"+outputPrefixName+theEntry.fileSuffix+systType+".root");
    if(!theEntry.sampleFile->IsOpen()) continue;
    files += " "+fileDirectory+"/"+outputPrefixName+theEntry.fileSuffix+systType+".root";

    TH1F *totEventInBaobab = (TH1F*) (theEntry.sampleFile)->Get("totEventInBaobab_tot");
    float norm = theEntry.InstrMETContribution*instLumi*theEntry.crossSection/totEventInBaobab->Integral(0, totEventInBaobab->GetNbinsX()+1);
    weights += " "+std::to_string(norm);

  }
  system("haddws "+files+" "+weights); //output = result.root

  //After harvesting, set all negative bins histo to 0. Instr.MET doesn't predict negative values
  TFile* file = TFile::Open("result.root");
  TFile *f_output = new TFile(fileDirectory+"/outputHZZ_InstrMET"+systType+".root","RECREATE");
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

  file->Close();
  f_output->Close();
  dataFile->Close();
  system("rm result.root");
  //Final normalization of the Instr.MET
  std::vector<MCentry> allMCsamples;
  TFile* dileptonFile = new TFile();
  MCentry signalEntry;

  bool isDatadriven = true;
  outputPrefixName = "outputHZZ_";
  takeHisto_HZZanalysis(allMCsamples, &dileptonFile, signalEntry, fileDirectory, isDatadriven);

  for (MCentry &theEntry: allMCsamples){
    theEntry.sampleFile = new TFile(fileDirectory+"/"+outputPrefixName+theEntry.fileSuffix+systType+".root");
  }

  TString theHistoName;
  double tot_evt = 0, Instr_evt = 0, all_MC_evt = 0;

  std::vector<TString> jetCat = {"_eq0jets","_geq1jets","_vbf"}; if(systType == "") jetCat.push_back("");
  std::vector<TString> lepCat = {"_ee","_mumu"}; if(systType == "") lepCat.push_back("_ll");

  for(unsigned int i = 0; i < jetCat.size(); i++){
    for(unsigned int j = 0; j < lepCat.size(); j++){

      theHistoName = "MET125_InstrMET_reweighting"+jetCat[i]+lepCat[j];
      TH1F *HZ_data = (TH1F*) dileptonFile->Get(theHistoName);
      tot_evt = HZ_data->Integral(0, HZ_data->GetNbinsX()+1);

      TH1F* MChistos[99]; //Only allow 99 MC processes
      int iteHisto=0;
      THStack *stackMCsamples = new THStack("stackMCsamples",theHistoName);
      for (MCentry theMCentry: allMCsamples){
        MChistos[iteHisto] = (TH1F*) (theMCentry.sampleFile)->Get(theHistoName);
        if (MChistos[iteHisto] == 0) continue;
        TH1F *totEventInBaobab = (TH1F*) (theMCentry.sampleFile)->Get("totEventInBaobab_tot");
        float norm = instLumi*theMCentry.crossSection/totEventInBaobab->Integral(0, totEventInBaobab->GetNbinsX()+1);
        if(theMCentry.crossSection != 0) MChistos[iteHisto]->Scale(norm);
        else{
          Instr_evt = MChistos[iteHisto]->Integral(0, MChistos[iteHisto]->GetNbinsX()+1);
        }
        stackMCsamples->Add(MChistos[iteHisto]);
        delete totEventInBaobab;
        iteHisto++;
      }
      all_MC_evt = ((TH1*)(stackMCsamples->GetStack()->Last()))->Integral(0, ((TH1*)(stackMCsamples->GetStack()->Last()))->GetNbinsX()+1);
      norma_factor[jetCat[i]+lepCat[j]] = (tot_evt - (all_MC_evt - Instr_evt)) / Instr_evt;
      std::cout << "Normalization factor for " << jetCat[i] << lepCat[j] << " = " << norma_factor[jetCat[i]+lepCat[j]] << std::endl;

    }
  }


  TFile *f_final = new TFile(fileDirectory+"/final.root","RECREATE");

  for (MCentry &theEntry: allMCsamples){
    if(theEntry.fileSuffix != "InstrMET") continue;
    std::cout<< "In " << theEntry.fileSuffix << std::endl;
    TIter listInstrMETPlots(theEntry.sampleFile->GetListOfKeys());
    TKey *keyInstrMETPlot;
    while ((keyInstrMETPlot = (TKey*)listInstrMETPlots())) {
      TString name = keyInstrMETPlot->GetTitle();
      TH1F *h_InstrMET = (TH1F*) theEntry.sampleFile->Get(name);
      for(unsigned int i = 0; i < jetCat.size(); i++){
        if(!name.Contains(jetCat[i])) continue;
        if(jetCat[i] == "" && (name.Contains("_eq0jets") || name.Contains("_geq1jets") || name.Contains("_vbf"))) continue;
        for(unsigned int j = 0; j < lepCat.size(); j++){
          if(!name.Contains(lepCat[j])) continue;
          //std::cout << "name = " << name << " and then " << jetCat[i] << lepCat[j] << std::endl;
          h_InstrMET->Scale(norma_factor[jetCat[i]+lepCat[j]]);
        }
      }
      f_final->cd();
      h_InstrMET->Write("");
    }
    delete keyPlot;
  }

  f_final->Close();
  system("mv "+fileDirectory+"/final.root "+fileDirectory+"/outputHZZ_InstrMET"+systType+".root");

  if(systType == "_final"){
    TFile *f_nominal = new TFile(fileDirectory+"/outputHZZ_InstrMET.root", "RECREATE");

    TFile* file_allSyst = TFile::Open(fileDirectory+"/outputHZZ_InstrMET"+systType+".root");
    TIter listPlots_allSyst(file_allSyst->GetListOfKeys());
    TKey *keyPlot_allSyst;
    while ((keyPlot_allSyst = (TKey*)listPlots_allSyst())) {
      if( ((TString) keyPlot_allSyst->GetTitle()).EndsWith("_up") || ((TString) keyPlot_allSyst->GetTitle()).EndsWith("_down") ) continue;
      TH1F *histo = (TH1F*) file_allSyst->Get(keyPlot_allSyst->GetTitle());
      f_nominal->cd();
      histo->Write("");
    }
    delete keyPlot;
    f_nominal->Write("");
  }

}

void systForMCProcess(TString suffix, std::string syst, std::pair<TString, TString> processWithTheSyst){

  std::cout << " Applying " << syst << " on " << processWithTheSyst.first << "..." << std::endl;
  //1. Construct Instr.MET without ZG, WG and W+jets
  TString base_path = std::string(getenv("HZZ2L2NU_BASE")) + "/";
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
    if( ((TString) theEntry.fileSuffix).Contains(processWithTheSyst.first)) continue;
    theEntry.sampleFile = new TFile(fileDirectory+"/"+outputPrefixName+theEntry.fileSuffix+".root");
    if(!theEntry.sampleFile->IsOpen()) continue;
    files += " "+fileDirectory+"/"+outputPrefixName+theEntry.fileSuffix+".root";

    TH1F *totEventInBaobab = (TH1F*) (theEntry.sampleFile)->Get("totEventInBaobab_tot");
    float norm = theEntry.InstrMETContribution*instLumi*theEntry.crossSection/totEventInBaobab->Integral(0, totEventInBaobab->GetNbinsX()+1);
    weights += " "+std::to_string(norm);
    delete theEntry.sampleFile;

  }
  system("haddws "+files+" "+weights); //output = result.root

  //2. Add the shapes up and down. Only for mT for now. Since pdf are only computed for mT this makes sense otherwise too many things are needed to implement.
  std::vector<TString> jetCat = {"_eq0jets","_geq1jets","_vbf"}; //if(systType == "") jetCat.push_back("");
  std::vector<TString> lepCat = {"_ee","_mumu"}; //if(systType == "") lepCat.push_back("_ll");
  TFile* file = TFile::Open("result.root");
  TFile *f_output = new TFile(fileDirectory+"/outputHZZ_InstrMET_"+processWithTheSyst.second+"_"+syst+".root","RECREATE");
  TIter listPlots(file->GetListOfKeys());
  TKey *keyPlot;
  while ((keyPlot = (TKey*)listPlots())) {
    TString name = keyPlot->GetTitle();
    TH1F *histo = (TH1F*) file->Get(name);
    if(!(name.Contains("mT_final"))) continue;
    for (MCentry theEntry: allSamples){
      if(theEntry.InstrMETContribution == 0) continue;
      if(! (((TString) theEntry.fileSuffix).Contains(processWithTheSyst.first))) continue;

      TFile* nominalFile = TFile::Open(fileDirectory+"/"+outputPrefixName+theEntry.fileSuffix+".root");
      theEntry.sampleFile = new TFile(fileDirectory+"/"+outputPrefixName+theEntry.fileSuffix+"_"+syst+".root");
      if(!nominalFile->IsOpen() || !theEntry.sampleFile->IsOpen()) continue;
      TH1F *totEventInBaobab = (TH1F*) (nominalFile)->Get("totEventInBaobab_tot");
      TH1F *h_processWithTheSyst = (TH1F*) (theEntry.sampleFile)->Get(name+"_"+syst);
      float norm = theEntry.InstrMETContribution*instLumi*theEntry.crossSection/totEventInBaobab->Integral(0, totEventInBaobab->GetNbinsX()+1);
      histo->Add(h_processWithTheSyst, norm); //minus sign already in norm.

      delete theEntry.sampleFile;
      delete nominalFile;
    }

    //Remove empty bins and normalize
    for(unsigned int i = 0; i < jetCat.size(); i++){
      if(!name.Contains(jetCat[i])) continue;
      if(jetCat[i] == "" && (name.Contains("_eq0jets") || name.Contains("_geq1jets") || name.Contains("_vbf"))) continue;
      for(unsigned int j = 0; j < lepCat.size(); j++){
        if(!name.Contains(lepCat[j])) continue;
        for(int i = 0; i <= histo->GetNbinsX(); i++){
          if(histo->GetBinContent(i) < 0) histo->SetBinContent(i, 0);
        }
        histo->Scale(norma_factor[jetCat[i]+lepCat[j]]);
      }
    }
    f_output->cd();
    histo->SetTitle(name+"_"+processWithTheSyst.second+"_"+syst);
    histo->SetName(name+"_"+processWithTheSyst.second+"_"+syst);
    histo->Write(name+"_"+processWithTheSyst.second+"_"+syst);
  }
  delete keyPlot;
  delete file;
  delete f_output;
  delete dataFile;
  system("rm result.root");

  //hadd to produce the final InstrMET.
  system("hadd "+fileDirectory+"/finalFile_InstrMET.root "+fileDirectory+"/outputHZZ_InstrMET_final.root "+fileDirectory+"/outputHZZ_InstrMET_"+processWithTheSyst.second+"_"+syst+".root");
  system("mv "+fileDirectory+"/finalFile_InstrMET.root "+fileDirectory+"/outputHZZ_InstrMET_final.root ");
}


std::vector<std::string> exec(const char* cmd) {
  std::array<char, 128> buffer;
  std::vector<std::string> result;
  std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
  if (!pipe) throw std::runtime_error("popen() failed!");
  while (!feof(pipe.get())) {
    if (fgets(buffer.data(), 128, pipe.get()) != nullptr){
      buffer.data()[std::strcspn(buffer.data(), "\r\n")] = 0; // works for LF, CR, CRLF, LFCR, ...  
      result.push_back(buffer.data());
    }
  }
  return result;
}

void constructSyst(TString suffix, std::string syst){

  std::cout << "Constructing Instr.MET systematic: " << syst << std::endl;

  std::vector<std::pair<TString,TString> > processWithTheSyst = {{"WJetsToLNu", "WJets"}, {"WGToLNuG", "WGamma"}, {"ZGTo2", "ZGamma"}};
  for(const auto process: processWithTheSyst) systForMCProcess(suffix, syst, process);

}


void harvestInstrMET(TString suffix, TString systType){
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(kTRUE); //To ensure that all histograms are created with the sum of weights
  TH2::SetDefaultSumw2(kTRUE); //To ensure that all histograms are created with the sum of weights

  //Cleaning
  TString base_path = std::string(getenv("HZZ2L2NU_BASE")) + "/";
  TString fileDirectory= base_path+"OUTPUTS/"+suffix+"/MERGED";
  system("rm "+fileDirectory+"/outputHZZ_InstrMET*root");

  //Nominal syst
  if(systType != "no" && systType != "") systType = "_final";
  else systType = "";
  systSuffixName = systType; //systSuffixName is used by samples.h

  constructNominalInstrMET(suffix, systType);

  //Check syst that were ran on for Instr.MET and compute syst:
  const char *cmd = (std::string("ls ")+fileDirectory.Data()+"/outputPhotonDatadriven_*{up,down}.root | rev | cut -d. -f2,3 | cut -d_ -f1,2 | rev |sort -u").c_str();
  std::vector<std::string> systList = exec(cmd);
  for(const auto s: systList) constructSyst(suffix, s);

}
