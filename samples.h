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

TString outputPrefixName;
TString systSuffixName;

struct MCentry{
  TString nameSample;
  TString legendEntry;
  TString fileSuffix;
  int color;
  float InstrMETContribution; //0 = this sample plays no role for the InstrMET construction; +1 = has to be added for the InstrMET construction; -1 = has to be substracted
  TFile *sampleFile;
  MCentry(TString theNameSample, TString theLegendEntry, TString theFileSuffix, int theColor, float theInstrMETContribution = 0)
    : nameSample(theNameSample), legendEntry(theLegendEntry), fileSuffix(theFileSuffix), color(theColor), InstrMETContribution(theInstrMETContribution)
  {}
  MCentry(){}
  MCentry(const MCentry& other) = default;
};

void takeHisto_HZZanalysis(std::vector<MCentry> & allMCsamples, TFile ** dataFile, std::vector<MCentry> & signalEntry, TString currentDirectory, bool MELA = false, bool isDatadriven = false){
  //MC samples
  allMCsamples.push_back(MCentry("WWTo2L2Nu",  "WW",   "WWTo2L2Nu",            595, 0));
  allMCsamples.push_back(MCentry("WZTo3LNu",   "WZ",   "WZTo3LNu",             590, 0));
  allMCsamples.push_back(MCentry("WZTo2L2Q",   "WZ",   "WZTo2L2Q",             590, 0));
  allMCsamples.push_back(MCentry("ZZTo4L",     "ZZ",   "ZZTo4L",               594, 0));
  allMCsamples.push_back(MCentry("ZZTo2L2Nu",  "ZZ",   "ZZTo2L2Nu",            594, 0));
  allMCsamples.push_back(MCentry("ZZTo2L2Q",   "ZZ",   "ZZTo2L2Q",             594, 0));
  allMCsamples.push_back(MCentry("ggZZ2mu2nu", "ZZ",   "GGToZZTo2Mu2Nu_BSI",   594, 0));
  allMCsamples.push_back(MCentry("ggZZ2e2nu",  "ZZ",   "GGToZZTo2E2Nu_BSI",    594, 0));
  allMCsamples.push_back(MCentry("TTbar",      "Top",  "TT",                   8, 0));
  allMCsamples.push_back(MCentry("TTbar",      "Top",  "TTWJetsToLNu",         8, 0));
  allMCsamples.push_back(MCentry("TTbar",      "Top",  "TTZToLLNuNu_M-10",     8, 0));
  allMCsamples.push_back(MCentry("singleTop",  "Top",  "ST_s-channel",         8, 0));
  allMCsamples.push_back(MCentry("singleTop",  "Top",  "ST_t-channel_top",     8, 0));
  allMCsamples.push_back(MCentry("singleTop",  "Top",  "ST_t-channel_antitop", 8, 0));
  allMCsamples.push_back(MCentry("singleTop",  "Top",  "ST_tW_antitop",        8, 0) );
  allMCsamples.push_back(MCentry("singleTop",  "Top",  "ST_tW_top",            8, 0) );
  allMCsamples.push_back(MCentry("ZZZ",        "VVV",  "ZZZ",                  869, 0));
  allMCsamples.push_back(MCentry("WZZ",        "VVV",  "WZZ",                  869, 0));
  allMCsamples.push_back(MCentry("WWZ",        "VVV",  "WWZ",                  869, 0));
  allMCsamples.push_back(MCentry("WWW",        "VVV",  "WWW",                  869, 0));
  allMCsamples.push_back(MCentry("WJetsToLNu",               "W", "WJetsToLNu_LO",            623, 0));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-100To200",   "W", "WJetsToLNu_HT-100To200",   623, 0));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-200To400",   "W", "WJetsToLNu_HT-200To400",   623, 0));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-400To600",   "W", "WJetsToLNu_HT-400To600",   623, 0));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-600To800",   "W", "WJetsToLNu_HT-600To800",   623, 0));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-800To1200",  "W", "WJetsToLNu_HT-800To1200",  623, 0));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-1200To2500", "W", "WJetsToLNu_HT-1200To2500", 623, 0));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-2500ToInf",  "W", "WJetsToLNu_HT-2500ToInf",  623, 0));
  if(isDatadriven) allMCsamples.push_back(MCentry("InstrMET",        "Instr. #it{p}_{T}^{miss}",    "InstrMET",   833, 0));
  else allMCsamples.push_back(MCentry("DY",        "DY",    "DYJetsToLL_M-50",   833, 0)); //Warning: if you change this line please also change it in WeightsAndDatadriven/InstrMET/macroToComputeClosureTestWeights.C

  //data
  delete *dataFile;
  TFile* tmp = new TFile(currentDirectory+"/"+outputPrefixName+"Data"+systSuffixName+".root"); 
  //TFile* tmp = new TFile(currentDirectory+"/outputHZZ_DoubleMuon-all.root"); 
  //TFile* tmp = new TFile(currentDirectory+"/outputHZZ_DoubleEG-all.root"); 
  *dataFile = tmp;

  //signal
  //Temporary. At the end, we will need all the mass points.
  if(MELA){
    signalEntry.push_back(MCentry("GluGluHToZZTo2L2Nu_M800",  "ggH800W5",  "GluGluHToZZTo2L2Nu_M800_W5_SOnly",  879, 0));
    signalEntry.push_back(MCentry("GluGluHToZZTo2L2Nu_M1500", "ggH1500W5", "GluGluHToZZTo2L2Nu_M1500_W5_SOnly", 602, 0));
  }
  else{
    // signalEntry.push_back(MCentry("GluGluHToZZTo2L2Nu_M800",  "ggH800",  "GluGluHToZZTo2L2Nu_M800",  879, 0));
  }
}

void takeHisto_NRB(std::vector<MCentry> & allMCsamples, TFile ** dataFile, std::vector<MCentry> & signalEntry, TString currentDirectory, bool MELA = false,  bool isDatadriven = false){
  //MC samples
  allMCsamples.push_back(MCentry("WWTo2L2Nu",  "WW",   "WWTo2L2Nu",            595, 0));
  allMCsamples.push_back(MCentry("WZTo3LNu",   "WZ",   "WZTo3LNu",             590, 0));
  allMCsamples.push_back(MCentry("WZTo2L2Q",   "WZ",   "WZTo2L2Q",             590, 0));
  allMCsamples.push_back(MCentry("ZZTo4L",     "ZZ",   "ZZTo4L",               594, 0));
  allMCsamples.push_back(MCentry("ZZTo2L2Nu",  "ZZ",   "ZZTo2L2Nu",            594, 0));
  allMCsamples.push_back(MCentry("ZZTo2L2Q",   "ZZ",   "ZZTo2L2Q",             594, 0));
  allMCsamples.push_back(MCentry("ZZTo2L2Nu", "ZZ#rightarrow Z#tau#tau",    "ZZToTauTau2Nu",         592, 0));
  allMCsamples.push_back(MCentry("ZZTo2L2Q",  "ZZ#rightarrow Z#tau#tau",    "ZZToTauTau2Q",          592, 0));
  allMCsamples.push_back(MCentry("ggZZ2mu2nu", "ZZ",   "GGToZZTo2Mu2Nu_BSI",   594, 0));
  allMCsamples.push_back(MCentry("ggZZ2e2nu",  "ZZ",   "GGToZZTo2E2Nu_BSI",    594, 0));
  allMCsamples.push_back(MCentry("TTbar",      "Top",  "TT",                   8, 0));
  allMCsamples.push_back(MCentry("TTbar",      "Top",  "TTWJetsToLNu",         8, 0));
  allMCsamples.push_back(MCentry("TTbar",      "Top",  "TTZToLLNuNu_M-10",     8, 0));
  allMCsamples.push_back(MCentry("singleTop",  "Top",  "ST_s-channel",         8, 0));
  allMCsamples.push_back(MCentry("singleTop",  "Top",  "ST_t-channel_top",     8, 0));
  allMCsamples.push_back(MCentry("singleTop",  "Top",  "ST_t-channel_antitop", 8, 0));
  allMCsamples.push_back(MCentry("singleTop",  "Top",  "ST_tW_antitop",        8, 0) );
  allMCsamples.push_back(MCentry("singleTop",  "Top",  "ST_tW_top",            8, 0) );
  allMCsamples.push_back(MCentry("ZZZ",        "VVV",  "ZZZ",                  869, 0));
  allMCsamples.push_back(MCentry("WZZ",        "VVV",  "WZZ",                  869, 0));
  allMCsamples.push_back(MCentry("WWZ",        "VVV",  "WWZ",                  869, 0));
  allMCsamples.push_back(MCentry("WWW",        "VVV",  "WWW",                  869, 0));
  allMCsamples.push_back(MCentry("WJetsToLNu",               "W", "WJetsToLNu_LO",            623, -1));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-100To200",   "W", "WJetsToLNu_HT-100To200",   623, -1));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-200To400",   "W", "WJetsToLNu_HT-200To400",   623, -1));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-400To600",   "W", "WJetsToLNu_HT-400To600",   623, -1));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-600To800",   "W", "WJetsToLNu_HT-600To800",   623, -1));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-800To1200",  "W", "WJetsToLNu_HT-800To1200",  623, -1));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-1200To2500", "W", "WJetsToLNu_HT-1200To2500", 623, -1));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-2500ToInf",  "W", "WJetsToLNu_HT-2500ToInf",  623, -1));
  allMCsamples.push_back(MCentry("DY",         "Z#rightarrow ee/#mu#mu",   "DYJetsToLL_M-50",       834, 0));
  allMCsamples.push_back(MCentry("DY",         "Z#rightarrow #tau#tau",    "DYJetsToTauTau_M-50",   833, 0));


  //data
  delete *dataFile;
  TFile* tmp = new TFile(currentDirectory+"/"+outputPrefixName+"Data"+systSuffixName+".root");
  //TFile* tmp = new TFile(currentDirectory+"/outputHZZ_DoubleMuon-all.root"); 
  //TFile* tmp = new TFile(currentDirectory+"/outputHZZ_DoubleEG-all.root"); 
  *dataFile = tmp;
  if(MELA){
    signalEntry.push_back(MCentry("GluGluHToZZTo2L2Nu_M800",  "ggH800W5",  "GluGluHToZZTo2L2Nu_M800_W5_SOnly",  879, 0));
    signalEntry.push_back(MCentry("GluGluHToZZTo2L2Nu_M1500", "ggH1500W5", "GluGluHToZZTo2L2Nu_M1500_W5_SOnly", 602, 0));
  }
  else{
    // signalEntry.push_back(MCentry("GluGluHToZZTo2L2Nu_M800",  "ggH800",  "GluGluHToZZTo2L2Nu_M800",  879, 0));
    // signalEntry.push_back(MCentry("GluGluHToZZTo2L2Nu_M1500", "ggH1500", "GluGluHToZZTo2L2Nu_M1500", 602, 0));
  }
}

void takeHisto_InstrMET(std::vector<MCentry> & allMCsamples, TFile ** dataFile, TString currentDirectory){
  //MC samples
  allMCsamples.push_back(MCentry("TGJets", "Top+#gamma", "TGJets", 8, -1)); 
  allMCsamples.push_back(MCentry("TTGJets", "Top+#gamma", "TTGJets", 8, -1)); 
  allMCsamples.push_back(MCentry("ZGTo2LG", "Z#gamma #rightarrow ll#gamma", "ZGTo2LG", 635, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-100To200", "W#rightarrow l#nu", "WJetsToLNu_HT-100To200", 623, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-1200To2500", "W#rightarrow l#nu", "WJetsToLNu_HT-1200To2500", 623, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-200To400", "W#rightarrow l#nu", "WJetsToLNu_HT-200To400", 623, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-2500ToInf", "W#rightarrow l#nu", "WJetsToLNu_HT-2500ToInf", 623, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-400To600", "W#rightarrow l#nu", "WJetsToLNu_HT-400To600", 623, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-600To800", "W#rightarrow l#nu", "WJetsToLNu_HT-600To800", 623, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-800To1200", "W#rightarrow l#nu", "WJetsToLNu_HT-800To1200", 623, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu", "W#rightarrow l#nu", "WJetsToLNu_LO", 623, -1)); 
  //allMCsamples.push_back(MCentry("ZNuNuGJets_MonoPhoton_PtG-130", "Z#gamma #rightarrow #nu#nu#gamma", "ZNuNuGJets_MonoPhoton_PtG-130", 800, -1)); // This is an other sample than AN2018_003_v3
  //allMCsamples.push_back(MCentry("ZNuNuGJets_MonoPhoton_PtG-40to130", "Z#gamma #rightarrow #nu#nu#gamma", "ZNuNuGJets_MonoPhoton_PtG-40to130", 800, -1)); // This is an other sample than AN2018_003_v3
  allMCsamples.push_back(MCentry("ZGTo2NuG_PtG-130", "Z#gamma #rightarrow #nu#nu#gamma", "ZGTo2NuG_PtG-130", 800, -1)); 
  allMCsamples.push_back(MCentry("ZGTo2NuG", "Z#gamma #rightarrow #nu#nu#gamma", "ZGTo2NuG", 800, -1)); 

  allMCsamples.push_back(MCentry("WGToLNuG", "W#gamma #rightarrow l#nu#gamma", "WGToLNuG", 52, -1)); 

  allMCsamples.push_back(MCentry("QCD_HT1000to1500", "QCD, HT>100", "QCD_HT1000to1500", 21, 0)); 
  allMCsamples.push_back(MCentry("QCD_HT100to200", "QCD, HT>100", "QCD_HT100to200", 21, 0)); //Sample with low stats 
  allMCsamples.push_back(MCentry("QCD_HT1500to2000", "QCD, HT>100", "QCD_HT1500to2000", 21, 0)); 
  allMCsamples.push_back(MCentry("QCD_HT2000toInf", "QCD, HT>100", "QCD_HT2000toInf", 21, 0)); 
  allMCsamples.push_back(MCentry("QCD_HT200to300", "QCD, HT>100", "QCD_HT200to300", 21, 0)); 
  allMCsamples.push_back(MCentry("QCD_HT300to500", "QCD, HT>100", "QCD_HT300to500", 21, 0)); 
  allMCsamples.push_back(MCentry("QCD_HT500to700", "QCD, HT>100", "QCD_HT500to700", 21, 0)); 
  allMCsamples.push_back(MCentry("QCD_HT700to1000", "QCD, HT>100", "QCD_HT700to1000", 21, 0)); 
  //allMCsamples.push_back(MCentry("QCD_HT50to100", "QCD, HT>50", "QCD_HT50to100", 21, 0)); //Only three events... and with big weight! 

  //allMCsamples.push_back(MCentry("QCD_Pt-120to170_EMEnriched", "QCD_EMEnr", "QCD_Pt-120to170_EMEnriched", 24, 0));
  //allMCsamples.push_back(MCentry("QCD_Pt-170to300_EMEnriched", "QCD_EMEnr", "QCD_Pt-170to300_EMEnriched", 24, 0));
  //allMCsamples.push_back(MCentry("QCD_Pt-20to30_EMEnriched", "QCD_EMEnr", "QCD_Pt-20to30_EMEnriched", 24, 0));
  //allMCsamples.push_back(MCentry("QCD_Pt-300toInf_EMEnriched", "QCD_EMEnr", "QCD_Pt-300toInf_EMEnriched", 24, 0));
  //allMCsamples.push_back(MCentry("QCD_Pt-30to50_EMEnriched", "QCD_EMEnr", "QCD_Pt-30to50_EMEnriched", 24, 0)); //Sample with low stats
  //allMCsamples.push_back(MCentry("QCD_Pt-50to80_EMEnriched", "QCD_EMEnr", "QCD_Pt-50to80_EMEnriched", 24, 0));
  //allMCsamples.push_back(MCentry("QCD_Pt-80to120_EMEnriched", "QCD_EMEnr", "QCD_Pt-80to120_EMEnriched", 24, 0));
  //allMCsamples.push_back(MCentry("QCD_Pt-20toInf_MuEnrichedPt15", "QCD_MuEnr", "QCD_Pt-20toInf_MuEnrichedPt15", 25, 0)); //Sample with low stats

  //LO samples - k-factors are applied in the code as a function of photon pT
  allMCsamples.push_back(MCentry("GJets_HT-100To200", "#gamma+jets", "GJets_HT-100To200", 93, 0)); //Warning: if you change these lines please also change it in WeightsAndDatadriven/InstrMET/macroToComputeClosureTestWeights.C
  allMCsamples.push_back(MCentry("GJets_HT-200To400", "#gamma+jets", "GJets_HT-200To400", 93, 0));
  allMCsamples.push_back(MCentry("GJets_HT-400To600", "#gamma+jets", "GJets_HT-400To600", 93, 0)); 
  allMCsamples.push_back(MCentry("GJets_HT-40To100", "#gamma+jets", "GJets_HT-40To100", 93, 0)); 
  allMCsamples.push_back(MCentry("GJets_HT-600ToInf", "#gamma+jets", "GJets_HT-600ToInf", 93, 0)); 

  //NLO samples
  //allMCsamples.push_back(MCentry("GJets_Pt-20To100", "#gamma+jets", "GJets_Pt-20To100", 93, 0)); 
  //allMCsamples.push_back(MCentry("GJets_Pt-100To200", "#gamma+jets", "GJets_Pt-100To200", 93, 0)); 
  //allMCsamples.push_back(MCentry("GJets_Pt-200To500", "#gamma+jets", "GJets_Pt-200To500", 93, 0)); 
  //allMCsamples.push_back(MCentry("GJets_Pt-500To1000", "#gamma+jets", "GJets_Pt-500To1000", 93, 0)); 
  //allMCsamples.push_back(MCentry("GJets_Pt-1000To2000", "#gamma+jets", "GJets_Pt-1000To2000", 93, 0)); 
  //allMCsamples.push_back(MCentry("GJets_Pt-2000To5000", "#gamma+jets", "GJets_Pt-2000To5000", 93, 0)); 


  //data
  delete *dataFile;
  TFile* tmp = new TFile(currentDirectory+"/"+outputPrefixName+"Data"+systSuffixName+".root"); 
  //TFile* tmp = new TFile(currentDirectory+"/outputInstrMET_SinglePhoton-all.root"); 
  *dataFile = tmp;
}
