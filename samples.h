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

float instLumi= 35920.;
TString outputPrefixName;

struct MCentry{
  TString nameSample;
  TString legendEntry;
  TString fileSuffix;
  float crossSection; //set to 0 for datadriven processes
  int color;
  float InstrMETContribution; //0 = this sample plays no role for the InstrMET construction; +1 = has to be added for the InstrMET construction; -1 = has to be substracted
  TFile *sampleFile;
  MCentry(TString theNameSample, TString theLegendEntry, TString theFileSuffix, float theCrossSection, int theColor, float theInstrMETContribution = 0)
    : nameSample(theNameSample), legendEntry(theLegendEntry), fileSuffix(theFileSuffix), crossSection(theCrossSection), color(theColor), InstrMETContribution(theInstrMETContribution)
  {}
};

void takeHisto_HZZanalysis(std::vector<MCentry> & allMCsamples, TFile ** dataFile, TString currentDirectory, bool isDatadriven = false){
  //MC samples
  allMCsamples.push_back(MCentry("WWTo2L2Nu", "WW",    "WWTo2L2Nu",         12.178,    595, 0));
  allMCsamples.push_back(MCentry("WZTo3LNu",  "WZ",    "WZTo3LNu",          4.42965,   590, 0));
  allMCsamples.push_back(MCentry("WZTo2L2Q",  "WZ",    "WZTo2L2Q",          5.595,     590, 0));
  allMCsamples.push_back(MCentry("ZZTo4L",    "ZZ",    "ZZTo4L",            1.256,     594, 0));
  allMCsamples.push_back(MCentry("ZZTo2L2Nu", "ZZ",    "ZZTo2L2Nu",         0.564,     594, 0));
  allMCsamples.push_back(MCentry("ZZTo2L2Q",  "ZZ",    "ZZTo2L2Q",          3.22,      594, 0));
  allMCsamples.push_back(MCentry("TTbar",     "Top",   "TTJets_DiLept",     87.31,       8, 0));
  allMCsamples.push_back(MCentry("W",         "W",     "WJetsToLNu",        61526.7,   623, 0));
  if(isDatadriven) allMCsamples.push_back(MCentry("InstrMET",        "Instr. #it{p}_{T}^{miss}",    "InstrMET",   0,      833, 0));
  else allMCsamples.push_back(MCentry("DY",        "DY",    "DYJetsToLL_M-50",   5765,      833, 0)); //Warning: if you change this line please also change it in WeightsAndDatadriven/InstrMET/macroToComputeClosureTestWeights.C

  //data
  delete *dataFile;
  TFile* tmp = new TFile(currentDirectory+"/"+outputPrefixName+"Data.root"); 
  //TFile* tmp = new TFile(currentDirectory+"/outputHZZ_DoubleMuon-all.root"); 
  //TFile* tmp = new TFile(currentDirectory+"/outputHZZ_DoubleEG-all.root"); 
  *dataFile = tmp;
}

void takeHisto_InstrMET(std::vector<MCentry> & allMCsamples, TFile ** dataFile, TString currentDirectory){
  //MC samples
  allMCsamples.push_back(MCentry("TGJets", "Top+#gamma", "TGJets", 2.967, 8, -1)); 
  allMCsamples.push_back(MCentry("TTGJets", "Top+#gamma", "TTGJets", 3.697, 8, -1)); 
  allMCsamples.push_back(MCentry("ZGTo2LG", "Z#gamma #rightarrow ll#gamma", "ZGTo2LG", 117.864, 635, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-100To200", "W#rightarrow l#nu", "WJetsToLNu_HT-100To200", 1345*1.21, 623, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-1200To2500", "W#rightarrow l#nu", "WJetsToLNu_HT-1200To2500", 1.329*1.21, 623, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-200To400", "W#rightarrow l#nu", "WJetsToLNu_HT-200To400", 359.7*1.21, 623, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-2500ToInf", "W#rightarrow l#nu", "WJetsToLNu_HT-2500ToInf", 0.03216*1.21, 623, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-400To600", "W#rightarrow l#nu", "WJetsToLNu_HT-400To600", 48.91*1.21, 623, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-600To800", "W#rightarrow l#nu", "WJetsToLNu_HT-600To800", 12.05*1.21, 623, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-800To1200", "W#rightarrow l#nu", "WJetsToLNu_HT-800To1200", 5.501*1.21, 623, -1)); 
  allMCsamples.push_back(MCentry("WJetsToLNu", "W#rightarrow l#nu", "WJetsToLNu", 61526.7, 623, -1)); 
  //allMCsamples.push_back(MCentry("ZNuNuGJets_MonoPhoton_PtG-130", "Z#gamma #rightarrow #nu#nu#gamma", "ZNuNuGJets_MonoPhoton_PtG-130", 0.223, 800, -1)); // This is an other sample than AN2018_003_v3
  //allMCsamples.push_back(MCentry("ZNuNuGJets_MonoPhoton_PtG-40to130", "Z#gamma #rightarrow #nu#nu#gamma", "ZNuNuGJets_MonoPhoton_PtG-40to130", 2.816, 800, -1)); // This is an other sample than AN2018_003_v3
  allMCsamples.push_back(MCentry("ZGTo2NuG_PtG-130", "Z#gamma #rightarrow #nu#nu#gamma", "ZGTo2NuG_PtG-130", 0.2768, 800, -1)); 
  allMCsamples.push_back(MCentry("ZGTo2NuG", "Z#gamma #rightarrow #nu#nu#gamma", "ZGTo2NuG", 28.04, 800, -1)); 

  allMCsamples.push_back(MCentry("WGToLNuG", "W#gamma #rightarrow l#nu#gamma", "WGToLNuG", 489, 52, -1)); 

  allMCsamples.push_back(MCentry("QCD_HT1000to1500", "QCD, HT>100", "QCD_HT1000to1500", 1207, 21, 0)); 
  allMCsamples.push_back(MCentry("QCD_HT100to200", "QCD, HT>100", "QCD_HT100to200", 27990000, 21, 0)); //Sample with low stats 
  allMCsamples.push_back(MCentry("QCD_HT1500to2000", "QCD, HT>100", "QCD_HT1500to2000", 119.9, 21, 0)); 
  allMCsamples.push_back(MCentry("QCD_HT2000toInf", "QCD, HT>100", "QCD_HT2000toInf", 25.24, 21, 0)); 
  allMCsamples.push_back(MCentry("QCD_HT200to300", "QCD, HT>100", "QCD_HT200to300", 1712000, 21, 0)); 
  allMCsamples.push_back(MCentry("QCD_HT300to500", "QCD, HT>100", "QCD_HT300to500", 347700, 21, 0)); 
  allMCsamples.push_back(MCentry("QCD_HT500to700", "QCD, HT>100", "QCD_HT500to700", 32100, 21, 0)); 
  allMCsamples.push_back(MCentry("QCD_HT700to1000", "QCD, HT>100", "QCD_HT700to1000", 6831, 21, 0)); 
  //allMCsamples.push_back(MCentry("QCD_HT50to100", "QCD, HT>50", "QCD_HT50to100", 246300000.0, 21, 0)); //Only three events... and with big weight! 

  //allMCsamples.push_back(MCentry("QCD_Pt-120to170_EMEnriched", "QCD_EMEnr", "QCD_Pt-120to170_EMEnriched", 477000*0.132, 24, 0));
  //allMCsamples.push_back(MCentry("QCD_Pt-170to300_EMEnriched", "QCD_EMEnr", "QCD_Pt-170to300_EMEnriched", 114000*0.165, 24, 0));
  //allMCsamples.push_back(MCentry("QCD_Pt-20to30_EMEnriched", "QCD_EMEnr", "QCD_Pt-20to30_EMEnriched", 557600000*0.0096, 24, 0));
  //allMCsamples.push_back(MCentry("QCD_Pt-300toInf_EMEnriched", "QCD_EMEnr", "QCD_Pt-300toInf_EMEnriched", 9000*0.15, 24, 0));
  //allMCsamples.push_back(MCentry("QCD_Pt-30to50_EMEnriched", "QCD_EMEnr", "QCD_Pt-30to50_EMEnriched", 136000000*0.073, 24, 0)); //Sample with low stats
  //allMCsamples.push_back(MCentry("QCD_Pt-50to80_EMEnriched", "QCD_EMEnr", "QCD_Pt-50to80_EMEnriched", 19800000*0.146, 24, 0));
  //allMCsamples.push_back(MCentry("QCD_Pt-80to120_EMEnriched", "QCD_EMEnr", "QCD_Pt-80to120_EMEnriched", 2800000*0.125, 24, 0));
  //allMCsamples.push_back(MCentry("QCD_Pt-20toInf_MuEnrichedPt15", "QCD_MuEnr", "QCD_Pt-20toInf_MuEnrichedPt15", 720648000*0.00042, 25, 0)); //Sample with low stats

  //LO samples - k-factors are applied in the code as a function of photon pT
  allMCsamples.push_back(MCentry("GJets_HT-100To200", "#gamma+jets", "GJets_HT-100To200", 9238, 93, 0)); //Warning: if you change these lines please also change it in WeightsAndDatadriven/InstrMET/macroToComputeClosureTestWeights.C
  allMCsamples.push_back(MCentry("GJets_HT-200To400", "#gamma+jets", "GJets_HT-200To400", 2305, 93, 0));
  allMCsamples.push_back(MCentry("GJets_HT-400To600", "#gamma+jets", "GJets_HT-400To600", 274.4, 93, 0)); 
  allMCsamples.push_back(MCentry("GJets_HT-40To100", "#gamma+jets", "GJets_HT-40To100", 20790, 93, 0)); 
  allMCsamples.push_back(MCentry("GJets_HT-600ToInf", "#gamma+jets", "GJets_HT-600ToInf", 93.46, 93, 0)); 

  //NLO samples
  //allMCsamples.push_back(MCentry("GJets_Pt-20To100", "#gamma+jets", "GJets_Pt-20To100", 137800.0, 93, 0)); 
  //allMCsamples.push_back(MCentry("GJets_Pt-100To200", "#gamma+jets", "GJets_Pt-100To200", 1024.0, 93, 0)); 
  //allMCsamples.push_back(MCentry("GJets_Pt-200To500", "#gamma+jets", "GJets_Pt-200To500", 68.66, 93, 0)); 
  //allMCsamples.push_back(MCentry("GJets_Pt-500To1000", "#gamma+jets", "GJets_Pt-500To1000", 1.014, 93, 0)); 
  //allMCsamples.push_back(MCentry("GJets_Pt-1000To2000", "#gamma+jets", "GJets_Pt-1000To2000", 0.02092, 93, 0)); 
  //allMCsamples.push_back(MCentry("GJets_Pt-2000To5000", "#gamma+jets", "GJets_Pt-2000To5000", 0.00007476, 93, 0)); 


  //data
  delete *dataFile;
  TFile* tmp = new TFile(currentDirectory+"/"+outputPrefixName+"Data.root"); 
  //TFile* tmp = new TFile(currentDirectory+"/outputInstrMET_SinglePhoton-all.root"); 
  *dataFile = tmp;
}
