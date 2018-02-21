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

#define VERBOSE false

float instLumi;

struct MCentry{
  TString nameSample;
  TString legendEntry;
  TString fileSuffix;
  float crossSection;
  int color;
  TFile *sampleFile;
  MCentry(TString theNameSample, TString theLegendEntry, TString theFileSuffix, float theCrossSection, int theColor)
    : nameSample(theNameSample), legendEntry(theLegendEntry), fileSuffix(theFileSuffix), crossSection(theCrossSection), color(theColor)
  {}
};

void takeHisto_HZZanalysis(std::vector<MCentry> & allMCsamples, TFile ** dataFile, TString currentDirectory){
  //MC samples
  allMCsamples.push_back(MCentry("WWTo2L2Nu", "WW",    "WWTo2L2Nu",         12.178,    595));
  allMCsamples.push_back(MCentry("WZTo3LNu",  "WZ",    "WZTo3LNu",          4.42965,   590));
  allMCsamples.push_back(MCentry("WZTo2L2Q",  "WZ",    "WZTo2L2Q",          5.595,     590));
  allMCsamples.push_back(MCentry("ZZTo4L",    "ZZ",    "ZZTo4L",            1.256,     594));
  allMCsamples.push_back(MCentry("ZZTo2L2Nu", "ZZ",    "ZZTo2L2Nu",         0.564,     594));
  allMCsamples.push_back(MCentry("ZZTo2L2Q",  "ZZ",    "ZZTo2L2Q",          3.22,      594));
  allMCsamples.push_back(MCentry("TTbar",     "Top",   "TTJets_DiLept",     87.31,     8));
  allMCsamples.push_back(MCentry("W",         "W",     "WJetsToLNu",        61526.7,   623));
  allMCsamples.push_back(MCentry("DY",        "DY",    "DYJetsToLL_M-50",   5765,      833));

  //data
  delete *dataFile;
  TFile* tmp = new TFile(currentDirectory+"/output_Data.root"); 
  //TFile* tmp = new TFile(currentDirectory+"/output_DoubleMuon-all.root"); 
  //TFile* tmp = new TFile(currentDirectory+"/output_DoubleEG-all.root"); 
  *dataFile = tmp;
}

void takeHisto_InstrMET(std::vector<MCentry> & allMCsamples, TFile ** dataFile, TString currentDirectory){
  //MC samples
  allMCsamples.push_back(MCentry("TGJets", "Top+#gamma", "TGJets", 2.967, 8));
  allMCsamples.push_back(MCentry("TTGJets", "Top+#gamma", "TTGJets", 3.697, 8));
  allMCsamples.push_back(MCentry("ZGTo2LG", "Z#gamma #rightarrow ll#gamma", "ZGTo2LG", 117.864, 635));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-100To200", "W#rightarrow l#nu", "WJetsToLNu_HT-100To200", 1345, 623));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-1200To2500", "W#rightarrow l#nu", "WJetsToLNu_HT-1200To2500", 1.329, 623));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-200To400", "W#rightarrow l#nu", "WJetsToLNu_HT-200To400", 359.7, 623));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-2500ToInf", "W#rightarrow l#nu", "WJetsToLNu_HT-2500ToInf", 0.03216, 623));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-400To600", "W#rightarrow l#nu", "WJetsToLNu_HT-400To600", 48.91, 623));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-600To800", "W#rightarrow l#nu", "WJetsToLNu_HT-600To800", 12.05, 623));
  allMCsamples.push_back(MCentry("WJetsToLNu_HT-800To1200", "W#rightarrow l#nu", "WJetsToLNu_HT-800To1200", 5.501, 623));
  allMCsamples.push_back(MCentry("WJetsToLNu", "W#rightarrow l#nu", "WJetsToLNu", 61526.7, 623));
  allMCsamples.push_back(MCentry("ZNuNuGJets_MonoPhoton_PtG-130", "Z#gamma #rightarrow #nu#nu#gamma", "ZNuNuGJets_MonoPhoton_PtG-130", 0.223, 800));
  allMCsamples.push_back(MCentry("ZNuNuGJets_MonoPhoton_PtG-40to130", "Z#gamma #rightarrow #nu#nu#gamma", "ZNuNuGJets_MonoPhoton_PtG-40to130", 2.816, 800));
  allMCsamples.push_back(MCentry("WGToLNuG", "W#gamma #rightarrow l#nu#gamma", "WGToLNuG", 489, 52));
  allMCsamples.push_back(MCentry("QCD_HT1000to1500", "QCD, HT>100", "QCD_HT1000to1500", 1207, 21));
  allMCsamples.push_back(MCentry("QCD_HT100to200", "QCD, HT>100", "QCD_HT100to200", 27990000, 21)); //Sample with low stats
  allMCsamples.push_back(MCentry("QCD_HT1500to2000", "QCD, HT>100", "QCD_HT1500to2000", 119.9, 21));
  allMCsamples.push_back(MCentry("QCD_HT2000toInf", "QCD, HT>100", "QCD_HT2000toInf", 25.24, 21));
  allMCsamples.push_back(MCentry("QCD_HT200to300", "QCD, HT>100", "QCD_HT200to300", 1712000, 21));
  allMCsamples.push_back(MCentry("QCD_HT300to500", "QCD, HT>100", "QCD_HT300to500", 347700, 21));
  allMCsamples.push_back(MCentry("QCD_HT500to700", "QCD, HT>100", "QCD_HT500to700", 32100, 21));
  allMCsamples.push_back(MCentry("QCD_HT700to1000", "QCD, HT>100", "QCD_HT700to1000", 6831, 21));
  allMCsamples.push_back(MCentry("QCD_Pt-120to170_EMEnriched", "QCD_EMEnr", "QCD_Pt-120to170_EMEnriched", 477000*0.132, 24));
  allMCsamples.push_back(MCentry("QCD_Pt-170to300_EMEnriched", "QCD_EMEnr", "QCD_Pt-170to300_EMEnriched", 114000*0.165, 24));
  allMCsamples.push_back(MCentry("QCD_Pt-20to30_EMEnriched", "QCD_EMEnr", "QCD_Pt-20to30_EMEnriched", 557600000*0.0096, 24));
  allMCsamples.push_back(MCentry("QCD_Pt-20toInf_MuEnrichedPt15", "QCD_EMEnr", "QCD_Pt-20toInf_MuEnrichedPt15", 720648000*0.00042, 24)); //Sample with low stats
  allMCsamples.push_back(MCentry("QCD_Pt-300toInf_EMEnriched", "QCD_EMEnr", "QCD_Pt-300toInf_EMEnriched", 9000*0.15, 24));
  allMCsamples.push_back(MCentry("QCD_Pt-30to50_EMEnriched", "QCD_EMEnr", "QCD_Pt-30to50_EMEnriched", 136000000*0.073, 24)); //Sample with low stats
  allMCsamples.push_back(MCentry("QCD_Pt-50to80_EMEnriched", "QCD_EMEnr", "QCD_Pt-50to80_EMEnriched", 19800000*0.146, 24));
  allMCsamples.push_back(MCentry("QCD_Pt-80to120_EMEnriched", "QCD_EMEnr", "QCD_Pt-80to120_EMEnriched", 2800000*0.125, 24));
  allMCsamples.push_back(MCentry("GJets_HT-100To200", "#gamma+jets", "GJets_HT-100To200", 9226.0, 390));
  allMCsamples.push_back(MCentry("GJets_HT-200To400", "#gamma+jets", "GJets_HT-200To400", 2300.0, 390));
  allMCsamples.push_back(MCentry("GJets_HT-400To600", "#gamma+jets", "GJets_HT-400To600", 277.4, 390));
  allMCsamples.push_back(MCentry("GJets_HT-40To100", "#gamma+jets", "GJets_HT-40To100", 20730.0, 390));
  allMCsamples.push_back(MCentry("GJets_HT-600ToInf", "#gamma+jets", "GJets_HT-600ToInf", 93.38, 390));

  //data
  delete *dataFile;
  TFile* tmp = new TFile(currentDirectory+"/output_Data.root"); 
  //TFile* tmp = new TFile(currentDirectory+"/output_SinglePhoton-all.root"); 
  *dataFile = tmp;

}

void doMetFilterEfficiencyPlots(TH1F* MZ_data, THStack * stackMCsamples){
  MZ_data->Scale(1.0/MZ_data->GetBinContent(MZ_data->GetSize()-2));

  TList *histKeys = stackMCsamples->GetHists();
  TIter next(histKeys);
  TObject* object = 0;
  TObjArray * list = stackMCsamples->GetStack();
  double totalEvent = ((TH1F*)list->Last())->GetBinContent(((TH1F*)list->Last())->GetSize()-2);
  while ((object = next()))
  {
    ((TH1F*)object)->Scale(1.0/totalEvent);
  }
  stackMCsamples->Modified(); 

}

void drawTheHisto(TFile *dataFile, std::vector<MCentry> allMCsamples, TString theHistoName, TString suffix, TString typeObject){
  gROOT->SetBatch();
  if(typeObject.Contains("TH1")) typeObject = "TH1";
  else if(typeObject.Contains("TH2")) typeObject = "TH2";
  else cout << "/!\\ WARNING /!\\ You have an histogram that is not a TH1 or a TH2 and I don't know how to draw it... so I don't." << endl;

  if(VERBOSE) cout<< "In draw the histo for "<<theHistoName<<endl;
  TH1F *MZ_data = (TH1F*) dataFile->Get(theHistoName);
  TH1F *totEventInBaobab_tot_data = (TH1F*) dataFile->Get("totEventInBaobab_tot");
  if(VERBOSE) cout << "the tot events data =" << totEventInBaobab_tot_data->Integral() << endl;


  TCanvas *c0 = new TCanvas("c0","canvas",600,800);
  TPad *pad =new TPad("up","up",0,0.25,1,1);
  pad->SetNumber(1);
  pad->SetGridx();
  pad->SetGridy();
  pad->SetTicky();
  pad->SetBottomMargin(0.006);
  pad->Draw();

  c0->cd();
  TPad *pad2 =new TPad("down","down",0,0,1,0.25);
  pad2->SetNumber(2);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->Draw();

  c0->cd();
  pad->cd();
  TLegend *t = new TLegend(0.69,0.66,0.99,0.99);
  t->SetLineColor(0);
  t->AddEntry(MZ_data, "Data", "l");

  if(typeObject== "TH1"){
    MZ_data->SetMarkerColor(kBlack);
    MZ_data->SetLineColor(kBlack);
  }

  TH1F* MChistos[99]; //Only allow 99 MC processes
  int iteHisto=0;
  TString lastLegend = "";
  THStack *stackMCsamples = new THStack("stackMCsamples",theHistoName);
  for (MCentry theMCentry: allMCsamples){
    if(VERBOSE) cout << "doing " << theMCentry.nameSample << endl;
    MChistos[iteHisto] = (TH1F*) (theMCentry.sampleFile)->Get(theHistoName);
    if (MChistos[iteHisto] == 0) continue;
    if(VERBOSE) cout << "found" << endl;
    TH1F *totEventInBaobab = (TH1F*) (theMCentry.sampleFile)->Get("totEventInBaobab_tot");
    float norm = instLumi*theMCentry.crossSection/totEventInBaobab->Integral();
    if(VERBOSE) cout << "scale is " << norm << endl;
    MChistos[iteHisto]->Scale(norm);
    if(typeObject== "TH1") MChistos[iteHisto]->SetLineColor(theMCentry.color);
    else if(typeObject== "TH2") MChistos[iteHisto]->SetLineColor(kBlack);
    MChistos[iteHisto]->SetFillColor(theMCentry.color);
    stackMCsamples->Add(MChistos[iteHisto]);
    if (lastLegend !=theMCentry.legendEntry){
      t->AddEntry(MChistos[iteHisto], theMCentry.legendEntry, "F");
      lastLegend = theMCentry.legendEntry;
    }
    delete totEventInBaobab;
    iteHisto++;
  }

  if(theHistoName == "metFilters_tot") doMetFilterEfficiencyPlots(MZ_data, stackMCsamples);

  c0->cd();
  pad->cd();

  if(typeObject== "TH1"){
    MZ_data->Draw("E1:same");
    stackMCsamples->Draw("HIST:same");
    MZ_data->Draw("E1:same");
  }
  else if(typeObject== "TH2"){
    stackMCsamples->Draw("");
    MZ_data->Draw("LEGO:same");
  }
  t->Draw();

  pad2->cd();
  TH1F *ratio = (TH1F*) MZ_data->Clone("ratio");
  ratio->Sumw2();
  ratio->Divide(MZ_data, ((TH1F*)stackMCsamples->GetStack()->Last()), 1,1);
  if(typeObject== "TH1"){
    ratio->SetMaximum(1.3);
    ratio->SetMinimum(0.7);
  }
  else if(typeObject== "TH2"){
    ratio->SetMaximum(1.5);
    ratio->SetMinimum(0.5);

  }
  ratio->SetTitle("");
  ratio->GetYaxis()->SetTitle("");
  ratio->GetYaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetLabelSize(0.1);
  ratio->GetXaxis()->SetLabelOffset(0.02);
  if(typeObject== "TH1") ratio->Draw("E1");
  else if(typeObject== "TH2") ratio->Draw("LEGO");
  TString outputDir = "OUTPUTS/"+suffix+"/PLOTS/";
  c0->Print(outputDir+theHistoName+".png");
  c0->Print(outputDir+theHistoName+".root");
  pad->cd();
  if(typeObject== "TH1") pad->cd()->SetLogy();
  else if(typeObject== "TH2") pad->cd()->SetLogz();
  MZ_data->SetMinimum(0.0001);
  MZ_data->Draw("E1:same");
  c0->Print(outputDir+theHistoName+"_log.png");
  delete c0;
  //  for (int i=0 ; i < iteHisto ; i++){  delete MChistos[i];}
  //  delete stackMCsamples;
}


void dataMCcomparison(TString analysisType, TString suffix){
  TString currentDirectory="OUTPUTS/"+suffix+"/MERGED";
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(kTRUE); //To ensure that all histograms are created with the sum of weights
  TH2::SetDefaultSumw2(kTRUE); //To ensure that all histograms are created with the sum of weights
  gErrorIgnoreLevel = kWarning;//kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal;
  if(VERBOSE) gErrorIgnoreLevel = kPrint;
  //gStyle->SetOptTitle(0);

  //This will be the order to draw, so smallest XS should be first line
  std::vector<MCentry> allMCsamples;
  TFile* dataFile = new TFile();

  if(analysisType == "HZZanalysis"){
    takeHisto_HZZanalysis(allMCsamples, &dataFile, currentDirectory);
    instLumi=16916.98*0.88;//FIXME 16916.16 = int. lumi of the double muon sample and the 0.88 is a SF put by hand;
  }  
  else if(analysisType == "InstrMET"){
    takeHisto_InstrMET(allMCsamples, &dataFile, currentDirectory);
    instLumi= 35866.932; //FIXME This value comes from the old code, from the baobabs it is: 16916.39+20092.7;//Lumi of runs GH + BtoF for singlePhoton samples;
  }

  for (MCentry &theEntry: allMCsamples){
    theEntry.sampleFile = new TFile(currentDirectory+"/output_"+theEntry.fileSuffix+".root");
  }

  TIter listPlots(dataFile->GetListOfKeys());
  TKey *keyPlot;
  while ((keyPlot = (TKey*)listPlots())) {
    TString typeObject = keyPlot->GetClassName();
    TString nomObject = keyPlot->GetTitle();
    if (nomObject.Contains("totEventInBaobab")) continue;
    if(VERBOSE) cout << "Type:" << typeObject << " and title:" << nomObject << endl;
    drawTheHisto(dataFile, allMCsamples, nomObject, suffix, typeObject);
  }

  //drawTheHisto(dataFile, allMCsamples, "M_Z_tot_mumu", suffix);
  //  drawTheHisto(dataFile, allMCsamples, "eventflow_tot", suffix);
  //  drawTheHisto(dataFile, allMCsamples, "MET_beforeMETcut", suffix);
  //drawTheHisto(dataFile, DYfile, "M_Z_tot_mumu_eq0jets", suffix);
  //drawTheHisto(dataFile, DYfile, "jetCategory_tot");




}
