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
#include "../Tools/CMS_lumi.C"

#define VERBOSE true

TString base_path = std::string(getenv("CMSSW_BASE")) + "/src/shears/HZZ2l2nu/";
float instLumi;
TString outputPrefixName;

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

void takeHisto_LO(std::vector<MCentry> & allMCsamples){
  //GJets_HT
  allMCsamples.push_back(MCentry("GJets_HT-100To200", "#gamma+jets LO", "GJets_HT-100To200", 9238, 21)); 
  allMCsamples.push_back(MCentry("GJets_HT-200To400", "#gamma+jets LO", "GJets_HT-200To400", 2305, 21));
  allMCsamples.push_back(MCentry("GJets_HT-400To600", "#gamma+jets LO", "GJets_HT-400To600", 274.4, 21)); 
  allMCsamples.push_back(MCentry("GJets_HT-40To100", "#gamma+jets LO", "GJets_HT-40To100", 20790, 21)); 
  allMCsamples.push_back(MCentry("GJets_HT-600ToInf", "#gamma+jets LO", "GJets_HT-600ToInf", 93.46, 21)); 
}

void takeHisto_NLO(std::vector<MCentry> & allMCsamples){
  //GJets_Pt
  allMCsamples.push_back(MCentry("GJets_Pt-20To100", "#gamma+jets NLO", "GJets_Pt-20To100", 137800.0, 93)); 
  allMCsamples.push_back(MCentry("GJets_Pt-100To200", "#gamma+jets NLO", "GJets_Pt-100To200", 1024.0, 93)); 
  allMCsamples.push_back(MCentry("GJets_Pt-200To500", "#gamma+jets NLO", "GJets_Pt-200To500", 68.66, 93)); 
  allMCsamples.push_back(MCentry("GJets_Pt-500To1000", "#gamma+jets NLO", "GJets_Pt-500To1000", 1.014, 93)); 
  allMCsamples.push_back(MCentry("GJets_Pt-1000To2000", "#gamma+jets NLO", "GJets_Pt-1000To2000", 0.02092, 93)); 
  allMCsamples.push_back(MCentry("GJets_Pt-2000To5000", "#gamma+jets NLO", "GJets_Pt-2000To5000", 0.00007476, 93)); 
}


void make_axis(TAxis* & xaxis, TAxis* & yaxis, int fontType, int pixelFontSize){
  yaxis->SetLabelFont(fontType);
  yaxis->SetLabelSize(pixelFontSize);
  //yaxis->SetNdivisions(505);
  yaxis->SetTitleOffset(2.1); //1.4
  yaxis->SetTitleFont(fontType);
  yaxis->SetTitleSize(pixelFontSize);

}

void writeCondition(ofstream& outputFile, TString bin_name, double binLow, double binUp, int bin, int lastBin = -1, int numberOfTabulations =0)
{
  TString condition;
  if(bin == 1) condition = "if";
  else condition = "else if";
  TString tabs = "";
  for(int i = 0; i < numberOfTabulations; i++) tabs += "\t";


  if(bin != lastBin) outputFile << tabs << condition << "(" << bin_name <<" >= " << binLow << " && " << bin_name << " < " << binUp << ")";
  else outputFile << "\telse";
}

void writeEfficiency(ofstream& outputFile, double weight)
{
  outputFile << " weight *=" << weight << ";\n";
}

void drawTheHisto(std::vector<MCentry> LO_MCsamples, std::vector<MCentry> NLO_MCsamples, TString theHistoName, TString suffix){

  if(VERBOSE) cout<< "In draw the histo for "<<theHistoName<<endl;

  //////////////////////////////////////////
  /////////// Configuration part ///////////
  //////////////////////////////////////////

  int fontType = 43; //precision 3 font, so label size is expressed in pixel now
  int pixelFontSize = 22;

  int W = 600;
  int H = 800;
  int W_ref = 600;
  int H_ref = 800;

  //references for T, B, L, R
  float T = 0.09*H_ref;
  float B = 0.30*H_ref;
  float L = 0.14*W_ref;
  float R = 0.04*W_ref;

  //////////////////////////////////////////
  //////// End of Configuration part ///////
  //////////////////////////////////////////


  TCanvas *c0 = new TCanvas("c0","canvas",W,H);
  //Fix the size of the pad: xmin,ymin,xmax,ymax
  TPad *pad = new TPad("up","up",0,0.27,1,1);
  pad->SetLeftMargin( L/W );
  pad->SetRightMargin( R/W );
  pad->SetTopMargin( T/H );
  pad->SetBottomMargin(0);
  pad->SetGridx();
  pad->SetGridy();
  pad->SetTicky();
  pad->Draw();

  c0->cd();
  TPad *pad2 =new TPad("down","down",0,0,1,0.27);
  //pad2->SetNumber(2);
  pad2->SetLeftMargin( L/W );
  pad2->SetRightMargin( R/W );
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin( B/H );
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->Draw();

  c0->cd();
  pad->cd();

  CMS_lumi(c0, 16, 0, true);
  TLegend *t = new TLegend(0.38,0.70,0.95,0.90);
  t->SetLineColor(0);
  t->SetBorderSize(1);
  t->SetNColumns(3);


  //TH1F* LO_MChistos[99]; //Only allow 99 MC processes
  std::vector<TH1F*> LO_MChistos(30); //Only allow 99 MC processes
  int iteHisto=0;
  bool firstPass = true;
  TH1F *stack_LO_MCsamples = new TH1F("hitst", "hist", 100, 0, 100); //= new TH1F("stack_LO_MCsamples",theHistoName);
  for (MCentry theMCentry: LO_MCsamples){
    if(VERBOSE) cout << "doing " << theMCentry.nameSample << endl;
    LO_MChistos[iteHisto] = (TH1F*) (theMCentry.sampleFile)->Get(theHistoName);
    if (LO_MChistos[iteHisto] == 0) continue;
    if(VERBOSE) cout << "found" << endl;
    TH1F *totEventInBaobab = (TH1F*) (theMCentry.sampleFile)->Get("totEventInBaobab_tot");
    float norm = instLumi*theMCentry.crossSection/totEventInBaobab->Integral();
    if(VERBOSE) cout << "scale is " << norm << endl;
    LO_MChistos[iteHisto]->Scale(norm);
    if(firstPass){
      stack_LO_MCsamples = LO_MChistos[iteHisto];
      stack_LO_MCsamples->SetLineColor(theMCentry.color);
      stack_LO_MCsamples->SetFillColor(theMCentry.color);
      t->AddEntry(stack_LO_MCsamples, theMCentry.legendEntry, "F");
      firstPass = false;
    }
    else{
      stack_LO_MCsamples->Add(LO_MChistos[iteHisto]);
    }
    delete totEventInBaobab;
    iteHisto++;
  }

  //TH1F* NLO_MChistos[99]; //Only allow 99 MC processes
  std::vector<TH1F*> NLO_MChistos(30);
  iteHisto=0;
  firstPass = true;
  TH1F *stack_NLO_MCsamples = new TH1F("hitst2", "hist2", 100, 0, 100); //= new TH1F("stack_LO_MCsamples",theHistoName);
  for (MCentry theMCentryB: NLO_MCsamples){
    if(VERBOSE) cout << "doing " << theMCentryB.nameSample << endl;
    NLO_MChistos[iteHisto] = (TH1F*) (theMCentryB.sampleFile)->Get(theHistoName);
    if (NLO_MChistos[iteHisto] == 0) continue;
    if(VERBOSE) cout << "found" << endl;
    TH1F *totEventInBaobab = (TH1F*) (theMCentryB.sampleFile)->Get("totEventInBaobab_tot");
    float norm = instLumi*theMCentryB.crossSection/totEventInBaobab->Integral();
    if(VERBOSE) cout << "scale is " << norm << endl;
    NLO_MChistos[iteHisto]->Scale(norm);
    if(firstPass){
      stack_NLO_MCsamples = NLO_MChistos[iteHisto];
      stack_NLO_MCsamples->SetLineColor(theMCentryB.color);
      stack_NLO_MCsamples->SetFillColor(theMCentryB.color);
      t->AddEntry(stack_NLO_MCsamples, theMCentryB.legendEntry, "F");
      firstPass = false;
    }
    else{
      stack_NLO_MCsamples->Add(NLO_MChistos[iteHisto]);
    }
    delete totEventInBaobab;
    iteHisto++;
  }

  c0->cd();
  pad->cd();

  stack_LO_MCsamples->Draw("");
  stack_NLO_MCsamples->Draw("same");

  t->Draw();
  TAxis* xaxis = new TAxis();
  TAxis* yaxis = new TAxis();
  stack_LO_MCsamples->SetMaximum(stack_LO_MCsamples->GetMaximum()*1.5);
  xaxis = stack_LO_MCsamples->GetXaxis();
  yaxis = stack_LO_MCsamples->GetYaxis();
  make_axis(xaxis, yaxis, fontType, pixelFontSize);
  gPad->RedrawAxis();

  TString theLeptonCategoryText, theJetCategoryText;
  theLeptonCategoryText = "#gamma";

  if(theHistoName.Contains("eq0jet")) theJetCategoryText = "0 jet";
  if(theHistoName.Contains("geq1jet")) theJetCategoryText = "#geq1 jet";
  if(theHistoName.Contains("vbf")) theJetCategoryText = "VBF";
  //if(theHistoName.Contains("vbf")) theJetCategoryText = "VBF-tagged";

  TLatex latex0;
  latex0.SetNDC();
  latex0.SetTextAngle(0);
  latex0.SetTextColor(kBlack);

  latex0.SetTextFont(42);
  latex0.SetTextSize(0.05*6/5.);
  latex0.DrawLatex(0.155,0.86,"2l2#nu");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);

  latex.SetTextFont(42);
  latex.SetTextSize(0.05*6/5.);
  latex.DrawLatex(0.155,0.805,"#it{"+theLeptonCategoryText+"}");

  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextAngle(0);
  latex2.SetTextColor(kBlack);

  latex2.SetTextSize(0.05);
  latex2.DrawLatex(0.215,0.804,theJetCategoryText);

  c0->cd();
  pad2->cd();
  TH1F *ratio = (TH1F*) stack_NLO_MCsamples->Clone("ratio");
  ratio->Sumw2();
  ratio->Divide(stack_LO_MCsamples);
  ratio->SetMaximum(2.0);
  ratio->SetMinimum(0.5);
  ratio->SetTitle("");
  ratio->GetYaxis()->SetTitle("#frac{NLO}{LO}");
  ratio->GetYaxis()->SetLabelFont(fontType);
  ratio->GetYaxis()->SetLabelSize(pixelFontSize);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetTitleFont(fontType);
  ratio->GetYaxis()->SetTitleSize(pixelFontSize);
  ratio->GetYaxis()->SetTitleOffset(2.1);
  ratio->GetXaxis()->SetLabelFont(fontType);
  ratio->GetXaxis()->SetLabelSize(pixelFontSize);
  ratio->GetXaxis()->SetTitleFont(fontType);
  ratio->GetXaxis()->SetTitleSize(pixelFontSize);
  ratio->GetXaxis()->SetTitleOffset(3.6); //4

  ratio->Draw("E1");
  TLine *line = new TLine(0,1.716910,1500,-0.11459);
  line->SetLineColor(kRed);
  line->Draw();

  TString outputDir = base_path + "WeightsAndDatadriven/";
  c0->Print(outputDir+theHistoName+".png");
  c0->Print(outputDir+theHistoName+".root");
  pad->cd();
  pad->cd()->SetLogy();
  stack_LO_MCsamples->SetMaximum(stack_LO_MCsamples->GetMaximum()*100.);
  stack_LO_MCsamples->SetMinimum(0.05);
  stack_LO_MCsamples->Draw("E1:same");
  c0->Print(outputDir+theHistoName+"_log.png");


	ofstream outputFile;
  outputFile.open(outputDir + "GJetsNLOWeights.C");
  TString bin_X_name = "boson.Pt()"; //The name you want for your X axis

	int nBins = ratio->GetNbinsX();
  double binLow_X = 0;
  double binUp_X = 0;
  for(int bin_X = 1; bin_X <= nBins; bin_X++){
    binLow_X = ratio->GetXaxis()->GetBinLowEdge(bin_X);
    binUp_X = binLow_X + ratio->GetXaxis()->GetBinWidth(bin_X);
    writeCondition(outputFile, bin_X_name, binLow_X, binUp_X, bin_X);
 		writeEfficiency(outputFile, ratio->GetBinContent(bin_X));
  }
  outputFile.close();




	delete c0;



}


void computeGJetsNLOWeight(TString suffix){
  TString currentDirectory= base_path + "OUTPUTS/"+suffix+"/MERGED";
  gROOT->SetBatch();
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(kTRUE); //To ensure that all histograms are created with the sum of weights
  TH2::SetDefaultSumw2(kTRUE); //To ensure that all histograms are created with the sum of weights
  gErrorIgnoreLevel = kWarning;//kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal;
  if(VERBOSE) gErrorIgnoreLevel = kPrint;
  //gStyle->SetOptTitle(0);

  //This will be the order to draw, so smallest XS should be first line
  std::vector<MCentry> LO_MCsamples;
  std::vector<MCentry> NLO_MCsamples;

  takeHisto_LO(LO_MCsamples);
  takeHisto_NLO(NLO_MCsamples);
  instLumi= 35866.932; //FIXME This value comes from the old code, from the baobabs it is: 16916.39+20092.7;//Lumi of runs GH + BtoF for singlePhoton samples;
  outputPrefixName = "outputInstrMET_";

  for (MCentry &theEntry: LO_MCsamples) theEntry.sampleFile = new TFile(currentDirectory+"/"+outputPrefixName+theEntry.fileSuffix+".root");
  for (MCentry &theEntry: NLO_MCsamples) theEntry.sampleFile = new TFile(currentDirectory+"/"+outputPrefixName+theEntry.fileSuffix+".root");

  drawTheHisto(LO_MCsamples, NLO_MCsamples, "pT_Z_SpikesRemovedAndSigmaCut_passQCDveto_ZMetPhi", suffix);

}
