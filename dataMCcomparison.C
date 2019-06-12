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
#include "Tools/CMS_lumi.C"
#include "samples.h"

#define VERBOSE true


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

void updateListOfPlots(std::map<TString, TString> & listOfHisto, TFile * file){
  TIter listPlots(file->GetListOfKeys());
  TKey *keyPlot;
  while ((keyPlot = (TKey*)listPlots())) {
    if( ((TString) keyPlot->GetTitle()).EndsWith("_up") || ((TString) keyPlot->GetTitle()).EndsWith("_down") ) continue;
    listOfHisto.insert(std::pair<TString, TString>(keyPlot->GetTitle(), keyPlot->GetClassName())); //Title and Type of the plot. This is done only for one occurence of the Title
  }
  delete keyPlot;
}

void progressbar(float progress){
  int barWidth = 70;

  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
}

void make_axis(TAxis* & xaxis, TAxis* & yaxis, int fontType, int pixelFontSize){
  yaxis->SetLabelFont(fontType);
  yaxis->SetLabelSize(pixelFontSize);
  //yaxis->SetNdivisions(505);
  yaxis->SetTitleOffset(2.1); //1.4
  yaxis->SetTitleFont(fontType);
  yaxis->SetTitleSize(pixelFontSize);

}

void drawTheHisto(TFile *dataFile, std::vector<MCentry> allSignals, std::vector<MCentry> allMCsamples, TString theHistoName, TString suffix, TString typeObject, TString analysisType){
  gROOT->SetBatch();
  if(typeObject.Contains("TH1")) typeObject = "TH1";
  else if(typeObject.Contains("TH2")) typeObject = "TH2";
  else{
    cout << "/!\\ WARNING /!\\ You have an histogram that is not a TH1 or a TH2 and I don't know how to draw it... so I don't." << endl;
    return;
  }
  if(VERBOSE) cout<< "In draw the histo for "<<theHistoName<<endl;
  bool dataExist = dataFile->GetListOfKeys()->Contains(theHistoName);
  TH1F *MZ_data = (TH1F*) dataFile->Get(theHistoName);
  TH1F *totEventInBaobab_tot_data = (TH1F*) dataFile->Get("totEventInBaobab_tot");
  if(VERBOSE) cout << "the tot events data =" << totEventInBaobab_tot_data->Integral() << endl;


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
  pad2->SetNumber(2);
  pad2->SetLeftMargin( L/W );
  pad2->SetRightMargin( R/W );
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin( B/H );
  pad2->SetGridx();
  pad2->SetGridy();
  if(dataExist) pad2->Draw();

  c0->cd();
  pad->cd();

  CMS_lumi(c0, 16, 0, true);
  //TLegend *t = new TLegend(0.69,0.66,0.96,0.91);
  TLegend *t = new TLegend(0.38,0.70,0.95,0.90);
  t->SetLineColor(0);
  t->SetBorderSize(1);
  t->SetNColumns(3);
  if(dataExist) t->AddEntry(MZ_data, "Data", "l");

  if(typeObject== "TH1"){
    if(dataExist) MZ_data->SetMarkerColor(kBlack);
    if(dataExist) MZ_data->SetLineColor(kBlack);
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
    if(theMCentry.crossSection != 0) MChistos[iteHisto]->Scale(norm);
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
  // Part for plotting the signal at 800 GeV
  TH1F *signalHisto[allSignals.size()];
  int iteSignal=0;
  for (MCentry signalEntry: allSignals){
    if(VERBOSE) cout << "doing signal" << signalEntry.nameSample << endl;
    signalHisto[iteSignal] = (TH1F*) (signalEntry.sampleFile)->Get(theHistoName);
    if (signalHisto[iteSignal] == 0) continue;
    if(VERBOSE) cout << "found" << endl;
    TH1F *totEventInBaobab = (TH1F*) (signalEntry.sampleFile)->Get("totEventInBaobab_tot");
    float norm = instLumi*signalEntry.crossSection/totEventInBaobab->Integral();
    if(VERBOSE) cout << "scale is " << norm << endl;
    if(VERBOSE) cout << "normalization before is " << signalHisto[iteSignal]->Integral() << endl;
    if(signalEntry.crossSection != 0) signalHisto[iteSignal]->Scale(norm);
    if(VERBOSE) cout << "normalization after is " << signalHisto[iteSignal]->Integral() << endl;
    if(typeObject== "TH1"){
      signalHisto[iteSignal]->SetLineColor(signalEntry.color);
      signalHisto[iteSignal]->SetLineWidth(2);
      signalHisto[iteSignal]->SetLineStyle(1);
    }
    else if(typeObject== "TH2") signalHisto[iteSignal]->SetLineColor(kBlack);
    t->AddEntry(signalHisto[iteSignal], signalEntry.legendEntry, "L");
    delete totEventInBaobab;
    iteSignal++;
  }

  if(iteHisto==0){
    std::cout << "No MC for this plot, not drawing it : " << theHistoName << std::endl;
    return;
  }
  if(theHistoName == "metFilters_tot") doMetFilterEfficiencyPlots(MZ_data, stackMCsamples);

  c0->cd();
  pad->cd();

  if(typeObject== "TH1"){
    if(dataExist){
      MZ_data->SetTitle("");
      MZ_data->Draw("E1:same");
      stackMCsamples->Draw("HIST:same");
      for (int iteSignal = 0; iteSignal < allSignals.size(); iteSignal++) {
        signalHisto[iteSignal]->Draw("HIST:same");
      }
      MZ_data->Draw("E1:same");
    }
    else{
      stackMCsamples->SetTitle("");
      stackMCsamples->Draw("HIST");
    }
  }
  else if(typeObject== "TH2"){
    stackMCsamples->SetTitle("");
    stackMCsamples->Draw("");
    if(dataExist) MZ_data->Draw("LEGO:same");
  }
  if(VERBOSE) cout << "Plot drawn. " << endl;
  t->Draw();
  if(VERBOSE) cout << "Legend drawn. " << endl;
  TAxis* xaxis = new TAxis();
  TAxis* yaxis = new TAxis();
  if(dataExist){
    MZ_data->SetMaximum(MZ_data->GetMaximum()*1.5);
    xaxis = MZ_data->GetXaxis();
    yaxis = MZ_data->GetYaxis();
  }
  else{
    stackMCsamples->SetMaximum(stackMCsamples->GetMaximum()*1.5);
    xaxis = stackMCsamples->GetXaxis();
    yaxis = stackMCsamples->GetYaxis();
  }
  make_axis(xaxis, yaxis, fontType, pixelFontSize);
  gPad->RedrawAxis();

  TString theLeptonCategoryText = "", theJetCategoryText = "";
  if(theHistoName.Contains("ee")) theLeptonCategoryText = "ee";
  else if(theHistoName.Contains("mumu")) theLeptonCategoryText = "#mu#mu";
  else if(theHistoName.Contains("ll")) theLeptonCategoryText = "ll";
  else if(theHistoName.Contains("gamma")) theLeptonCategoryText = "#gamma";

  if(analysisType == "InstrMET") theLeptonCategoryText = "#gamma";

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
  //if(theJetCategoryText.Contains("VBF")) latex2.DrawLatex(0.155,0.68,theJetCategoryText);
  //else latex2.DrawLatex(0.215,0.754,theJetCategoryText);
  latex2.DrawLatex(0.23,0.804,theJetCategoryText);

  pad2->cd();
  if(dataExist){
    TH1F *ratio = (TH1F*) MZ_data->Clone("ratio");
    ratio->Sumw2();
    ratio->Divide(MZ_data, ((TH1F*)stackMCsamples->GetStack()->Last()), 1,1);
    if(typeObject== "TH1"){
      //ratio->SetMaximum(1.3);
      //ratio->SetMinimum(0.7);
      ratio->SetMaximum(3.0);
      ratio->SetMinimum(0.);
    }
    else if(typeObject== "TH2"){
      ratio->SetMaximum(1.5);
      ratio->SetMinimum(0.5);

    }
    ratio->SetTitle("");
    ratio->GetYaxis()->SetTitle("#frac{Data}{#Sigma Bkg.}");
    //ratio->GetYaxis()->SetLabelSize(0.08);
    //ratio->GetXaxis()->SetTitleSize(0.12);
    //ratio->GetXaxis()->SetLabelSize(0.1);
    //ratio->GetXaxis()->SetLabelOffset(0.02);
    ratio->GetYaxis()->SetLabelFont(fontType);
    ratio->GetYaxis()->SetLabelSize(pixelFontSize);
    //ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->CenterTitle();
    ratio->GetYaxis()->SetTitleFont(fontType);
    ratio->GetYaxis()->SetTitleSize(pixelFontSize);
    ratio->GetYaxis()->SetTitleOffset(2.1);
    ratio->GetXaxis()->SetLabelFont(fontType);
    ratio->GetXaxis()->SetLabelSize(pixelFontSize);
    ratio->GetXaxis()->SetTitleFont(fontType);
    ratio->GetXaxis()->SetTitleSize(pixelFontSize);
    ratio->GetXaxis()->SetTitleOffset(3.6); //4


    if(typeObject== "TH1") ratio->Draw("E1");
    else if(typeObject== "TH2") ratio->Draw("LEGO");
  }
  TString outputDir = "OUTPUTS/"+suffix+"/PLOTS/";
  c0->Print(outputDir+theHistoName+".png");
  c0->Print(outputDir+theHistoName+".root");
  pad->cd();
  if(typeObject== "TH1") pad->cd()->SetLogy();
  else if(typeObject== "TH2") pad->cd()->SetLogz();
  if(dataExist){
    MZ_data->SetMaximum(MZ_data->GetMaximum()*100000.);
    MZ_data->SetMinimum(0.0002);
    MZ_data->Draw("E1:same");
  }
  c0->Print(outputDir+theHistoName+"_log.png");
  delete c0;
  //  for (int i=0 ; i < iteHisto ; i++){  delete MChistos[i];}
  //  delete stackMCsamples;
}


void dataMCcomparison(TString analysisType, TString suffix, TString doMELA){
  TString currentDirectory="OUTPUTS/"+suffix+"/MERGED";
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(kTRUE); //To ensure that all histograms are created with the sum of weights
  TH2::SetDefaultSumw2(kTRUE); //To ensure that all histograms are created with the sum of weights
  gErrorIgnoreLevel = kWarning;//kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal;
  if(VERBOSE) gErrorIgnoreLevel = kPrint;
  //gStyle->SetOptTitle(0);

  systSuffixName = ""; //we don't look at syst in the dataMCcomparison script.
  std::vector<MCentry> allMCsamples;
  TFile* dataFile = new TFile();
  std::vector<MCentry> allSignals;
  bool MELA = false;
  if (doMELA == "--doMelaReweight") MELA = true;
  if(analysisType == "HZZanalysis"){
    outputPrefixName = "outputHZZ_";
    takeHisto_HZZanalysis(allMCsamples, &dataFile, allSignals, currentDirectory, MELA);
  }  
  else if(analysisType == "InstrMET"){
    outputPrefixName = "outputInstrMET_";
    takeHisto_InstrMET(allMCsamples, &dataFile, currentDirectory);
  }
  else if(analysisType == "HZZdatadriven"){
    bool isDatadriven = true;
    outputPrefixName = "outputHZZ_";
    
    takeHisto_HZZanalysis(allMCsamples, &dataFile, allSignals, currentDirectory,  MELA, isDatadriven);
  }
  else if(analysisType == "NRB"){
    bool isDatadriven = true;
    outputPrefixName = "outputNRB_";
    takeHisto_NRB(allMCsamples, &dataFile, allSignals, currentDirectory, MELA, isDatadriven);
  }

  for (MCentry &theEntry: allMCsamples){
    theEntry.sampleFile = new TFile(currentDirectory+"/"+outputPrefixName+theEntry.fileSuffix+".root");
  }
  for (MCentry &signalEntry: allSignals){
    signalEntry.sampleFile = new TFile(currentDirectory+"/"+outputPrefixName+signalEntry.fileSuffix+".root");
  }
  //make list of histo from data and MC
  std::map<TString, TString> listOfHisto; //A map containing the name of the histo. First element is the name of the histo and the second is its type
  updateListOfPlots(listOfHisto, dataFile);
  for (MCentry theMCentry: allMCsamples) updateListOfPlots(listOfHisto, theMCentry.sampleFile);

  //Remove unwanted histograms here:
  listOfHisto.erase("totEventInBaobab");
  listOfHisto.erase("totEventInBaobab_tot");


  float index = 0;
  for (std::pair<TString, TString> element : listOfHisto) {
    TString histoName = element.first;
    TString typeName = element.second;
    if(VERBOSE) cout << "Type:" << typeName << " and title:" << histoName << endl;
    drawTheHisto(dataFile, allSignals, allMCsamples, histoName, suffix, typeName, analysisType);
    index++;
    if(!VERBOSE) progressbar( index/(1.*listOfHisto.size()));
  }
  std::cout << "\nDone." << std::endl;

  //drawTheHisto(dataFile, allMCsamples, "eventflow_tot", suffix, "TH1", analysisType);




}
