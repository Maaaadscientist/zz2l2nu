#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TMath.h>
#include <math.h>
#include <vector>
#include <TFile.h>
#include <iostream>
#include "TROOT.h"
#include <TGaxis.h>
#include <TLine.h>
#include <algorithm> 
#include <map>
#include <TStyle.h>
#include "../../Tools/CMS_lumi.C"
#include <type_traits>

//General options
#define VERBOSE true
#define HIDE_WARNING true
#define DEBUG_HISTOS true

int addQCD = 0;

struct MCentry{
  TString nameSample;
  TString legendEntry;
  TString fileSuffix;
  int color;
  TFile *sampleFile;
  MCentry(TString theNameSample, TString theLegendEntry, TString theFileSuffix, int theColor)
    : nameSample(theNameSample), legendEntry(theLegendEntry), fileSuffix(theFileSuffix), color(theColor)
  {}
};

//FIXME This has to be updated with the latest value in dataMCComparison.C script
void takeHisto_GJets(std::vector<MCentry> & allMCsamples){
  //GJets_HT
  allMCsamples.push_back(MCentry("GJets_HT-100To200", "#gamma+jets LO", "GJets_HT-100To200", 93));
  allMCsamples.push_back(MCentry("GJets_HT-200To400", "#gamma+jets LO", "GJets_HT-200To400", 93));
  allMCsamples.push_back(MCentry("GJets_HT-400To600", "#gamma+jets LO", "GJets_HT-400To600", 93));
  allMCsamples.push_back(MCentry("GJets_HT-40To100", "#gamma+jets LO", "GJets_HT-40To100", 93));
  allMCsamples.push_back(MCentry("GJets_HT-600ToInf", "#gamma+jets LO", "GJets_HT-600ToInf", 93));
  //QCD_HT
  if(addQCD ==1){
    allMCsamples.push_back(MCentry("QCD_HT1000to1500", "QCD, HT>100", "QCD_HT1000to1500", 21));
    allMCsamples.push_back(MCentry("QCD_HT100to200", "QCD, HT>100", "QCD_HT100to200", 21)); //Sample with low stats 
    allMCsamples.push_back(MCentry("QCD_HT1500to2000", "QCD, HT>100", "QCD_HT1500to2000", 21));
    allMCsamples.push_back(MCentry("QCD_HT2000toInf", "QCD, HT>100", "QCD_HT2000toInf", 21));
    allMCsamples.push_back(MCentry("QCD_HT200to300", "QCD, HT>100", "QCD_HT200to300", 21));
    allMCsamples.push_back(MCentry("QCD_HT300to500", "QCD, HT>100", "QCD_HT300to500", 21));
    allMCsamples.push_back(MCentry("QCD_HT500to700", "QCD, HT>100", "QCD_HT500to700", 21));
    allMCsamples.push_back(MCentry("QCD_HT700to1000", "QCD, HT>100", "QCD_HT700to1000", 21));
  }
}

void takeHisto_DY(std::vector<MCentry> & allMCsamples){
  //DY
  allMCsamples.push_back(MCentry("DY",        "DY",    "DYJetsToLL_M-50",   833));
}

template<class T> T* newHisto();

  template<>
TH1F* newHisto<TH1F>()
{
  return new TH1F("hist", "hist", 100, 0, 100);
}
  template<>
TH2F* newHisto<TH2F>()
{
  return new TH2F("hist", "hist", 100, 0, 100, 100, 0, 100);
}
  template<>
TProfile* newHisto<TProfile>()
{
  return new TProfile("hist", "hist", 100, 0, 100, 0, 100);
}

template<class T>
T* getHistoFromMCentries(std::vector<MCentry> MCentries, TString theHistoName){
  std::vector<T*> MChistos(99); //Only allow 99 MC processes
  int iteHisto=0;
  bool firstPass = true;
  T *stack_MCsamples = newHisto<T>(); 
  for (MCentry theMCentry: MCentries){
    MChistos[iteHisto] = (T*) (theMCentry.sampleFile)->Get(theHistoName);
    if (MChistos[iteHisto] == 0) continue;
    if(firstPass){
      stack_MCsamples = MChistos[iteHisto];
      stack_MCsamples->SetLineColor(theMCentry.color);
      stack_MCsamples->SetFillColor(theMCentry.color);
      firstPass = false;
    }
    else{
      stack_MCsamples->Add(MChistos[iteHisto]);
    }
    iteHisto++;
  }

  return stack_MCsamples;
}

void make_axis(TAxis* & xaxis, TAxis* & yaxis, int fontType, int pixelFontSize){
  yaxis->SetLabelFont(fontType);
  yaxis->SetLabelSize(pixelFontSize);
  //yaxis->SetNdivisions(505);
  yaxis->SetTitleOffset(2.1); //1.4
  yaxis->SetTitleFont(fontType);
  yaxis->SetTitleSize(pixelFontSize);

}
template<class Type>
void drawHistoClosure(Type* HZZ, Type* InstrMET, TString base_path, TFile *f_output, TString theHistoName, TString jetCat, TString lepCat, TString suffix  = ""){

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

  TCanvas *c0 = new TCanvas(theHistoName+suffix,theHistoName+suffix,W,H);
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
  TLegend *t = new TLegend(0.7,0.70,0.95,0.90);
  t->SetLineColor(0);
  t->SetBorderSize(1);
  //t->SetNColumns(3);
  t->AddEntry(HZZ, "Z+jets", "F");
  if(addQCD ==1) t->AddEntry(InstrMET, "#gamma/QCD reweighted", "F");
  else t->AddEntry(InstrMET, "#gamma reweighted", "F");

  c0->cd();
  pad->cd();

  InstrMET->SetStats(0);
  HZZ->SetStats(0);

  InstrMET->Draw("E1");
  HZZ->Draw("E1:same");
  t->Draw("same");

  TAxis* xaxis = new TAxis();
  TAxis* yaxis = new TAxis();
  InstrMET->SetMaximum(InstrMET->GetMaximum()*1.5);
  xaxis = InstrMET->GetXaxis();
  yaxis = InstrMET->GetYaxis();
  make_axis(xaxis, yaxis, fontType, pixelFontSize);
  gPad->RedrawAxis();

  TString theLeptonCategoryText, theJetCategoryText;

  if(lepCat.Contains("ee")) theLeptonCategoryText = "ee";
  else if(lepCat.Contains("mumu")) theLeptonCategoryText = "#mu#mu";
  else if(lepCat.Contains("ll")) theLeptonCategoryText = "ll";

  if(jetCat.Contains("eq0jet")) theJetCategoryText = "0 jet";
  else if(jetCat.Contains("geq1jet")) theJetCategoryText = "#geq1 jet";
  else if(jetCat.Contains("vbf")) theJetCategoryText = "VBF";

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
  latex2.DrawLatex(0.245,0.804,theJetCategoryText);


  c0->cd();
  pad2->cd();
  TH1F *ratio = (TH1F*) HZZ->Clone("ratio");
  ratio->Sumw2();
  ratio->Divide(InstrMET);
  ratio->SetMaximum(2.0);
  ratio->SetMinimum(0.5);
  ratio->SetTitle("");
  ratio->GetYaxis()->SetTitle("#frac{Z+jets}{#gamma+jets}");
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
  TString outputDir = base_path + "WeightsAndDatadriven/InstrMET/";
  c0->Print(outputDir+theHistoName+suffix+jetCat+lepCat+".png");

  c0->cd();
  pad->cd();
  pad->cd()->SetLogy();
  InstrMET->SetMaximum(InstrMET->GetMaximum()*100.);
  InstrMET->SetMinimum(0.05);
  HZZ->Draw("E1:same");
  c0->Print(outputDir+theHistoName+suffix+jetCat+lepCat+"_log.png");

  f_output->cd();
  c0->Write();

}

void step1_weight_NVtx(TString base_path, std::vector<MCentry> DY_MCsamples, std::vector<MCentry> GJets_MCsamples){
  if(VERBOSE) std::cout<< "Launching step1 of the computation of the weights for Closure Test: NVtx" << std::endl; 
  // Initialization: opening and creating ROOT files
  TFile *f_output = new TFile(base_path+"WeightsAndDatadriven/InstrMET/closureTest_weight_NVtx.root","RECREATE");

  std::vector<TString> jetCat = {"_eq0jets","_geq1jets","_vbf", ""};
  std::vector<TString> lepCat = {"_ee","_mumu", "_ll"};
  std::vector<std::vector<TH1F*> > HZZ_nvtx(jetCat.size());
  std::vector<TH1F*> InstrMET_nvtx;
  std::vector<TH1F*> results_1D;

  for(unsigned int i = 0; i < jetCat.size(); i++){
    InstrMET_nvtx.push_back( getHistoFromMCentries<TH1F>(GJets_MCsamples, "reco-vtx_MET125_InstrMET_reweighting"+jetCat[i]+"_gamma"));
    //InstrMET_nvtx[i]->Scale(1./InstrMET_nvtx[i]->Integral());
    if(DEBUG_HISTOS) InstrMET_nvtx[i]->Write();
    for(unsigned int j = 0; j < lepCat.size(); j++){
      HZZ_nvtx[i].push_back( (TH1F*) getHistoFromMCentries<TH1F>(DY_MCsamples, "reco-vtx_MET125_InstrMET_reweighting"+jetCat[i]+lepCat[j]));
      //HZZ_nvtx[i][j]->Scale(1./HZZ_nvtx[i][j]->Integral());
      results_1D.push_back( (TH1F*) HZZ_nvtx[i][j]->Clone("WeightHisto"+jetCat[i]+lepCat[j]+"_AllBins"));
      if(results_1D.back()->GetEntries() != 0 ) results_1D.back()->Divide(InstrMET_nvtx[i]);
      results_1D.back()->Write();
      if(DEBUG_HISTOS) HZZ_nvtx[i][j]->Write(); //Normalized and used to compute weights
    }
  }
  if(VERBOSE) std::cout<< "Step1 of reweighting done. Weights are available under: closureTest_weight_NVtx.root" << std::endl; 
}

void step2_weight_pt(TString base_path,  std::vector<MCentry> DY_MCsamples, std::vector<MCentry> GJets_MCsamples) {
  if(VERBOSE) std::cout<< "Launching step2 of the computation of the weights for Closure Test: Pt" << std::endl; 

  // Initialization: opening and creating ROOT files
  TFile *f_output = new TFile(base_path+"WeightsAndDatadriven/InstrMET/closureTest_weight_pt.root","RECREATE");

  std::vector<TString> jetCat = {"_eq0jets","_geq1jets","_vbf", ""};
  std::vector<TString> lepCat = {"_ee","_mumu", "_ll"};

  std::vector<std::vector<TH1F*> > HZZ_pt(jetCat.size());
  std::vector<std::vector<TH1F*> > InstrMET_pt(jetCat.size());
  std::vector<TH1F*> results_1D;
    for(unsigned int i = 0; i < jetCat.size(); i++){
      for(unsigned int j = 0; j < lepCat.size(); j++){
        InstrMET_pt[i].push_back( getHistoFromMCentries<TH1F>(GJets_MCsamples, "pT_Boson_MET125_InstrMET_reweightingAfter"+lepCat[j]+jetCat[i]));
        //InstrMET_pt[i][j]->Scale(1./InstrMET_pt[i][j]->Integral()); //In order to have a good renormallization between gamma and dilepton data, we don't scale their weights
        if(DEBUG_HISTOS) InstrMET_pt[i][j]->Write();
        HZZ_pt[i].push_back( getHistoFromMCentries<TH1F>(DY_MCsamples, "pT_Boson_MET125_InstrMET_reweighting"+jetCat[i]+lepCat[j]));
        //HZZ_pt[i][j]->Scale(1./HZZ_pt[i][j]->Integral()); //In order to have a good renormallization between gamma and dilepton data, we don't scale their weights
        results_1D.push_back( (TH1F*) HZZ_pt[i][j]->Clone("WeightHisto"+jetCat[i]+lepCat[j]+"_AllBins"));
        if(results_1D.back()->GetEntries() != 0 ) results_1D.back()->Divide(InstrMET_pt[i][j]);
        results_1D.back()->Write();
        if(DEBUG_HISTOS) HZZ_pt[i][j]->Write(); //Normalized and used to compute weights
      }
    }
  if(VERBOSE) std::cout<< "Step2 of reweighting done. Weights are available under: closureTest_weight_pt.root" << std::endl;

}

void step3_weight_lineshape(TString base_path, std::vector<MCentry> DY_MCsamples) {
  if(VERBOSE) std::cout<< "Launching step3 of the computation of the weights for Closure Test: mass lineshape" << std::endl;

  TFile *f_output = new TFile(base_path+"WeightsAndDatadriven/InstrMET/closureTest_lineshape_mass.root","RECREATE");

  std::vector<TH1F*> results_1D;

  std::vector<TString> jetCat = {"_eq0jets","_geq1jets","_vbf", ""};
  std::vector<TString> lepCat = {"_ee","_mumu", "_ll"};

  for(unsigned int i = 0; i < jetCat.size(); i++){
    for(unsigned int j = 0; j < lepCat.size(); j++){
      results_1D.push_back( (TH1F*) ((TH1F*) getHistoFromMCentries<TH1F>(DY_MCsamples, "M_Boson_MET125_InstrMET_reweighting"+jetCat[i]+lepCat[j]))->Clone("WeightHisto"+jetCat[i]+lepCat[j]+"_AllBins") );
      results_1D.back()->Write();
    }
  }
  if(VERBOSE) std::cout<< "Step3 of reweighting done. Weights are available under: closureTest_lineshape_mass.root" << std::endl;

}

void step4_closure(TString base_path,  std::vector<MCentry> DY_MCsamples, std::vector<MCentry> GJets_MCsamples) {
  if(VERBOSE) std::cout<< "Launching step4 of the computation of the weights: Closure Test" << std::endl;

  TFile *f_output = new TFile(base_path+"WeightsAndDatadriven/InstrMET/closureTestResults.root","RECREATE");

  std::vector<TString> plotsForClosure_TH1 = {"reco-vtx", "pT_Boson", "M_Boson", "MET", "mT", "mT_closure", "DeltaPhi_MET_Boson", "DeltaPhi_MET_Jet", "METoverPt_zoom", "eta_Boson", "pT_jet0", "nJets", "selJetsHT", "MET_HT300", "MET_Pt0-300", "MET_Pt300-400", "MET_Pt400-600", "MET_Pt600-Inf", "rho", "MET_HT0-100", "MET_HT100-200", "MET_HT200-300", "MET_HT300-500", "MET_HT500-1000", "MET_HT1000-Inf", "MET_dPhiMetBosons0p5-1", "MET_dPhiMetBosons1-2", "MET_dPhiMetBosons2-Inf"};
  std::vector<TString> plotsForClosure_TProfile = {"METvsBosonPt", "METvsMT", "METvsDPhiMETBos", "METvsDPhiMETJet", "METvsJetPt", "METvsNJets", "METvsBosonEta", "METvsHT", "HTvsBosonEta"};
  std::vector<TString> jetCat = {"_eq0jets","_geq1jets","_vbf", ""};
  std::vector<TString> lepCat = {"_ee","_mumu", "_ll"};

  std::vector<std::vector<std::vector<TH1F*> > > HZZ(plotsForClosure_TH1.size(), std::vector<std::vector<TH1F*> >(jetCat.size()/*, std::vector<TH1F*>*/)) ;
  std::vector<std::vector<std::vector<TH1F*> > > InstrMET(plotsForClosure_TH1.size(), std::vector<std::vector<TH1F*> >(jetCat.size()/*, std::vector<TH1F*>*/)) ;
  for(unsigned int h = 0; h < plotsForClosure_TH1.size(); h++){
    for(unsigned int i = 0; i < jetCat.size(); i++){
      for(unsigned int j = 0; j < lepCat.size(); j++){
        InstrMET[h][i].push_back( getHistoFromMCentries<TH1F>(GJets_MCsamples, plotsForClosure_TH1[h]+"_InstrMET_AllWeightsAndLineshapeApplied"+lepCat[j]+jetCat[i]+"_gamma"));
        HZZ[h][i].push_back( getHistoFromMCentries<TH1F>(DY_MCsamples, plotsForClosure_TH1[h]+"_InstrMET_reweighting"+jetCat[i]+lepCat[j]));
        InstrMET[h][i][j]->Scale(HZZ[h][i][j]->Integral()/InstrMET[h][i][j]->Integral());
        if(DEBUG_HISTOS) InstrMET[h][i][j]->Write();
        if(DEBUG_HISTOS) HZZ[h][i][j]->Write();
        drawHistoClosure<TH1F>(HZZ[h][i][j], InstrMET[h][i][j], base_path, f_output, plotsForClosure_TH1[h], jetCat[i], lepCat[j]);
      }
    }
  }

  std::vector<std::vector<std::vector<TProfile*> > > HZZ_TProfile(plotsForClosure_TProfile.size(), std::vector<std::vector<TProfile*> >(jetCat.size()/*, std::vector<TProfile*>*/)) ;
  std::vector<std::vector<std::vector<TProfile*> > > InstrMET_TProfile(plotsForClosure_TProfile.size(), std::vector<std::vector<TProfile*> >(jetCat.size()/*, std::vector<TProfile*>*/)) ;
  for(unsigned int h = 0; h < plotsForClosure_TProfile.size(); h++){
    for(unsigned int i = 0; i < jetCat.size(); i++){
      for(unsigned int j = 0; j < lepCat.size(); j++){
        InstrMET_TProfile[h][i].push_back( getHistoFromMCentries<TProfile>(GJets_MCsamples, plotsForClosure_TProfile[h]+"_InstrMET_AllWeightsAndLineshapeApplied"+lepCat[j]+jetCat[i]+"_gamma"));
        if(DEBUG_HISTOS) InstrMET_TProfile[h][i][j]->Write();
        HZZ_TProfile[h][i].push_back( getHistoFromMCentries<TProfile>(DY_MCsamples, plotsForClosure_TProfile[h]+"_InstrMET_reweighting"+jetCat[i]+lepCat[j]));
        if(DEBUG_HISTOS) HZZ_TProfile[h][i][j]->Write();
        drawHistoClosure<TProfile>(HZZ_TProfile[h][i][j], InstrMET_TProfile[h][i][j], base_path, f_output, plotsForClosure_TProfile[h], jetCat[i], lepCat[j]);
      }
    }
  }

  std::vector<TString> plots_final = {"M_Boson", "pT_Boson", "MET", "mT", "nJets", "eta_Boson", "DeltaPhi_MET_Jet", "DeltaPhi_MET_Boson", "reco-vtx", "METphi"};
  std::vector<std::vector<std::vector<TH1F*> > > HZZ_final(plots_final.size(), std::vector<std::vector<TH1F*> >(jetCat.size()/*, std::vector<TH1F*>*/)) ;
  std::vector<std::vector<std::vector<TH1F*> > > InstrMET_final(plots_final.size(), std::vector<std::vector<TH1F*> >(jetCat.size()/*, std::vector<TH1F*>*/)) ;
  for(unsigned int h = 0; h < plots_final.size(); h++){
    for(unsigned int i = 0; i < jetCat.size(); i++){
      for(unsigned int j = 0; j < lepCat.size(); j++){
        InstrMET_final[h][i].push_back( getHistoFromMCentries<TH1F>(GJets_MCsamples, plots_final[h]+"_final"+lepCat[j]+jetCat[i]+"_gamma"));
        HZZ_final[h][i].push_back( getHistoFromMCentries<TH1F>(DY_MCsamples, plots_final[h]+"_final"+jetCat[i]+lepCat[j]));
        InstrMET_final[h][i][j]->Scale(HZZ_final[h][i][j]->Integral()/InstrMET_final[h][i][j]->Integral());
        if(DEBUG_HISTOS) InstrMET_final[h][i][j]->Write();
        if(DEBUG_HISTOS) HZZ_final[h][i][j]->Write();
        drawHistoClosure<TH1F>(HZZ_final[h][i][j], InstrMET_final[h][i][j], base_path, f_output, plots_final[h], jetCat[i], lepCat[j], "_final");
      }
    }
  }


  if(VERBOSE) std::cout<< "Closure test done." << std::endl;
}

void macroToComputeClosureTestWeights(int reweightingStep, int addQCD){
  gROOT->SetBatch();
  TH1::SetDefaultSumw2(kTRUE); //To ensure that all histograms are created with the sum of weights
  if(HIDE_WARNING) gErrorIgnoreLevel=kError;

  TString GJets_suffix = "";
  if(reweightingStep == 1){
    GJets_suffix = "closureTest_PhotonMC_NoWeight";
  }
  else if(reweightingStep == 2){
    GJets_suffix = "closureTest_PhotonMC_NVtx_WeightApplied";
  }
  else if(reweightingStep == 3){
    GJets_suffix = "closureTest_PhotonMC_NVtx_WeightApplied"; //Not really needed actually since step 3 doesn't look at GJets... But like this it's implemented if needed.
  }
  else if(reweightingStep == 4){
    GJets_suffix = "closureTest_PhotonMC_AllWeightsAndLineshapeApplied";
  }

  //Define all paths
  TString base_path = std::string(getenv("HZZ2L2NU_BASE")) + "/";
  TString DY_path = base_path + "OUTPUTS/closureTest_DiLeptonMC/MERGED/"; // The path of the DiLepton MC used for reweighting
  TString GJets_path = base_path + "OUTPUTS/"+GJets_suffix+"/MERGED/"; // The path of the GJets MC used for reweighting

  //Open and create ROOT files
  std::vector<MCentry> DY_MCsamples;
  std::vector<MCentry> GJets_MCsamples;

  takeHisto_DY(DY_MCsamples);
  takeHisto_GJets(GJets_MCsamples);
  for (MCentry &theEntry: DY_MCsamples) theEntry.sampleFile = new TFile(DY_path+"outputHZZ_"+theEntry.fileSuffix+".root");
  for (MCentry &theEntry: GJets_MCsamples) theEntry.sampleFile = new TFile(GJets_path+"outputInstrMET_"+theEntry.fileSuffix+".root");

  if(reweightingStep == 1) step1_weight_NVtx(base_path, DY_MCsamples, GJets_MCsamples);
  else if (reweightingStep == 2) step2_weight_pt(base_path, DY_MCsamples, GJets_MCsamples);
  else if (reweightingStep == 3) step3_weight_lineshape(base_path, DY_MCsamples);
  else if (reweightingStep == 4) step4_closure(base_path, DY_MCsamples, GJets_MCsamples);
  else std::cout<< "You launched this macro without a correct reweighting step: "<< reweightingStep << std::endl;
}
