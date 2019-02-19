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
#include "CMS_lumi.C"
#include <sstream>
#include <iomanip>

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
      std::ostringstream out;
          out << std::setprecision(n) << a_value;
              return out.str();
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

void make_axis(TAxis* & xaxis, TAxis* & yaxis, int fontType, int pixelFontSize){
  yaxis->SetLabelFont(fontType);
  yaxis->SetLabelSize(pixelFontSize);
  //yaxis->SetNdivisions(505);
  yaxis->SetTitleOffset(2.1); //1.4
  yaxis->SetTitleFont(fontType);
  yaxis->SetTitleSize(pixelFontSize);

}


void draw(TH1F* h_nominal, std::vector<TH1F*> h_syst, TString theHistoName, TString outDir, FILE* yields){

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


  TCanvas *c = new TCanvas("c","canvas",W,H);
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

  c->cd();
  TPad *pad2 =new TPad("down","down",0,0,1,0.27);
 // pad2->SetNumber(2);
  pad2->SetLeftMargin( L/W );
  pad2->SetRightMargin( R/W );
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin( B/H );
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->Draw();

  c->cd();
  pad->cd();

  CMS_lumi(c, 16, 0, true);
  //TLegend *t = new TLegend(0.69,0.66,0.96,0.91);
  TLegend *t = new TLegend(0.38,0.70,0.95,0.90);
  t->SetLineColor(0);
  t->SetBorderSize(1);
  t->SetNColumns(3);
  t->AddEntry(h_nominal, "Nominal", "l");

  h_nominal->SetMarkerColor(kBlack);
  h_nominal->SetLineColor(kBlack);


  //End of the setup of the canvas
  c->cd();
  pad->cd();

  h_nominal->Draw("E1");
  h_nominal->SetTitle("");
  TString textToPrint = theHistoName+"\n";
  std::cout << theHistoName << std::endl;
  textToPrint += "Nominal : \t\t" + to_string_with_precision(h_nominal->Integral(), 4) + "\n";
  std::cout<< "Nominal : \t\t" << h_nominal->Integral() << std::endl; 

  TString name;
  std::string legendName;
  std::string lepCat = (theHistoName.Contains("_ee")) ? "_ee": "_mumu";
  int offset = (theHistoName.Contains("_ee")) ? 4: 6;
  gStyle->SetPalette(1);
  int h_syst_color = 0;
  for(const auto h_s: h_syst){
    name = h_s->GetName();
    legendName = std::string(name);
    std::size_t posBegin = legendName.find(lepCat)+offset;
    if(name.Contains("_up")){
      std::size_t posEnd = legendName.find("_up");
      std::string shortLegend = legendName.substr(posBegin, posEnd-posBegin);
      t->AddEntry(h_s, shortLegend.c_str(), "f");
    }
    if(name.Contains("ZGamma")) h_syst_color = 432; //kCyan
    else if(name.Contains("WGamma")) h_syst_color = 616; //kMagenta
    else if(name.Contains("WJets"))  h_syst_color = 416; //kGreen
 
    if(name.Contains("QCDscale")) h_syst_color += 0;
    else if(name.Contains("alphaS")) h_syst_color += 1;
    else if(name.Contains("pdf")) h_syst_color += 2;

    h_s->SetLineColor(h_syst_color);
    h_s->Draw("HIST same");
    textToPrint += legendName.substr(posBegin) + " : \t\t" + to_string_with_precision(h_s->Integral(), 4) + "\t (" + to_string_with_precision(100.*(h_s->Integral()-h_nominal->Integral())/h_nominal->Integral(), 4) + "%)\n";
    std::cout << legendName.substr(posBegin) << " : \t\t" << h_s->Integral() << " (" << 100.*(h_s->Integral()-h_nominal->Integral())/h_nominal->Integral() << "%)" << std::endl;
  }
  textToPrint += "--------------------------------------\n";
  fprintf(yields, "%s\n", textToPrint.Data());

  t->Draw();
  TAxis* xaxis = new TAxis();
  TAxis* yaxis = new TAxis();
  h_nominal->SetMaximum(h_nominal->GetMaximum()*3);
  xaxis = h_nominal->GetXaxis();
  yaxis = h_nominal->GetYaxis();
  make_axis(xaxis, yaxis, fontType, pixelFontSize);
  gPad->RedrawAxis();

  TString theLeptonCategoryText = "", theJetCategoryText = "";
  if(theHistoName.Contains("ee")) theLeptonCategoryText = "ee";
  else if(theHistoName.Contains("mumu")) theLeptonCategoryText = "#mu#mu";

  if(theHistoName.Contains("eq0jet")) theJetCategoryText = "0 jet";
  if(theHistoName.Contains("geq1jet")) theJetCategoryText = "#geq1 jet";
  if(theHistoName.Contains("vbf")) theJetCategoryText = "VBF";

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

  c->cd();
  pad2->cd();
  for(const auto h_s: h_syst){
    
  TH1F *ratio = (TH1F*) h_nominal->Clone("ratio");
    ratio->Sumw2();
    ratio->Divide(h_nominal, h_s, 1,1);
    name = h_s->GetName();
    if(name.Contains("ZGamma")) h_syst_color = 432; //kCyan
    else if(name.Contains("WGamma")) h_syst_color = 616; //kMagenta
    else if(name.Contains("WJets"))  h_syst_color = 416; //kGreen
 
    if(name.Contains("QCDscale")) h_syst_color += 0;
    else if(name.Contains("alphaS")) h_syst_color += 1;
    else if(name.Contains("pdf")) h_syst_color += 2;

    ratio->SetLineColor(h_syst_color);

    
    //ratio->SetMaximum(1.3);
    //ratio->SetMinimum(0.7);
    ratio->SetMaximum(3.0);
    ratio->SetMinimum(0.);
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


    ratio->Draw("HIST:same");
  }
  c->Print(outDir+theHistoName+".png");
  c->Print(outDir+theHistoName+".root");
  c->Print(outDir+theHistoName+".pdf");

  delete c;
}


void InstrMET_syst_study(std::string suffix){
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(kTRUE); //To ensure that all histograms are created with the sum of weights
  TH2::SetDefaultSumw2(kTRUE); //To ensure that all histograms are created with the sum of weights

  //Cleaning
  TString base_path = std::string(getenv("HZZ2L2NU_BASE")) + "/";
  TString fileDirectory= base_path+"OUTPUTS/"+suffix+"/MERGED";
  TString outDir = base_path+"OUTPUTS/"+suffix+"/PLOTS/SYST_InstrMET/";

  TFile* file = TFile::Open(fileDirectory+"/outputHZZ_InstrMET.root");

  //Check syst that were ran on for Instr.MET and compute syst:
  const char *cmd = ("ls "+std::string(getenv("HZZ2L2NU_BASE")) + "/OUTPUTS/"+suffix+"/OUTPUTS/outputPhotonDatadriven_*{up,down}_0.root | rev | cut -d_ -f2,3 | cut -d_ -f1,2 | rev |sort -u").c_str();
  std::vector<std::string> systList = exec(cmd);

  std::vector<TString> process = {"_WJets", "_WGamma", "_ZGamma"};
  std::vector<TString> jetCat = {"_eq0jets","_geq1jets","_vbf"};
  std::vector<TString> lepCat = {"_ee","_mumu"};
  system("mkdir -p "+base_path+"OUTPUTS/"+suffix+"/PLOTS");
  system("mkdir -p "+base_path+"OUTPUTS/"+suffix+"/PLOTS/SYST_InstrMET");
  system("rm -rf "+base_path+"OUTPUTS/"+suffix+"/PLOTS/SYST_InstrMET/*");

  FILE* yields = fopen(outDir+"yields.txt","w");
  std::vector<TH1F*> h_syst;
  for(const auto jet: jetCat){
    for(const auto lep: lepCat){
      TH1F* h_nominal = (TH1F*) file->Get("mT_final"+jet+lep);
      h_syst.clear();
      for(const auto proc: process){
        for(const auto syst: systList){
          h_syst.push_back( (TH1F*) file->Get("mT_final"+jet+lep+proc+"_"+syst));
        }
      }
      draw(h_nominal, h_syst, "mT_final"+jet+lep, outDir, yields);
    }
  }
  file->Close();

}
