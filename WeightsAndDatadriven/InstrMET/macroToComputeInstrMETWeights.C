#include <TH1F.h>
#include <TH2F.h>
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

//You can tune easily the options here
//Reweighting options
#define DO_NVTX_VS_PT_REWEIGHTING true // if true, do the NVt reweighting by bin of trigger threshold (pT)

//General options
#define VERBOSE true
#define HIDE_WARNING true
#define DEBUG_HISTOS true



void step1_weight_NVtx_vs_pt(TString base_path, TFile *f_HZZ) {
  if(VERBOSE) std::cout<< "Launching step1 of the computation of the weights for Instr. MET: NVtx vs Pt" << std::endl; 
  // Initialization: opening and creating ROOT files
  TString InstrMET_noWeight_path = base_path + "OUTPUTS/InstrMET_dataOnly_forNvtxWeights/merged/"; // Path of photon data without any reweighting for Instr. MET

  TFile *f_InstrMET = TFile::Open(InstrMET_noWeight_path+"outputInstrMET_Data.root");
  TFile *f_output = new TFile(base_path+"WeightsAndDatadriven/InstrMET/InstrMET_weight_NVtx.root","RECREATE");

  // Main Loop of Step1
  std::vector<TString> jetCat = {"_eq0jets","_geq1jets","_vbf"}; jetCat.push_back(""); //Adding all jet categories
  std::vector<TString> lepCat = {"_ee","_mumu", "_ll"};
  std::vector<std::vector<TH1F*> > HZZ_nvtx(jetCat.size());

  //If you use 2D reweighting (DO_NVTX_VS_PT_REWEIGHTING == true)
  std::vector<TH2F*> InstrMET_zpt_vs_nvtx;
  std::vector<std::vector<TH2F*> > results_2D(jetCat.size());

  //If you use 1D reweighting (DO_NVTX_VS_PT_REWEIGHTING == false)
  std::vector<TH1F*> InstrMET_nvtx;
  std::vector<TH1F*> results_1D;

  double weight, error;

  for(unsigned int i = 0; i < jetCat.size(); i++){
    if(!DO_NVTX_VS_PT_REWEIGHTING) {
      InstrMET_nvtx.push_back( (TH1F*) f_InstrMET->Get("reco-vtx_MET125_InstrMET_reweighting"+jetCat[i]+"_gamma"));
      InstrMET_nvtx[i]->Scale(1./InstrMET_nvtx[i]->Integral());
      if(DEBUG_HISTOS) InstrMET_nvtx[i]->Write();
    }
    if(DO_NVTX_VS_PT_REWEIGHTING){
      InstrMET_zpt_vs_nvtx.push_back( (TH2F*) f_InstrMET->Get("nvtxvsBosonPt_2D_MET125_InstrMET_reweighting"+jetCat[i]+"_gamma"));
      if(DEBUG_HISTOS) InstrMET_zpt_vs_nvtx[i]->Write(); //Not normalized, just to see the raw distribution
    }
    for(unsigned int j = 0; j < lepCat.size(); j++){

      if(DO_NVTX_VS_PT_REWEIGHTING){
        std::cout << "In jetCat " << jetCat[i] << " and lepCat " << lepCat[j] << std::endl;
        HZZ_nvtx[i].push_back( (TH1F*) f_HZZ->Get("reco-vtx_MET125_InstrMET_reweighting"+jetCat[i]+lepCat[j]));
        HZZ_nvtx[i][j]->Scale(1./HZZ_nvtx[i][j]->Integral());
        results_2D[i].push_back( (TH2F*) InstrMET_zpt_vs_nvtx[i]->Clone("InstrMET_weight_zpt_vs_nvtx"+jetCat[i]+lepCat[j]));

        for(int pt_bin = 0; pt_bin <= results_2D[i][j]->GetNbinsX(); pt_bin++){ //Starting at 1 and going at < (instead of <=) to avoid under- and over-flow bins.
          TH1D* temp_projection = InstrMET_zpt_vs_nvtx[i]->ProjectionY("zpt_vs_nvtx_InstrMET_reweighting"+jetCat[i]+"_gamma_bin"+std::to_string(pt_bin), pt_bin, pt_bin);
          temp_projection->Scale(1./temp_projection->Integral());
          if(DEBUG_HISTOS && j ==0) temp_projection->Write(); //Normalized and used to compute weights. The j ==0 is there to only write this plot once, since it's the same for ee or mumu (since it's a gamma plot)
          TH1D* temp_WeightHisto = ((TH1D*) HZZ_nvtx[i][j]->Clone("WeightHisto"+jetCat[i]+lepCat[j]+"_bin"+std::to_string(pt_bin)));
          if(temp_WeightHisto->GetEntries() != 0 ) temp_WeightHisto->Divide(temp_projection);
          if(DEBUG_HISTOS) temp_WeightHisto->Write(); //The weights histo - 1 histo per pt bin
          for(int nvtx_bin = 0; nvtx_bin <= results_2D[i][j]->GetNbinsY(); nvtx_bin++){
            if(InstrMET_zpt_vs_nvtx[i]->GetBinContent(pt_bin, nvtx_bin) > 0){
              weight=0;
              error=0;
              weight=1.*temp_WeightHisto->GetBinContent(nvtx_bin);
              error=1.*temp_WeightHisto->GetBinError(nvtx_bin);
              if(HZZ_nvtx[i][j]->GetBinContent(nvtx_bin) == 0) std::cout<< " --- /!\\ ERROR /!\\ --- \n The binning chosen for the Nvtx is given 0 value for Z+jets while it has a positive value for Photon+jets (here for Nvtx bin " << nvtx_bin <<". Please change the binning!" << std::endl;
            }
            else{
              weight=0;
              error=0;
            }
            results_2D[i][j]->SetBinContent(pt_bin, nvtx_bin, weight);
            results_2D[i][j]->SetBinError(pt_bin, nvtx_bin, error);
          }
        }
        results_2D[i][j]->Write();
      }
      else{
        HZZ_nvtx[i].push_back( (TH1F*) f_HZZ->Get("reco-vtx_MET125_InstrMET_reweighting"+jetCat[i]+lepCat[j]));
        HZZ_nvtx[i][j]->Scale(1./HZZ_nvtx[i][j]->Integral());
        results_1D.push_back( (TH1F*) HZZ_nvtx[i][j]->Clone("WeightHisto"+jetCat[i]+lepCat[j]+"_AllBins"));
        if(results_1D.back()->GetEntries() != 0 ) results_1D.back()->Divide(InstrMET_nvtx[i]);
        results_1D.back()->Write();

      }
      if(DEBUG_HISTOS) HZZ_nvtx[i][j]->Write(); //Normalized and used to compute weights


    }

  }
  if(VERBOSE) std::cout<< "Step1 of reweighting done. Weights are available under: InstrMET_weight_NVtx.root" << std::endl; 

}

void step2_weight_pt(TString base_path, TFile *f_HZZ) {
  if(VERBOSE) std::cout<< "Launching step2 of the computation of the weights for Instr. MET: Pt" << std::endl; 

  // Initialization: opening and creating ROOT files
  TString InstrMET_NVtxWeight_path = base_path + "OUTPUTS/InstrMET_dataOnly_forPtWeights/merged/"; // Path of photon data already reweighted for the # of vertices

  TFile *f_InstrMET = TFile::Open(InstrMET_NVtxWeight_path+"outputInstrMET_Data.root");
  TFile *f_output = new TFile(base_path+"WeightsAndDatadriven/InstrMET/InstrMET_weight_pt.root","RECREATE");

  std::vector<TString> jetCat = {"_eq0jets", "_geq1jets", "_vbf", ""};
  std::vector<TString> lepCat = {"_ee", "_mumu", "_ll"};

  std::vector<std::vector<TH1F*> > HZZ_nvtx(jetCat.size());
  std::vector<std::vector<TH1F*> > InstrMET_pt(jetCat.size());
  std::vector<TH1F*> results_1D;
  for(unsigned int i = 0; i < jetCat.size(); i++){
    for(unsigned int j = 0; j < lepCat.size(); j++){
      InstrMET_pt[i].push_back( (TH1F*) f_InstrMET->Get("pT_Boson_MET125_InstrMET_reweightingAfter"+lepCat[j]+jetCat[i]));
      //InstrMET_pt[i][j]->Scale(1./InstrMET_pt[i][j]->Integral()); //In order to have a good renormallization between gamma and dilepton data, we don't scale their weights
      if(DEBUG_HISTOS) InstrMET_pt[i][j]->Write();
      HZZ_nvtx[i].push_back( (TH1F*) f_HZZ->Get("pT_Boson_MET125_InstrMET_reweighting"+jetCat[i]+lepCat[j]));
      //HZZ_nvtx[i][j]->Scale(1./HZZ_nvtx[i][j]->Integral()); //In order to have a good renormallization between gamma and dilepton data, we don't scale their weights
      results_1D.push_back( (TH1F*) HZZ_nvtx[i][j]->Clone("WeightHisto"+jetCat[i]+lepCat[j]+"_AllBins"));
      if(results_1D.back()->GetEntries() != 0 ) results_1D.back()->Divide(InstrMET_pt[i][j]);
      results_1D.back()->Write();
      if(DEBUG_HISTOS) HZZ_nvtx[i][j]->Write(); //Normalized and used to compute weights
    }
  }
  //HZZ_zpt.push_back( (TH1F*) f_HZZ->Get("pT_Boson_InstrMET_reweighting__"+jetCat[i]+"__"+lepCat[j]));
  if(VERBOSE) std::cout<< "Step2 of reweighting done. Weights are available under: InstrMET_weight_pt.root" << std::endl;

}

void step3_weight_lineshape(TString base_path, TFile *f_HZZ) {
  if(VERBOSE) std::cout<< "Launching step3 of the computation of the weights for Instr. MET: mass lineshape" << std::endl;

  TFile *f_output = new TFile(base_path+"WeightsAndDatadriven/InstrMET/InstrMET_lineshape_mass.root","RECREATE");

  std::vector<TH1F*> results_1D;

  std::vector<TString> jetCat = {"_eq0jets", "_geq1jets", "_vbf", ""};
  std::vector<TString> lepCat = {"_ee", "_mumu", "_ll"};

  for(unsigned int i = 0; i < jetCat.size(); i++){
    for(unsigned int j = 0; j < lepCat.size(); j++){
      results_1D.push_back( (TH1F*) ((TH1F*) f_HZZ->Get("M_Boson_MET125_InstrMET_reweighting"+jetCat[i]+lepCat[j]))->Clone("WeightHisto"+jetCat[i]+lepCat[j]+"_AllBins") );
      results_1D.back()->Write();
    }
  }
  if(VERBOSE) std::cout<< "Step3 of reweighting done. Weights are available under: InstrMET_lineshape_mass.root" << std::endl;
}

void macroToComputeInstrMETWeights(int reweightingStep) {
  gROOT->SetBatch();
  TH1::SetDefaultSumw2(kTRUE); //To ensure that all histograms are created with the sum of weights
  if(HIDE_WARNING) gErrorIgnoreLevel=kError;

  //Define all paths
  TString base_path = std::string(getenv("HZZ2L2NU_BASE")) + "/";
  //TString HZZ_path = base_path + "OUTPUTS/computeInstrMET_DiLeptonData/MERGED/"; // The path of the DiLepton data used for reweighting
  TString HZZ_path = base_path + "OUTPUTS/HZZ_dataOnly_with2DHisto/merged/"; // The path of the DiLepton data used for reweighting

  //Open and create ROOT files
  TFile *f_HZZ = TFile::Open(HZZ_path+"outputHZZ_Data.root");

  if(reweightingStep == 1) step1_weight_NVtx_vs_pt(base_path, f_HZZ);
  else if (reweightingStep == 2) step2_weight_pt(base_path, f_HZZ);
  else if (reweightingStep == 3) step3_weight_lineshape(base_path, f_HZZ);
  else std::cout<< "You launched this macro without a correct reweighting step: "<< reweightingStep << std::endl;
}
