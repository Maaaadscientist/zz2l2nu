#define MeanWeightsComputation_cxx
#include "MeanWeightsComputation.h"
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <string>

#define DO_RESAMPLING_VALIDATION 0

void MeanWeightsComputation::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;

  TFile *fout = new TFile("jackknife_validation_2017.root","recreate");
  TFile *fWeights = new TFile("meanWeights_2017.root","recreate");

  std::vector<TH1F*> lambdas;
  std::vector<TH1F*> mean_weights;
  std::vector<TH1F*> lambdas_perBin;
  std::vector<TH1F*> mean_weights_perBin;
  std::vector<TH1F*> lambdas_tot;
  std::vector<TH1F*> mean_weights_tot;
  TString jet_cat_name[3] = {"eq0jets","geq1jets","vbf"};
  for (int iJetCat = 0 ; iJetCat < 3 ; iJetCat++) {
    lambdas_tot.push_back(new TH1F("lambdas_tot_"+jet_cat_name[iJetCat],"#lambda, combined", n_mT_binning_-1, mT_binning_));
    mean_weights_tot.push_back(new TH1F("mean_weights_tot_"+jet_cat_name[iJetCat],"mean weights, combined", n_mT_binning_-1, mT_binning_));
  }

  // Alternative histograms for validation, using the jackknife method.
  TH2F* lambdas_tot_jackknife = new TH2F("lambdas_tot_jackknife","#lambda, combined", n_mT_binning_-1, mT_binning_, 3000, 0, 300);
  TH2F* mean_weights_tot_jackknife = new TH2F("mean_weights_tot_jackknife","mean weights, combined", n_mT_binning_-1, mT_binning_, 2000, 0, 2);

  std::map<std::pair<int, int>, float> mean;
  std::map<std::pair<int, int>, float> variance;
  std::map<int, float> mean_tot[3];  // 3 objects, 1 for each jet category
  std::map<int, float> variance_tot[3];
  for (int i = 0 ; i < int(pT_thresholds_.size()) ; i++) {
    if (i != int(pT_thresholds_.size())-1) {
      TString name_lambdas = TString("lambdas_"+std::to_string(i));
      TString title_lambdas = TString("#lambda, "+std::to_string(pT_thresholds_[i])+ " < pT < "+std::to_string(pT_thresholds_[i+1]));
      TString name_weights = TString("mean_weights_"+std::to_string(i));
      TString title_weights = TString("mean weights, "+std::to_string(pT_thresholds_[i])+ " < pT < "+std::to_string(pT_thresholds_[i+1]));
      lambdas.push_back(new TH1F(name_lambdas, title_lambdas, n_mT_binning_-1, mT_binning_));
      mean_weights.push_back(new TH1F(name_weights, title_weights, n_mT_binning_-1, mT_binning_));
    }
    for (int j = 0 ; j < int(mT_thresholds_.size()) ; j++) {
      mean[std::make_pair(i,j)] = 0;
      variance[std::make_pair(i,j)] = 0;
      for (int k = 0 ; k < 3 ; k++) {
        mean_tot[k][j] = 0;
        variance_tot[k][j] = 0;
      }
    }
  }
  for (int i = 0 ; i < int(mT_thresholds_.size()) ; i++) {
    if (i != int(mT_thresholds_.size())-1) {
      TString name_lambdas = TString("lambdas_jackknife_"+std::to_string(i));
      TString title_lambdas = TString("#lambda, 100 jackknife samples, "+std::to_string(mT_thresholds_[i])+ " < mT < "+std::to_string(mT_thresholds_[i+1])+";#lambda;Number of jackknife replicas");
      TString name_weights = TString("mean_weights_jackknife_"+std::to_string(i));
      TString title_weights = TString("mean weights, 100 jackknife samples, "+std::to_string(mT_thresholds_[i])+ " < mT < "+std::to_string(mT_thresholds_[i+1])+";Mean weight;Number of jackknife replicas");
      lambdas_perBin.push_back(new TH1F(name_lambdas, title_lambdas, 50, 0, 0));
      mean_weights_perBin.push_back(new TH1F(name_weights, title_weights, 50, 0, 0));
    }
  }

  std::vector<std::map<int, float>> mean_tot_jackknife;
  std::vector<std::map<int, float>> variance_tot_jackknife;
  for (int i = 0 ; i < 100 ; i++) {
    mean_tot_jackknife.push_back(std::map<int, float>());
    variance_tot_jackknife.push_back(std::map<int, float>());
    for (int j = 0 ; j < int(mT_thresholds_.size()) ; j++) {
      mean_tot_jackknife[i][j] = 0;
      variance_tot_jackknife[i][j] = 0;
    }
  }

  // First loop, to get the estimation
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry % 100000 == 0) std::cout << "Event " << jentry << std::endl;

    if (ptmiss < 125) continue;

    std::pair<int, int> bins;
    bins = find_thresholds_binning(photon_pt, mT);
    //std::cout << "pT = " << photon_pt << ", bin pT = " << bins.first << ", mT = " << mT << ", bin mT = " << bins.second << ", weight = " << photon_reweighting << " * " << trigger_weight << " = " << photon_reweighting*trigger_weight << std::endl;
    mean[bins] += photon_reweighting*trigger_weight;
    variance[bins] += pow(photon_reweighting*trigger_weight, 2);
    if (jet_cat > 2) continue;
    mean_tot[jet_cat][bins.second] += photon_reweighting*trigger_weight;
    variance_tot[jet_cat][bins.second] += pow(photon_reweighting*trigger_weight, 2);
  }

  // Second loop, done 100 times, for jackkife method
  if (DO_RESAMPLING_VALIDATION) {
    for (int iJack = 0 ; iJack < 100 ; iJack++) {
      std::cout << "Jackknife method, iteration " << iJack << "..." << std::endl;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        if (jentry % 1000000 == 0) std::cout << "Event " << jentry << std::endl;
        if (jentry % 100 == iJack) continue; // Don't keep 1 in 100 events.

        if (ptmiss < 125) continue;

        std::pair<int, int> bins;
        bins = find_thresholds_binning(photon_pt, mT);
        mean_tot_jackknife[iJack][bins.second] += photon_reweighting*trigger_weight;
        variance_tot_jackknife[iJack][bins.second] += pow(photon_reweighting*trigger_weight, 2);
      }
    }
  }

  for (int i = 0 ; i < int(lambdas.size()) ; i++) {
    for (int j = 0 ; j < int(lambdas[i]->GetNbinsX()) ; j++) {
      //std::cout << "mean[" << i << "," << j << "] = " << mean[std::make_pair(i, j)] << std::endl;
      //std::cout << "variance[" << i << "," << j << "] = " << variance[std::make_pair(i, j)] << std::endl;
      //std::cout << "so bin "<< j+1 << " of lambda[" << i << "] = " << pow(mean[std::make_pair(i,j)],2)/variance[std::make_pair(i,j)] << std::endl;
      if (variance[std::make_pair(i,j)] != 0) {
        lambdas[i]->SetBinContent(j+1, pow(mean[std::make_pair(i,j)],2)/variance[std::make_pair(i,j)]);
        lambdas[i]->SetBinError(j+1, mean[std::make_pair(i,j)]/sqrt(variance[std::make_pair(i,j)]));
      }
      else {
        lambdas[i]->SetBinContent(j+1,0);
      }
      if (mean[std::make_pair(i,j)] != 0) {
        mean_weights[i]->SetBinContent(j+1, variance[std::make_pair(i,j)]/mean[std::make_pair(i,j)]);
        mean_weights[i]->SetBinError(j+1, 0);
      }
      else {
        mean_weights[i]->SetBinContent(j+1,0);
      }
      lambdas[i]->SetBinContent(0,0);
      mean_weights[i]->SetBinContent(0,0);
    }
  }
  for (int iJetCat = 0 ; iJetCat < 3 ; iJetCat++) {
    for (int j = 0 ; j < int(lambdas_tot[iJetCat]->GetNbinsX() + 1) ; j++) {
      if (variance_tot[iJetCat][j] != 0) {
        lambdas_tot[iJetCat]->SetBinContent(j+1, pow(mean_tot[iJetCat][j],2)/variance_tot[iJetCat][j]);
        lambdas_tot[iJetCat]->SetBinError(j+1, mean_tot[iJetCat][j]/sqrt(variance_tot[iJetCat][j]));
      }
      else {
        lambdas_tot[iJetCat]->SetBinContent(j+1,0);
      }
      if (mean_tot[iJetCat][j] != 0) {
        mean_weights_tot[iJetCat]->SetBinContent(j+1, variance_tot[iJetCat][j]/mean_tot[iJetCat][j]);
        mean_weights_tot[iJetCat]->SetBinError(j+1, 0);
      }
      else {
        mean_weights_tot[iJetCat]->SetBinContent(j+1,0);
      }
    }
    lambdas_tot[iJetCat]->SetBinContent(0,0);
    mean_weights_tot[iJetCat]->SetBinContent(0,0);
  }
  for (int iJack = 0 ; iJack < 100 ; iJack++) {
    for (int j = 0 ; j < int(lambdas_tot_jackknife->GetNbinsX() + 1) ; j++) {
      float binLowEdge = lambdas_tot_jackknife->GetXaxis()->GetBinLowEdge(j+1);
      if (variance_tot_jackknife[iJack][j] != 0) {
        lambdas_tot_jackknife->Fill(binLowEdge, pow(mean_tot_jackknife[iJack][j],2)/variance_tot_jackknife[iJack][j]*100./99.);
        lambdas_perBin[j]->Fill(pow(mean_tot_jackknife[iJack][j],2)/variance_tot_jackknife[iJack][j]*100./99.);
      }
      else {
        lambdas_tot_jackknife->Fill(binLowEdge,0);
      }
      if (mean_tot_jackknife[iJack][j] != 0) {
        mean_weights_tot_jackknife->Fill(binLowEdge, variance_tot_jackknife[iJack][j]/mean_tot_jackknife[iJack][j]*99./100.);
        mean_weights_perBin[j]->Fill(variance_tot_jackknife[iJack][j]/mean_tot_jackknife[iJack][j]*99./100.);
      }
      else {
        mean_weights_tot_jackknife->Fill(binLowEdge,0);
      }
    }
    lambdas_tot_jackknife->Fill(-1,0);
    mean_weights_tot_jackknife->Fill(-1,0);
  }

  fout->cd();
  for (auto hist : lambdas) 
    hist->Write();
  for (auto hist : mean_weights)
    hist->Write();
  for (auto hist : lambdas_perBin)
    hist->Write();
  for (auto hist : mean_weights_perBin)
    hist->Write();
  for (int iJetCat = 0 ; iJetCat < 3 ; iJetCat++) {
    lambdas_tot[iJetCat]->Write();
    mean_weights_tot[iJetCat]->Write();
  }
  lambdas_tot_jackknife->Write();
  mean_weights_tot_jackknife->Write();
  fout->Close();

  fWeights->cd();
  for (int iJetCat = 0 ; iJetCat < 3 ; iJetCat++) {
    mean_weights_tot[iJetCat]->Write();
  }
  fWeights->Close();
}

std::pair<int, int> MeanWeightsComputation::find_thresholds_binning(float pt, float mt)
{
  std::vector<float>::iterator it_pt;
  it_pt = lower_bound(pT_thresholds_.begin(), pT_thresholds_.end(), pt);
  int ptBin = it_pt - pT_thresholds_.begin() - 1;
  std::vector<float>::iterator it_mt;
  it_mt = lower_bound(mT_thresholds_.begin(), mT_thresholds_.end(), mt);
  int mtBin = it_mt - mT_thresholds_.begin() - 1;
  return std::make_pair(ptBin, mtBin);
}
