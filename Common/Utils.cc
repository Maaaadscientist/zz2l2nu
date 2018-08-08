#include "Utils.h"
#include <utility>
#include <TFile.h>
#include <fstream>
#include <list>

#define PI 3.141592654

namespace utils
{
  double deltaPhi (TLorentzVector v1, TLorentzVector v2)
  {
    return deltaPhi (v1.Phi(), v2.Phi());
  }

  double deltaPhi (float phi1, float phi2)
  {
    double dPhi = phi1 - phi2; 
    if (dPhi > PI) dPhi -= 2*PI;
    if (dPhi <= -PI) dPhi += 2*PI; //So that dPhi is always between -pi and +pi.
    return dPhi;
  }

  double deltaR(TLorentzVector v1, TLorentzVector v2)
  {
    double dEta = v1.Eta()-v2.Eta();
    double dPhi = deltaPhi(v1,v2);
    return sqrt(pow(dEta,2) + pow(dPhi,2));
  }

  double deltaR(float eta1, float phi1, float eta2, float phi2)
  {
    double dEta = eta1 - eta2; 
    double dPhi = deltaPhi(phi1, phi2);
    return sqrt(pow(dEta,2) + pow(dPhi,2));
  }

  double getPhotonEnergy (double pT, double eta)
  {
    double theta = 2*atan(exp(-eta));
    return fabs(pT/sin(theta));
  }

  bool passVBFcuts(std::vector<TLorentzVectorWithIndex> selJets, TLorentzVector boson)
  {
    if(selJets.size()>=2){
      float etamin=0., etamax=0;
      if(selJets[0].Eta()>selJets[1].Eta()) {etamax = selJets[0].Eta(); etamin = selJets[1].Eta();}
      else {etamax = selJets[1].Eta(); etamin = selJets[0].Eta();}
      bool centralJetVeto = true;
      if(selJets.size()>2){
        for(unsigned int i = 2 ; i < selJets.size() ; i++){
          if(selJets[i].Eta()>etamin && selJets[i].Eta()<etamax) centralJetVeto = false;
        }
      }
      bool centralBoson = (boson.Eta()>etamin && boson.Eta()<etamax);
      bool passDeltaEta = (fabs(selJets[0].Eta() - selJets[1].Eta())>4);
      bool passMjj = ((selJets[0]+selJets[1]).M()>500);
      if(centralJetVeto && centralBoson && passDeltaEta && passMjj) return true;
    }
    return false;
  }

  double photon_rhoCorrectedIso(double pfIso, double rho, double sceta, TString isoType){
    return TMath::Max(pfIso - rho * photonID_effArea(sceta, isoType), 0.);
  }

  double photonID_effArea(double sceta, TString isoType){
    double effArea(1.0);
    if(isoType=="chIso"){
      if(fabs(sceta)<1.0)                           effArea=0.0360;
      else if(fabs(sceta)>1.0 && fabs(sceta)<1.479) effArea=0.0377;
      else if(fabs(sceta)>1.479 && fabs(sceta)<2.0) effArea=0.0306;
      else if(fabs(sceta)>2.0 && fabs(sceta)<2.2)   effArea=0.0283;
      else if(fabs(sceta)>2.2 && fabs(sceta)<2.3)   effArea=0.0254;
      else if(fabs(sceta)>2.3 && fabs(sceta)<2.4)   effArea=0.0217;
      else                                          effArea=0.0167;
    }
    if(isoType=="nhIso"){
      if(fabs(sceta)<1.0)                           effArea=0.0597;
      else if(fabs(sceta)>1.0 && fabs(sceta)<1.479) effArea=0.0807;
      else if(fabs(sceta)>1.479 && fabs(sceta)<2.0) effArea=0.0629;
      else if(fabs(sceta)>2.0 && fabs(sceta)<2.2)   effArea=0.0197;
      else if(fabs(sceta)>2.2 && fabs(sceta)<2.3)   effArea=0.0184;
      else if(fabs(sceta)>2.3 && fabs(sceta)<2.4)   effArea=0.0284;
      else                                          effArea=0.0591;
    }
    if(isoType=="gIso"){
      if(fabs(sceta)<1.0)                           effArea=0.1210;
      else if(fabs(sceta)>1.0 && fabs(sceta)<1.479) effArea=0.1107;
      else if(fabs(sceta)>1.479 && fabs(sceta)<2.0) effArea=0.0699;
      else if(fabs(sceta)>2.0 && fabs(sceta)<2.2)   effArea=0.1056;
      else if(fabs(sceta)>2.2 && fabs(sceta)<2.3)   effArea=0.1457;
      else if(fabs(sceta)>2.3 && fabs(sceta)<2.4)   effArea=0.1719;
      else                                          effArea=0.1998;
    }
    return effArea;
  }

  bool passMetFilter(ULong64_t TrigMET, std::vector<std::pair<int, int> > & listMETFilter, bool isMC){
    //from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#ICHEP_2016
    std::vector<int> trigMET_position_all{14,  // primary vertex filter 
                                          8,  // beam halo filter.
                                          //see below: 3,  // HBHE noise filter 
                                          //see below: 4,  // HBHEiso noise filter 
                                          12,  // ECAL TP filter 
                                          15,  // ee badSC noise filter (not suggested for MC, but should be 100% according to David. To check) 
                                          24,  // badMuon (actually they went badPFMuon) 
                                          25   // badCharged hadron 
                                            };
    std::vector<int> trigMET_position_dataOnly{0, //duplicateMuons
                                               1, //badMuon
                                               2, //noBadMuon (i.e. no duplicate and no bad muons)
                                               3, //for the data, the condition on HBHE noise filter is inversed!
                                               4  //for the data, the condition on HBHEiso noise filter is inversed!
                                                };
    std::vector<int> trigMET_position_MCOnly{3, // HBHE noise filter (is also applied on data, but with the opposite value for the flag!)
                                             4  // HBHEiso noise filter (is also applied on data, but with the opposite value for the flag!)
                                              };

    bool passAllFilters = true;
    for(unsigned int i = 0; i < trigMET_position_all.size(); i++){
      if((TrigMET & (1<<trigMET_position_all[i]))) listMETFilter.push_back(std::make_pair(trigMET_position_all[i], 1));
      else{
        listMETFilter.push_back(std::make_pair(trigMET_position_all[i], 0));
        passAllFilters=false;
      }
    }
    if(!isMC){
      //For the 5 flags on data above, the logic is inversed: we ask the bin to be 0!
      for(unsigned int i = 0; i < trigMET_position_dataOnly.size(); i++){
        if(!(TrigMET & (1<<trigMET_position_dataOnly[i]))) listMETFilter.push_back(std::make_pair(trigMET_position_dataOnly[i], 1));
        else{
          listMETFilter.push_back(std::make_pair(trigMET_position_dataOnly[i], 0));
          passAllFilters=false;
        }
      }
    }
    else{
      for(unsigned int i = 0; i < trigMET_position_MCOnly.size(); i++){
        if((TrigMET & (1<<trigMET_position_MCOnly[i]))) listMETFilter.push_back(std::make_pair(trigMET_position_MCOnly[i], 1));
        else{
          listMETFilter.push_back(std::make_pair(trigMET_position_MCOnly[i], 0));
          passAllFilters=false;
        }
      }

    }
    return passAllFilters;
  }

  bool file_exist(std::string name){
    std::ifstream f(name.c_str());
    return f.good();
  }

  std::map<double, double> TH1toMap(TH1D *h_weight){
    std::map<double, double> myMap;
    int nBinsX = h_weight->GetNbinsX();
    for(int bin_X = 1; bin_X <= nBinsX; bin_X++){
      myMap[h_weight->GetXaxis()->GetBinLowEdge(bin_X)] = h_weight->GetBinContent(bin_X); 
    }
    return myMap;
  }

  std::map<double, std::pair<double, double> > TH1toMap(std::string fileName, std::string histoName){
    std::map<double, std::pair<double, double> > myMap; //label = bin low edge; pair = (bin content; bin error)
    TFile *f_weight = TFile::Open((TString) fileName);
    TH1D *h_weight = (TH1D*) f_weight->Get((TString) histoName);
    int nBinsX = h_weight->GetNbinsX();
    for(int bin_X = 1; bin_X <= nBinsX; bin_X++){
      myMap[h_weight->GetXaxis()->GetBinLowEdge(bin_X)].first = h_weight->GetBinContent(bin_X); 
      myMap[h_weight->GetXaxis()->GetBinLowEdge(bin_X)].second = h_weight->GetBinError(bin_X); 
      
    }
    f_weight->Close();
    return myMap;
  }

  void giveMassToPhoton(TLorentzVector & boson, TH1D* h_weight){
    double mass = 0;
    while(fabs(mass-91)>15){
      mass = h_weight->GetRandom();
    }
    boson.SetPtEtaPhiE(boson.Pt(), boson.Eta(), boson.Phi(), sqrt(pow(mass,2)+pow(boson.P(),2)) );
  }

  void loadInstrMETWeights(bool weight_NVtx_exist, bool weight_Pt_exist, bool weight_Mass_exist, std::map<TString, std::map<double, std::pair<double, double> > > & NVtxWeight_map, std::map<TString, std::map<double, std::pair<double, double> > > & PtWeight_map, std::map<TString, TH1D*> & LineshapeMassWeight_map, std::string weightFileType, std::string base_path, std::vector<std::string> v_jetCat){
    if(weight_NVtx_exist){
      std::cout << "NVtx weight file has been found! Some histo (called 'After_eeR' and 'After_mumuR') will have the NVtx reweighting applied :)" << std::endl;
      NVtxWeight_map["_ee"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root", "WeightHisto_ee_AllBins");
      NVtxWeight_map["_mumu"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root", "WeightHisto_mumu_AllBins");
      NVtxWeight_map["_ll"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root", "WeightHisto_ll_AllBins");
      if(weight_Pt_exist){
        std::cout << "Pt weight file has also been found! Some histo (called 'AfterPtR') will have both reweighting applied :)" << std::endl;
        PtWeight_map["_ee"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto_ee_AllBins");
        PtWeight_map["_mumu"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto_mumu_AllBins");
        PtWeight_map["_ll"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto_ll_AllBins");
        for(unsigned int i =0; i < v_jetCat.size(); i++){
          PtWeight_map["_ee"+v_jetCat[i]] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto"+v_jetCat[i]+"_ee_AllBins");
          PtWeight_map["_mumu"+v_jetCat[i]] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto"+v_jetCat[i]+"_mumu_AllBins");
          PtWeight_map["_ll"+v_jetCat[i]] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto"+v_jetCat[i]+"_ll_AllBins");
        }
      }
      if(weight_Mass_exist){
        std::cout << "Lineshape mass file has been found! Some histo (named 'andMassivePhoton') will have the lineshape applied :)" << std::endl;
        TFile *f_weight_lineshape = TFile::Open((TString) base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_lineshape_mass.root");
        LineshapeMassWeight_map["_ee"] = (TH1D*) ((TH1D*) f_weight_lineshape->Get("WeightHisto_ee_AllBins"))->Clone();
        LineshapeMassWeight_map["_mumu"] = (TH1D*) ((TH1D*) f_weight_lineshape->Get("WeightHisto_mumu_AllBins"))->Clone();
        LineshapeMassWeight_map["_ll"] = (TH1D*) ((TH1D*) f_weight_lineshape->Get("WeightHisto_ll_AllBins"))->Clone();
      }
    }
  }

  double getTheoryUncertainties(std::vector<double> *EvtWeights, TString syst){
    // if(syst == "pdf_up") return getPdfUncertainty(EvtWeights, true); // Now done in the main code.
    // else if(syst == "pdf_down") return getPdfUncertainty(EvtWeights, false);
    if(syst == "QCDscale_up") return getQCDScaleUncertainty(EvtWeights, true);
    else if(syst == "QCDscale_down") return getQCDScaleUncertainty(EvtWeights, false);
    else if(syst == "alphaS_up") return getAlphaUncertainty (EvtWeights, true);
    else if(syst == "alphaS_down") return getAlphaUncertainty (EvtWeights, false);
    else return 1.;
  }

  double getPdfUncertainty(std::vector<double> *EvtWeights, bool isUp){ //DEPRECATED. This function is not used anymore, only kept for comparison with this old method.
    if(EvtWeights->size() < 110) throw std::out_of_range("Vector of weights not filled properly."); //This happened randomly for some events in ZZ2l2v for 2016 MC.
    double squaredSum = 0.;
    for(int i = 10 ; i < 110 ; i++){ // Correspond to the PDF replicas giving different weights.
      squaredSum += (1.*EvtWeights->at(i)/EvtWeights->at(0) - 1.) * (1.*EvtWeights->at(i)/EvtWeights->at(0) - 1.);
    }
    double pdfFinalWeight = 1;
    if(isUp) pdfFinalWeight = 1. + sqrt(squaredSum/99.); // Standard deviation of the distribution. N-1 = 99 (we have 100 entries)
    else pdfFinalWeight = 1. - sqrt(squaredSum/99.);
    return pdfFinalWeight;
  }

  double getQCDScaleUncertainty(std::vector<double> *EvtWeights, bool isUp){
    std::vector<int> indexes = {2, 3, 4, 5, 7, 9}; // Correspond to id 1002, 1003,..., 1009, which account for the variations of mu_R and/or mu_F by a factor 0.5, 1 or 2. The 2 cases not considered are mu_R = 2 / mu_F = 0.5, and the reverse.
    if(EvtWeights->size() < indexes.size()) throw std::out_of_range("Vector of weights not filled properly."); //This happened randomly for some events in ZZ2l2v for 2016 MC.
    std::list<double> QCDScaleWeights;
    for(int i = 0 ; i < indexes.size() ; i++) QCDScaleWeights.push_back(EvtWeights->at(indexes[i]));
    double QCDFinalWeight = 1.;
    if(isUp) QCDFinalWeight = *std::max_element(QCDScaleWeights.begin(),QCDScaleWeights.end())/EvtWeights->at(0)/1.; // We take (conservatively) the biggest variation for the scale up
    else QCDFinalWeight = *std::min_element(QCDScaleWeights.begin(),QCDScaleWeights.end())/EvtWeights->at(0)/1.;
    return QCDFinalWeight;
  }

  double getAlphaUncertainty(std::vector<double> *EvtWeights, bool isUp){
    if(EvtWeights->size() < 112) throw std::out_of_range("Vector of weights not filled properly."); //This happened randomly for some events in ZZ2l2v for 2016 MC.
    double alphaWeight = 1.;
    double alphaUnc = fabs(0.5*(EvtWeights->at(110)-EvtWeights->at(111))/EvtWeights->at(0)); // Method used to symmetrize the uncertainty. There was a mysterious factor sqrt(0.75) in the old framework that I simply removed.
    if(isUp) alphaWeight = 1. + alphaUnc;
    else alphaWeight = 1. - alphaUnc;
    return alphaWeight;
  }

}
