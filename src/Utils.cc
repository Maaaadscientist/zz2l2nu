#include <Utils.h>

#include <cmath>
#include <fstream>
#include <utility>

#include <TFile.h>

#include <Logger.h>
#include <GenWeight.h>

#define PI 3.141592654

namespace utils {

double deltaPhi(TLorentzVector const &v1, TLorentzVector const &v2)
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


bool PassVbfCuts(std::vector<Jet> const &selJets, TLorentzVector const &boson) {
  if(selJets.size()>=2){
    float etamin=0., etamax=0;
    if(selJets[0].p4.Eta()>selJets[1].p4.Eta()) {etamax = selJets[0].p4.Eta(); etamin = selJets[1].p4.Eta();}
    else {etamax = selJets[1].p4.Eta(); etamin = selJets[0].p4.Eta();}
    bool centralJetVeto = true;
    if(selJets.size()>2){
      for(unsigned int i = 2 ; i < selJets.size() ; i++){
        if(selJets[i].p4.Eta()>etamin && selJets[i].p4.Eta()<etamax) centralJetVeto = false;
      }
    }
    bool centralBoson = (boson.Eta()>etamin && boson.Eta()<etamax);
    bool passDeltaEta = (fabs(selJets[0].p4.Eta() - selJets[1].p4.Eta())>4);
    bool passMjj = ((selJets[0].p4+selJets[1].p4).M()>500);
    if(centralJetVeto && centralBoson && passDeltaEta && passMjj) return true;
  }
  return false;
}

double photon_rhoCorrectedIso(double pfIso, double rho, double sceta,
                              TString const &isoType) {
  return TMath::Max(pfIso - rho * photonID_effArea(sceta, isoType), 0.);
}

double photonID_effArea(double sceta, TString const &isoType) {
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

/*
bool passMetFilter(ULong64_t TrigMET, std::vector<std::pair<int, int> > & listMETFilter, bool isMC){
  //from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#ICHEP_2016
  std::vector<int> trigMET_position_all{14,  // primary vertex filter 
                                        8,  // beam halo filter, Tight cut version.
                                        9,  // beam halo filter, superTight cut version.
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
*/

bool file_exist(std::string const &name) {
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

std::map<double, std::pair<double, double>> TH1toMap(
    std::string const &fileName, std::string const &histoName) {
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

void loadInstrMETWeights(
    bool weight_NVtx_exist, bool weight_Pt_exist, bool weight_Mass_exist,
    std::map<TString, std::map<double, std::pair<double, double>>> &NVtxWeight_map,
    std::map<TString, std::map<double, std::pair<double, double>>> &PtWeight_map,
    std::map<TString, TH1D*> &LineshapeMassWeight_map,
    std::string const &weightFileType, std::string const &base_path,
    std::vector<std::string> const &v_jetCat) {
  if(weight_NVtx_exist){
    LOG_DEBUG << "NVtx weight file has been found! Some histo (called 'After_eeR' and 'After_mumuR') will have the NVtx reweighting applied :)";
    NVtxWeight_map["_ee"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root", "WeightHisto_ee_AllBins");
    NVtxWeight_map["_mumu"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root", "WeightHisto_mumu_AllBins");
    NVtxWeight_map["_ll"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root", "WeightHisto_ll_AllBins");
    if(weight_Pt_exist){
      LOG_DEBUG << "Pt weight file has also been found! Some histo (called 'AfterPtR') will have both reweighting applied :)";
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
      LOG_DEBUG << "Lineshape mass file has been found! Some histo (named 'andMassivePhoton') will have the lineshape applied :)";
      TFile *f_weight_lineshape = TFile::Open((TString) base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_lineshape_mass.root");
      LineshapeMassWeight_map["_ee"] = (TH1D*) ((TH1D*) f_weight_lineshape->Get("WeightHisto_ee_AllBins"))->Clone();
      LineshapeMassWeight_map["_mumu"] = (TH1D*) ((TH1D*) f_weight_lineshape->Get("WeightHisto_mumu_AllBins"))->Clone();
      LineshapeMassWeight_map["_ll"] = (TH1D*) ((TH1D*) f_weight_lineshape->Get("WeightHisto_ll_AllBins"))->Clone();
    }
  }
}


double getTheoryUncertainties(GenWeight const &genWeight,
                              std::string_view syst) {
  if (syst == "QCDscale_up")
    return genWeight.EnvelopeMEScale(GenWeight::Var::Up);
  else if (syst == "QCDscale_down")
    return genWeight.EnvelopeMEScale(GenWeight::Var::Down);
  else if (syst == "alphaS_up")
    return genWeight.RelWeightAlphaS(GenWeight::Var::Up);
  else if (syst == "alphaS_down")
    return genWeight.RelWeightAlphaS(GenWeight::Var::Down);
  else
    return 1.;
}

}  // namespace utils

