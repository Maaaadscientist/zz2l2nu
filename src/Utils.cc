#include <Utils.h>

#include <cmath>
#include <fstream>
#include <utility>

#include <TFile.h>

#include <FileInPath.h>
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

std::map<std::pair<double, double>, std::pair<double, double>> TH2toMap(
    std::string const &fileName, std::string const &histoName) {
  std::map<std::pair<double,double>, std::pair<double, double> > myMap; //label = bin low edge; pair = (bin content; bin error)
  TFile *f_weight = TFile::Open((TString) fileName);
  TH2D *h_weight = (TH2D*) f_weight->Get((TString) histoName);
  int nBinsX = h_weight->GetNbinsX();
  int nBinsY = h_weight->GetNbinsY();
  for(int bin_X = 1; bin_X <= nBinsX; bin_X++){
    for(int bin_Y = 1; bin_Y <= nBinsY; bin_Y++){
      myMap[std::make_pair(h_weight->GetYaxis()->GetBinLowEdge(bin_Y),h_weight->GetXaxis()->GetBinLowEdge(bin_X))].first = h_weight->GetBinContent(bin_X,bin_Y); //to have (vtx,pT) to allow for comparisons with pairs. Could be updated later.
      myMap[std::make_pair(h_weight->GetYaxis()->GetBinLowEdge(bin_Y),h_weight->GetXaxis()->GetBinLowEdge(bin_X))].second = h_weight->GetBinError(bin_X,bin_Y);
    }
  }
  f_weight->Close();
  return myMap;
}

void giveMassToPhoton(TLorentzVector & boson, TH1* h_weight){
  double mass = 0;
  while(fabs(mass-91)>15){
    mass = h_weight->GetRandom();
  }
  boson.SetPtEtaPhiE(boson.Pt(), boson.Eta(), boson.Phi(), sqrt(pow(mass,2)+pow(boson.P(),2)) );
}

void loadInstrMETWeights(
    bool applyNvtxWeights, bool applyPtWeights, bool applyMassLineshape,
    std::map<TString, std::map<double, std::pair<double, double>>> &NVtxWeight_map,
    std::map<TString, std::map<double, std::pair<double, double>>> &PtWeight_map,
    std::map<TString, TH1*> &LineshapeMassWeight_map,
    std::vector<std::string> const &v_jetCat,
    Options const &options) {
  if (applyNvtxWeights){
    LOG_DEBUG << "NVtx weights will be applied.";
    std::string nvtxFile = FileInPath::Resolve(Options::NodeAs<std::string>(
      options.GetConfig(), {"photon_reweighting", "nvtx_reweighting"}));
    NVtxWeight_map["_ll"] = utils::TH1toMap(nvtxFile, "WeightHisto_ll");
  }
  if (applyPtWeights){
    LOG_DEBUG << "Pt weights will be applied.";
    if (!applyNvtxWeights) {
      LOG_WARN << "pT weights are supposed to be applied only on top of nvtx weights.";
    }
    std::string ptFile = FileInPath::Resolve(Options::NodeAs<std::string>(
      options.GetConfig(), {"photon_reweighting", "pt_reweighting"}));
    PtWeight_map["_ll"] = utils::TH1toMap(ptFile, "WeightHisto_ll");
    for (unsigned int i =0; i < v_jetCat.size(); i++){
      PtWeight_map["_ll"+v_jetCat[i]] = utils::TH1toMap(ptFile, "WeightHisto"+v_jetCat[i]+"_ll");
    }
  }
  if (applyMassLineshape){
    LOG_DEBUG << "Mass lineshape will be applied.";
    std::string massFile = FileInPath::Resolve(Options::NodeAs<std::string>(
      options.GetConfig(), {"photon_reweighting", "mass_lineshape"}));
    TFile *f_weight_lineshape = TFile::Open((TString) massFile);
    LineshapeMassWeight_map["_ll"] = (TH1D*) ((TH1D*) f_weight_lineshape->Get("WeightHisto_ll"))->Clone();
  }
}

void loadMeanWeights(
    bool applyMeanWeights,
    std::map<TString, std::map<double, double>> &meanWeight_map,
    std::vector<std::string> const &v_jetCat,
    Options const &options) {
  if (applyMeanWeights) {
    std::string meanWeightsFile = FileInPath::Resolve(Options::NodeAs<std::string>(
      options.GetConfig(), {"photon_reweighting", "mean_weights"}));
    TFile *file = TFile::Open((TString) meanWeightsFile);
    for (unsigned int i = 0 ; i < v_jetCat.size() ; i++) {
      TH1D *histo = (TH1D*) file->Get((TString)"mean_weights_tot"+v_jetCat[i]);
      meanWeight_map[v_jetCat[i]] = utils::TH1toMap(histo);
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

