#include "../Includes/SmartSelectionMonitor_hzz.h"

bool SmartSelectionMonitor_hzz::declareHistos(){ //FIXME: Later, will take an array as input for the binnings.
  addHistogram(new TH1F("totEventInBaobab",";Number of events in Baobab;Events",500,0,100));
  addHistogram(new TH1F("pile-up",";Number of PU events;Events",50,0,50));
  addHistogram(new TH1F("truth-pile-up",";Truth number of PU events;Events",50,0,50));
  addHistogram(new TH1F("reco-vtx",";Number of reco vtx;Events",50,0,50));
  TH1F *h =(TH1F*) addHistogram(new TH1F("eventflow",";;Events",10,0,10));
  h->GetXaxis()->SetBinLabel(1,"skimmed");
  h->GetXaxis()->SetBinLabel(2,"#geq 2 iso leptons");
  h->GetXaxis()->SetBinLabel(3,"|M-91|<15");
  h->GetXaxis()->SetBinLabel(4,"p_{T}>55");
  h->GetXaxis()->SetBinLabel(5,"3^{rd}-lepton veto");
  h->GetXaxis()->SetBinLabel(6,"b-veto");
  h->GetXaxis()->SetBinLabel(7,"#Delta #phi(jet,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(8,"#Delta #phi(Z,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(9,"E_{T}^{miss}>80");
  h->GetXaxis()->SetBinLabel(10,"E_{T}^{miss}>125");
  addHistogram(new TH1F("pT_mu",";p_{T} of muon;Events",200,0,800));
  addHistogram(new TH1F("pT_e",";p_{T} of electron;Events",200,0,800));
  addHistogram(new TH1F("nb_mu",";number of muons;Events",10,0,10));
  addHistogram(new TH1F("nb_e",";number of electrons;Events",10,0,10));

  Double_t METaxis[]={0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,125,150,175,200,250,300,400,500,600,700,800,900,1000};
  Double_t zptaxis[]= {0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360,375,390,405,435,465,495,525,555,585,615,675,735,795,855,975,1500};
  Double_t mTaxis[]={100,120,140,160,180,200,220,240,260,280,300,325,350,375,400,450,500,600,700,800,900,1000,1500,2000};
  Int_t nMETAxis=sizeof(METaxis)/sizeof(Double_t);
  Int_t nzptAxis=sizeof(zptaxis)/sizeof(Double_t);
  Int_t nmTAxis=sizeof(mTaxis)/sizeof(Double_t);
  addHistogram(new TH1F("M_Z",";M_{Z};Events",100,76,106));
  addHistogram(new TH1F("pT_Z",";p_{T,Z};Events",nzptAxis-1,zptaxis));
  addHistogram(new TH1F("MET",";Missing transverse energy (GeV);Events",nMETAxis-1,METaxis));
  addHistogram(new TH1F("mT",";m_{T};Events",nmTAxis-1,mTaxis));

  TH1F *hc = (TH1F*) addHistogram(new TH1F("jetCategory",";Jet Category;Events",3,0,3));
  hc->GetXaxis()->SetBinLabel(1,"= 0 jets");
  hc->GetXaxis()->SetBinLabel(2,"#geq 1 jets");
  hc->GetXaxis()->SetBinLabel(3,"vbf");
  addHistogram(new TH1F("nJets",";N_{jets #geq 30 GeV};Events",20,0,20)); //Jets using the same criteria as the ones for jet bin category
  return true;
}

bool SmartSelectionMonitor_hzz::declareHistos_InstrMET(){ //FIXME: Like stated above, later this will take an array as input for the binnings.
  addHistogram(new TH1F("totEventInBaobab",";Number of events in Baobab;Events",500,0,100));
  addHistogram(new TH1F("pile-up",";Number of PU events;Events",100,0,100));
  addHistogram(new TH1F("truth-pile-up",";Truth number of PU events;Events",100,0,100));
  addHistogram(new TH1F("reco-vtx",";Number of reco vtx;Events",100,0,100));
  TH1F *h =(TH1F*) addHistogram(new TH1F("eventflow",";;Events",12,0,12));
  h->GetXaxis()->SetBinLabel(1,"skimmed (ID, p_{T}, trigger, ...)"); //what is coming from the bonzai
  h->GetXaxis()->SetBinLabel(2,"prescale weight");
  h->GetXaxis()->SetBinLabel(3,"photon efficiency SF");
  h->GetXaxis()->SetBinLabel(4,"PU reweighting");
  h->GetXaxis()->SetBinLabel(5,"MET filters");
  h->GetXaxis()->SetBinLabel(6,"p_{T, #gamma}>55");
  h->GetXaxis()->SetBinLabel(7,"no lepton veto");
  h->GetXaxis()->SetBinLabel(8,"b-veto");
  h->GetXaxis()->SetBinLabel(9,"#Delta #phi(jet,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(10,"#Delta #phi(#gamma,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(11,"E_{T}^{miss}>80");
  h->GetXaxis()->SetBinLabel(12,"E_{T}^{miss}>125");
  addHistogram(new TH1F("pT_mu",";p_{T} of muon;Events",200,0,800));
  addHistogram(new TH1F("pT_e",";p_{T} of electron;Events",200,0,800));
  addHistogram(new TH1F("nb_mu",";number of muons;Events",10,0,10));
  addHistogram(new TH1F("nb_e",";number of electrons;Events",10,0,10));
  addHistogram(new TH1F("nb_pho",";number of photons;Events",10,0,10));

  Double_t METaxis[]={0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,125,150,175,200,250,300,400,500,600,700,800,900,1000};
  Double_t zptaxis[]= {0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360,375,390,405,435,465,495,525,555,585,615,675,735,795,855,975,1500};
  Double_t mTaxis[]={100,120,140,160,180,200,220,240,260,280,300,325,350,375,400,450,500,600,700,800,900,1000,1500,2000};
  Int_t nMETAxis=sizeof(METaxis)/sizeof(Double_t);
  Int_t nzptAxis=sizeof(zptaxis)/sizeof(Double_t);
  Int_t nmTAxis=sizeof(mTaxis)/sizeof(Double_t);
  addHistogram(new TH1F("M_Z",";M_{Z/#gamma};Events",100,76,106));
  addHistogram(new TH1F("pT_Z",";p_{T,Z/#gamma};Events",nzptAxis-1,zptaxis));
  addHistogram(new TH1F("pT_Z_unif",";p_{T,Z/#gamma};Events",300, 0, 1500));
  addHistogram(new TH1F("eta_Z",";#eta_{Z/#gamma};Events",80, -4, 4));
  addHistogram(new TH1F("MET",";Missing transverse energy (GeV);Events",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_unif",";Missing transverse energy (GeV);Events",200,0,1000));
  addHistogram(new TH1F("mT",";m_{T};Events",nmTAxis-1,mTaxis));
  addHistogram(new TH1F("mT_unif",";m_{T};Events",300,0,1500));

  TH1F *hc = (TH1F*) addHistogram(new TH1F("jetCategory",";Jet Category;Events",3,0,3));
  hc->GetXaxis()->SetBinLabel(1,"= 0 jets");
  hc->GetXaxis()->SetBinLabel(2,"#geq 1 jets");
  hc->GetXaxis()->SetBinLabel(3,"vbf");
  addHistogram(new TH1F("nJets",";N_{jets #geq 30 GeV};Events",20,0,20)); //Jets using the same criteria as the ones for jet bin category

  TH1F *h_metFilter =(TH1F*) addHistogram(new TH1F("metFilters",";;Events passing MET filters",26+1,0,26+1)); //We add +1 everywhere so we can start with bin 0
  h_metFilter->GetXaxis()->SetBinLabel(26+1,"all"); //1 is already taken so we take the last bin of the MET filters +1, i.e 25+1=26
  h_metFilter->GetXaxis()->SetBinLabel(0+1,"duplicateMuons");
  h_metFilter->GetXaxis()->SetBinLabel(1+1,"badMuons");
  h_metFilter->GetXaxis()->SetBinLabel(2+1,"noBadMuons");
  h_metFilter->GetXaxis()->SetBinLabel(14+1,"primary vertex filter");
  h_metFilter->GetXaxis()->SetBinLabel(8+1,"beam halo filter");
  h_metFilter->GetXaxis()->SetBinLabel(3+1,"HBHE noise filter");
  h_metFilter->GetXaxis()->SetBinLabel(4+1,"HBHEiso noise filter");
  h_metFilter->GetXaxis()->SetBinLabel(12+1,"ECAL TP filter");
  h_metFilter->GetXaxis()->SetBinLabel(15+1,"ee badSC noise filter");
  h_metFilter->GetXaxis()->SetBinLabel(24+1,"badPFMuon");
  h_metFilter->GetXaxis()->SetBinLabel(25+1,"badCharged hadron");

  return true;  
} 

bool SmartSelectionMonitor_hzz::fillHistoForAllCategories(TString name, double variable, evt currentEvt, TString tag, double weight){
  fillHisto(name, tag+currentEvt.s_lepCat+currentEvt.s_jetCat, variable, weight,true); //The 'true' means that we normalize by the bin width!
  fillHisto(name, tag+currentEvt.s_lepCat, variable, weight, true);
  fillHisto(name, tag+currentEvt.s_jetCat, variable, weight, true);
  fillHisto(name, tag, variable, weight, true);
  return true;
}

bool SmartSelectionMonitor_hzz::fillAnalysisHistos(evt currentEvt, TString tag, double weight){
  std::map<std::string, double> data;
  data["mT"] = currentEvt.transverseMass;
  data["M_Z"] = currentEvt.MZ;
  data["pT_Z"] = currentEvt.pTZ;
  data["MET"] = currentEvt.MET;
  data["nJets"] = currentEvt.nJets;
  for(std::map<std::string,double>::iterator it = data.begin() ; it != data.end() ; it++) fillHistoForAllCategories(it->first, it->second, currentEvt, tag, weight);
  return true;
}
bool SmartSelectionMonitor_hzz::fillAnalysisHistos_InstrMET(evt currentEvt, TString tag, double weight){
  std::map<std::string, double> data;
  data["mT"] = currentEvt.transverseMass; //To keep the shape, it should be normalized to the bin width. This is done in the data-MC comparison script
  data["mT_unif"] = currentEvt.transverseMass;
  data["M_Z"] = currentEvt.MZ;
  data["pT_Z"] = currentEvt.pTZ;
  data["pT_Z_unif"] = currentEvt.pTZ;
  data["eta_Z"] = currentEvt.etaZ;
  data["MET"] = currentEvt.MET;
  data["MET_unif"] = currentEvt.MET;
  data["nJets"] = currentEvt.nJets;
  for(std::map<std::string,double>::iterator it = data.begin() ; it != data.end() ; it++) fillHistoForAllCategories(it->first, it->second, currentEvt, tag, weight);
  return true;
}
