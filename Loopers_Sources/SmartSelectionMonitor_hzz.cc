#include "../Loopers_Includes/SmartSelectionMonitor_hzz.h"

bool SmartSelectionMonitor_hzz::declareHistos(){ //FIXME: Later, will take an array as input for the binnings.
  addHistogram(new TH1F("totEventInBaobab",";Number of events in Baobab;Events",500,0,100));
  addHistogram(new TH1F("pile-up",";Number of PU events;Events",100,0,100));
  addHistogram(new TH1F("truth-pile-up",";Truth number of PU events;Events",100,0,100));
  addHistogram(new TH1F("reco-vtx",";Number of reco vtx;Events",100,0,100)); //Use for Photon reweighting method, don't change binning if you don't know what you're doing

  Double_t VtxBins_ee_eq0jets[] = {0,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,38,41,70,100};   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  Double_t VtxBins_ee_geq1jets[] = {0,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35, 37, 39, 41, 43, 45, 47, 49, 51, 55, 60, 100};   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  Double_t VtxBins_ee_vbf[] = {0,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,34, 37, 40, 45, 55, 100};   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  Double_t VtxBins_mumu_eq0jets[] = {0,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35, 37, 39, 42, 45, 50, 55, 70, 100};   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  Double_t VtxBins_mumu_geq1jets[] = {0,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35, 37, 39, 41, 43, 45, 47, 49, 51, 55, 60, 100};   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  Double_t VtxBins_mumu_vbf[] = {0,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,34, 37, 40, 45, 55, 100};   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  Int_t NVtxBins_ee_eq0jets = sizeof(VtxBins_ee_eq0jets)/sizeof(VtxBins_ee_eq0jets[0]) - 1;
  Int_t NVtxBins_ee_geq1jets = sizeof(VtxBins_ee_geq1jets)/sizeof(VtxBins_ee_geq1jets[0]) - 1;
  Int_t NVtxBins_ee_vbf = sizeof(VtxBins_ee_vbf)/sizeof(VtxBins_ee_vbf[0]) - 1;
  Int_t NVtxBins_mumu_eq0jets = sizeof(VtxBins_mumu_eq0jets)/sizeof(VtxBins_mumu_eq0jets[0]) - 1;
  Int_t NVtxBins_mumu_geq1jets = sizeof(VtxBins_mumu_geq1jets)/sizeof(VtxBins_mumu_geq1jets[0]) - 1;
  Int_t NVtxBins_mumu_vbf = sizeof(VtxBins_mumu_vbf)/sizeof(VtxBins_mumu_vbf[0]) - 1;

  addHistogram(new TH1F("reco-vtx_ee_eq0jets",";Number of reco vtx;Events",NVtxBins_ee_eq0jets, VtxBins_ee_eq0jets));
  addHistogram(new TH1F("reco-vtx_ee_geq1jets",";Number of reco vtx;Events",NVtxBins_ee_geq1jets, VtxBins_ee_geq1jets));
  addHistogram(new TH1F("reco-vtx_ee_vbf",";Number of reco vtx;Events",NVtxBins_ee_vbf, VtxBins_ee_vbf));
  addHistogram(new TH1F("reco-vtx_mumu_eq0jets",";Number of reco vtx;Events",NVtxBins_mumu_eq0jets, VtxBins_mumu_eq0jets));
  addHistogram(new TH1F("reco-vtx_mumu_geq1jets",";Number of reco vtx;Events",NVtxBins_mumu_geq1jets, VtxBins_mumu_geq1jets));
  addHistogram(new TH1F("reco-vtx_mumu_vbf",";Number of reco vtx;Events",NVtxBins_mumu_vbf, VtxBins_mumu_vbf));

  TH1F *h =(TH1F*) addHistogram(new TH1F("eventflow",";;Events",10,0,10));
  h->GetXaxis()->SetBinLabel(1,"skimmed");
  h->GetXaxis()->SetBinLabel(2,"#geq 2 iso leptons");
  h->GetXaxis()->SetBinLabel(3,"|M-91|<15");
  h->GetXaxis()->SetBinLabel(4,"p_{T}>55");
  h->GetXaxis()->SetBinLabel(5,"#Delta #phi(Z,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(6,"3^{rd}-lepton veto");
  h->GetXaxis()->SetBinLabel(7,"b-veto");
  h->GetXaxis()->SetBinLabel(8,"#Delta #phi(jet,E_{T}^{miss})>0.5");
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
  addHistogram(new TH1F("pT_Z",";p_{T,Z};Events",nzptAxis-1,zptaxis));  //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  addHistogram(new TH1F("MET",";Missing transverse energy (GeV);Events",nMETAxis-1,METaxis));
  addHistogram(new TH1F("METphi",";#phi of missing transverse energy;Events",80,-4.,4.));
  addHistogram(new TH1F("mT",";m_{T};Events",nmTAxis-1,mTaxis));

  TH1F *hc = (TH1F*) addHistogram(new TH1F("jetCategory",";Jet Category;Events",3,0,3));
  hc->GetXaxis()->SetBinLabel(1,"= 0 jets");
  hc->GetXaxis()->SetBinLabel(2,"#geq 1 jets");
  hc->GetXaxis()->SetBinLabel(3,"vbf");
  addHistogram(new TH1F("nJets",";N_{jets #geq 30 GeV};Events",20,0,20)); //Jets using the same criteria as the ones for jet bin category

  addHistogram(new TH1F("pT_l1",";p_{T} of lepton 1;Events",16,10,250));
  addHistogram(new TH1F("pT_l2",";p_{T} of lepton 2;Events",16,10,250));
  addHistogram(new TH1F("eta_l1",";#eta of lepton 1;Events",20,-2.5,2.5));
  addHistogram(new TH1F("eta_l2",";#eta of lepton 2;Events",20,-2.5,2.5));

  addHistogram(new TH1F("runNumber",";run number",100,273158, 284044));
  return true;
}

bool SmartSelectionMonitor_hzz::declareHistos_InstrMET(){ 
  addHistogram(new TH1F("totEventInBaobab",";Number of events in Baobab;Events",500,0,100));
  addHistogram(new TH1F("pile-up",";Number of PU events;Events",100,0,100));
  addHistogram(new TH1F("truth-pile-up",";Truth number of PU events;Events",100,0,100));
  addHistogram(new TH1F("reco-vtx",";Number of reco vtx;Events",100,0,100));   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  TH1F *h =(TH1F*) addHistogram(new TH1F("eventflow",";;Events",15,0,15));
  h->GetXaxis()->SetBinLabel(1,"skimmed (ID, p_{T}, trigger, ...)"); //what is coming from the bonzai
  h->GetXaxis()->SetBinLabel(2,"prescale weight");
  h->GetXaxis()->SetBinLabel(3,"photon efficiency SF");
  h->GetXaxis()->SetBinLabel(4,"PU reweighting");
  h->GetXaxis()->SetBinLabel(5,"MET filters");
  h->GetXaxis()->SetBinLabel(6,"QCD/GJets double counting");
  h->GetXaxis()->SetBinLabel(7,"ZnunuG k-factor");
  h->GetXaxis()->SetBinLabel(8,"WJets double counting");
  h->GetXaxis()->SetBinLabel(9,"p_{T, #gamma}>55");
  h->GetXaxis()->SetBinLabel(10,"#Delta #phi(#gamma,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(11,"no lepton veto");
  h->GetXaxis()->SetBinLabel(12,"b-veto");
  h->GetXaxis()->SetBinLabel(13,"#Delta #phi(jet,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(14,"E_{T}^{miss}>80");
  h->GetXaxis()->SetBinLabel(15,"E_{T}^{miss}>125");
  addHistogram(new TH1F("pT_mu",";p_{T} of muon;Events",200,0,800));
  addHistogram(new TH1F("pT_e",";p_{T} of electron;Events",200,0,800));
  addHistogram(new TH1F("nb_mu",";number of muons;Events",10,0,10));
  addHistogram(new TH1F("nb_e",";number of electrons;Events",10,0,10));
  addHistogram(new TH1F("nb_pho",";number of photons;Events",10,0,10));

  Double_t METaxis[]={0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,125,150,175,200,250,300,400,500,600,700,800,900,1000};
  Double_t zptaxis[]= {0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360,375,390,405,435,465,495,525,555,585,615,675,735,795,855,975,1500};  //Use for Photon reweighting method, don't change binning if you don't know what you're doing
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

  //FIXME Compa old new
  addHistogram( new TH1F( "qt",         ";Transverse momentum [GeV];Events / GeV",3000,0,3000));
  addHistogram( new TH1F( "qtraw",      ";Transverse momentum [GeV];Events / GeV",3000,0,3000));
  addHistogram( new TH1F( "qtMet125",      ";Transverse momentum [GeV];Events / GeV",3000,0,3000));
  addHistogram( new TH1F( "qt_rebin",         ";Transverse momentum [GeV];Events / GeV",nzptAxis-1,zptaxis));
  addHistogram( new TH1F( "qtraw_rebin",      ";Transverse momentum [GeV];Events / GeV",nzptAxis-1,zptaxis));
  addHistogram( new TH1F( "qtMet125_rebin",      ";Transverse momentum [GeV];Events / GeV",nzptAxis-1,zptaxis));
  addHistogram( new TH1F( "nvtx",";Vertices;Events",100, 0, 100)); //50,0,50) );
  addHistogram( new TH1F( "nvtxraw",";Vertices;Events",100, 0, 100 )); //50,0,50) );
  addHistogram( new TH1F( "rho",";#rho;Events",100,0,50) );
  addHistogram(new TH1F("MET_phi",";<MET #phi>;Events", 25, -4, 4));
  addHistogram(new TH1F("DeltaPhi_MET_Phot",";#Delta #phi(#gamma,E_{T}^{miss});Events", 50, 0, 4));
  addHistogram(new TH1F("DeltaPhi_MET_Jet",";#Delta min(#phi(jet,E_{T}^{miss}));Events", 50, 0, 4));
  Double_t ZPtBins[] = {0,30,36,50,75,90,120,165,3000}; //which correspond to our trigger threshold. Starting at 165 the prescale is set to 1 so we stop binning there.
  Int_t NZPtBins = sizeof(ZPtBins)/sizeof(ZPtBins[0]) - 1;

  Double_t VtxBins[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100};
  Int_t NVtxBins = sizeof(VtxBins)/sizeof(VtxBins[0]) - 1;
  
  Double_t VtxBins_ee_eq0jets[] = {0,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,38,41,70,100};   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  Double_t VtxBins_ee_geq1jets[] = {0,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35, 37, 39, 41, 43, 45, 47, 49, 51, 55, 60, 100};   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  Double_t VtxBins_ee_vbf[] = {0,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,34, 37, 40, 45, 55, 100};   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  Double_t VtxBins_mumu_eq0jets[] = {0,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35, 37, 39, 42, 45, 50, 55, 70, 100};   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  Double_t VtxBins_mumu_geq1jets[] = {0,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35, 37, 39, 41, 43, 45, 47, 49, 51, 55, 60, 100};   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  Double_t VtxBins_mumu_vbf[] = {0,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,34, 37, 40, 45, 55, 100};   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  Int_t NVtxBins_ee_eq0jets = sizeof(VtxBins_ee_eq0jets)/sizeof(VtxBins_ee_eq0jets[0]) - 1;
  Int_t NVtxBins_ee_geq1jets = sizeof(VtxBins_ee_geq1jets)/sizeof(VtxBins_ee_geq1jets[0]) - 1;
  Int_t NVtxBins_ee_vbf = sizeof(VtxBins_ee_vbf)/sizeof(VtxBins_ee_vbf[0]) - 1;
  Int_t NVtxBins_mumu_eq0jets = sizeof(VtxBins_mumu_eq0jets)/sizeof(VtxBins_mumu_eq0jets[0]) - 1;
  Int_t NVtxBins_mumu_geq1jets = sizeof(VtxBins_mumu_geq1jets)/sizeof(VtxBins_mumu_geq1jets[0]) - 1;
  Int_t NVtxBins_mumu_vbf = sizeof(VtxBins_mumu_vbf)/sizeof(VtxBins_mumu_vbf[0]) - 1;

  addHistogram( new TH2F ("zpt_vs_nvtx_ee_eq0jets",";zpt;#vertices",NZPtBins, ZPtBins, NVtxBins_ee_eq0jets, VtxBins_ee_eq0jets) );   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  addHistogram( new TH2F ("zpt_vs_nvtx_ee_geq1jets",";zpt;#vertices",NZPtBins, ZPtBins, NVtxBins_ee_geq1jets, VtxBins_ee_geq1jets) );   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  addHistogram( new TH2F ("zpt_vs_nvtx_ee_vbf",";zpt;#vertices",NZPtBins, ZPtBins, NVtxBins_ee_vbf, VtxBins_ee_vbf) );   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  addHistogram( new TH2F ("zpt_vs_nvtx_mumu_eq0jets",";zpt;#vertices",NZPtBins, ZPtBins, NVtxBins_mumu_eq0jets, VtxBins_mumu_eq0jets) );   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  addHistogram( new TH2F ("zpt_vs_nvtx_mumu_geq1jets",";zpt;#vertices",NZPtBins, ZPtBins, NVtxBins_mumu_geq1jets, VtxBins_mumu_geq1jets) );   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  addHistogram( new TH2F ("zpt_vs_nvtx_mumu_vbf",";zpt;#vertices",NZPtBins, ZPtBins, NVtxBins_mumu_vbf, VtxBins_mumu_vbf) );   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  addHistogram( new TH2F ("zpt_vs_rho",";zpt;rho", NZPtBins, ZPtBins, NVtxBins, VtxBins) );
  
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
  data["METphi"] = currentEvt.METphi;
  data["nJets"] = currentEvt.nJets;
  if (tag=="inZpeak"){
    data["pT_l1"] = currentEvt.lep1pT;
    data["eta_l1"] = currentEvt.lep1eta;
    data["pT_l2"] = currentEvt.lep2pT;
    data["eta_l2"] = currentEvt.lep2eta;
    data["runNumber"] = currentEvt.runNumber;
    data["reco-vtx"] = currentEvt.nVtx;
  }
  for(std::map<std::string,double>::iterator it = data.begin() ; it != data.end() ; it++) fillHistoForAllCategories(it->first, it->second, currentEvt, tag, weight);
  return true;
}
bool SmartSelectionMonitor_hzz::fillAnalysisHistos_InstrMET(evt currentEvt, TString tag, double weight){
  std::map<std::string, double> data;
  data["mT"] = currentEvt.transverseMass;
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
