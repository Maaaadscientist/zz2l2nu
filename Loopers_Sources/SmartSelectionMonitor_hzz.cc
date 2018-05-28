#include "../Loopers_Includes/SmartSelectionMonitor_hzz.h"

bool SmartSelectionMonitor_hzz::declareHistos(){ //FIXME: Later, will take an array as input for the binnings.
  addHistogram(new TH1F("totEventInBaobab",";Number of events in Baobab;Events",500,0,100));
  addHistogram(new TH1F("pile-up",";Number of PU events;Events",100,0,100));
  addHistogram(new TH1F("truth-pile-up",";Truth number of PU events;Events",100,0,100));
  addHistogram(new TH1F("reco-vtx",";Number of reco vtx;Events",100,0,100)); //Use for Photon reweighting method, don't change binning if you don't know what you're doing

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
  addHistogram(new TH1F("mT",";m_{T} (GeV);Events",nmTAxis-1,mTaxis));

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

  addHistogram(new TH1F("custom_HT",";#Sigma p_{T}^{gen jets};Events",200,0,2000));
  addHistogram(new TH1F("selJetsHT",";#Sigma p_{T}^{sel jets};Events",200,0,2000));

  //Control histo for closure test
  addHistogram(new TH1F("DeltaPhi_MET_Phot",";#Delta #phi(#gamma,E_{T}^{miss});Events", 50, 0, 4));
  addHistogram(new TH1F("DeltaPhi_MET_Jet",";min(#Delta#phi(jet,E_{T}^{miss}));Events", 50, 0, 4));
  addHistogram(new TH1F("METoverPt_zoom",";MET/p^{boson}_{T};Events",40,0, 2));
  addHistogram(new TH1F("eta_Z",";#eta_{Z/#gamma};Events",80, -4, 4));
  addHistogram(new TH1F("pT_jet0",";p_{T} of jet 0;Events",200,0,800));
  addHistogram(new TH1F("MET_HT300",";Missing transverse energy (GeV);Events",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt0-300",";Missing transverse energy (GeV);Events (Pt0-300)",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt300-400",";Missing transverse energy (GeV);Events (Pt300-400)",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt400-600",";Missing transverse energy (GeV);Events (Pt400-600)",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt600-Inf",";Missing transverse energy (GeV);Events (Pt600-Inf)",nMETAxis-1,METaxis));

  

  //TProfile for closure
  addHistogram(new TProfile("METvsBosonPt",";boson p_{T} (GeV);MET profile (GeV)", nzptAxis-1, zptaxis, 0, 500));
  addHistogram(new TProfile("METvsMT",";m_{T} (GeV);MET profile (GeV)", nmTAxis-1, mTaxis, 0, 500));
  addHistogram(new TProfile("METvsDPhiMETBos",";#Delta #phi(boson,E_{T}^{miss});MET profile (GeV)", 20, 0, 4, 0, 500));
  addHistogram(new TProfile("METvsDPhiMETJet",";min(#Delta#phi(jet,E_{T}^{miss}));MET profile (GeV)",20, 0, 4, 0, 500));
  addHistogram(new TProfile("METvsJetPt",";leading jet p_{T} (GeV);MET profile (GeV)", 40, 0, 800, 0, 500));
  addHistogram(new TProfile("METvsNJets",";# jets;MET profile (GeV)", 10, 0, 10, 0, 500));
  addHistogram(new TProfile("METvsBosonEta",";boson #eta;MET profile (GeV)", 40, -4, 4, 0, 500));
  addHistogram(new TProfile("METvsHT",";HT;MET profile (GeV)", 15, 0, 1500, 0, 500));
  addHistogram(new TProfile("HTvsBosonEta",";boson #eta;HT profile (GeV)", 40, -4, 4, 0, 500));

  return true;
}

bool SmartSelectionMonitor_hzz::declareHistos_InstrMET(){ 
  addHistogram(new TH1F("totEventInBaobab",";Number of events in Baobab;Events",500,0,100));
  addHistogram(new TH1F("pile-up",";Number of PU events;Events",100,0,100));
  addHistogram(new TH1F("truth-pile-up",";Truth number of PU events;Events",100,0,100));
  addHistogram(new TH1F("reco-vtx",";Number of reco vtx;Events",100,0,100));   //Use for Photon reweighting method, don't change binning if you don't know what you're doing
  TH1F *h =(TH1F*) addHistogram(new TH1F("eventflow",";;Events",17,0,17));
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
  h->GetXaxis()->SetBinLabel(12,"Dilepton NVtx Reweighting");
  h->GetXaxis()->SetBinLabel(13,"Dilepton Pt Reweighting");
  h->GetXaxis()->SetBinLabel(14,"b-veto");
  h->GetXaxis()->SetBinLabel(15,"#Delta #phi(jet,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(16,"E_{T}^{miss}>80");
  h->GetXaxis()->SetBinLabel(17,"E_{T}^{miss}>125");
  addHistogram(new TH1F("custom_HT",";#Sigma p_{T}^{gen jets};Events",200,0,2000));
  addHistogram(new TH1F("selJetsHT",";#Sigma p_{T}^{sel jets};Events",200,0,2000));
  addHistogram(new TH1F("PFJetsHT",";#Sigma p_{T}^{PF tight jets > 5 GeV};Events",200,0,2000));
  addHistogram(new TH1F("pT_jet0",";p_{T} of jet 0;Events",200,0,800));
  addHistogram(new TH1F("pT_jet1",";p_{T} of jet 1;Events",200,0,800));
  addHistogram(new TH1F("pT_jet2",";p_{T} of jet 2;Events",200,0,800));
  addHistogram(new TH1F("pT_jet3",";p_{T} of jet 3;Events",200,0,800));

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
  addHistogram(new TH1F("mT",";m_{T} (GeV);Events",nmTAxis-1,mTaxis));
  addHistogram(new TH1F("mT_unif",";m_{T} (GeV);Events",300,0,1500));

  TH1F *hc = (TH1F*) addHistogram(new TH1F("jetCategory",";Jet Category;Events",3,0,3));
  hc->GetXaxis()->SetBinLabel(1,"= 0 jets");
  hc->GetXaxis()->SetBinLabel(2,"#geq 1 jets");
  hc->GetXaxis()->SetBinLabel(3,"vbf");
  addHistogram(new TH1F("nJets",";N_{jets #geq 30 GeV};Events",20,0,20)); //Jets using the same criteria as the ones for jet bin category
  addHistogram(new TH1F("nb_selPhoton",";N_{sel photons};Events",8,0,8)); //Jets using the same criteria as the ones for jet bin category

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

  //Compa old new
  addHistogram(new TH1F("qt_rebin",      ";Transverse momentum [GeV];Events / GeV",nzptAxis-1,zptaxis));
  addHistogram(new TH1F("nvtx",";Vertices;Events",100, 0, 100)); //50,0,50) );
  addHistogram(new TH1F("rho",";#rho;Events",100,0,50) );
  addHistogram(new TH1F("MET_phi",";<MET #phi>;Events", 25, -4, 4));
  addHistogram(new TH1F("DeltaPhi_MET_Phot",";#Delta #phi(#gamma,E_{T}^{miss});Events", 50, 0, 4));
  addHistogram(new TH1F("DeltaPhi_MET_Jet",";min(#Delta#phi(jet,E_{T}^{miss}));Events", 50, 0, 4));

  //Photon Study
  //R9
  addHistogram( new TH1F ("R9", ";R9;Events", 150,0,1.5));
  //Sigma ieta ieta
  addHistogram( new TH1F ("SigmaIetaIeta", ";#sigma_{i#eta i#eta};Events", 50, 0, 0.05));
  addHistogram( new TH1F ("SigmaIetaIeta_zoom", ";#sigma_{i#eta i#eta};Events", 40, 0, 0.02));
  //Isolation
  Double_t PfIsoChHadBins[] = {1, 2, 10, 50};
  Int_t NPfIsoChHadBins = sizeof(PfIsoChHadBins)/sizeof(PfIsoChHadBins[0]) - 1;
  addHistogram( new TH1F ("PfIsoChHad_showQCD", ";PfIsoChHad;Events", NPfIsoChHadBins, PfIsoChHadBins));
  addHistogram( new TH1F ("PfIsoChHad", ";PfIsoChHad;Events", 50, 0, 10));
  addHistogram( new TH1F ("PfIsoNeutralHad", ";PfIsoNeutralHad;Events", 100, 0, 20));
  addHistogram( new TH1F ("PfIsoPhot", ";PfIsoPhot;Events", 100, 0, 20));
  addHistogram( new TH1F ("IsoEcal", ";IsoEcal;Events", 150, 0, 60));
  addHistogram( new TH1F ("IsoHcal", ";IsoHcal;Events", 100, 0, 40));
  addHistogram( new TH1F ("IsoTk", ";IsoTk;Events", 100, 0, 40));
  addHistogram( new TH1F ("PfIsoPuChHad", ";PfIsoPuChHad;Events", 150, 0, 75));
  addHistogram( new TH1F ("PfIsoEcalClus", ";PfIsoEcalClus;Events", 100, 0, 20));
  addHistogram( new TH1F ("PfIsoHcalClus", ";PfIsoHcalClus;Events", 100, 0, 20));
  //HoverE
  addHistogram( new TH1F ("HoE", ";HoE;Events", 40, 0, 0.20));
  addHistogram( new TH1F ("HoE_zoom", ";HoE;Events", 30, 0, 0.03));
  //Rho corrected photon ID Iso
  addHistogram( new TH1F ("RhoCorrPfIsoChHad", ";Rho Corrected PfIsoChHad;Events", 50, 0, 10));
  addHistogram( new TH1F ("RhoCorrPfIsoChHad_zoom", ";Rho Corrected PfIsoChHad;Events", 50, 0, 0.5));
  addHistogram( new TH1F ("RhoCorrPfIsoNeutralHad", ";Rho Corrected PfIsoNeutralHad;Events", 100, 0, 20));
  addHistogram( new TH1F ("RhoCorrPfIsoNeutralHad_zoom", ";Rho Corrected PfIsoNeutralHad;Events", 100, 0, 20));
  addHistogram( new TH1F ("RhoCorrPfIsoPhot", ";Rho Corrected PfIsoPhot;Events", 100, 0, 20));
  addHistogram( new TH1F ("RhoCorrPfIsoPhot_zoom", ";Rho Corrected PfIsoPhot;Events", 100, 0, 20));


  //More MET variables
  addHistogram(new TH1F("METsigx2",";Significance x^2 (C(0,0)) of the Missing transverse energy (GeV);Events",500,0, 1000));
  addHistogram(new TH1F("METsigxy",";Significance x-y (C(0,1)) of the Missing transverse energy (GeV);Events",500,0, 1000));
  addHistogram(new TH1F("METsigy2",";Significance y^2 (C(1,1)) of the Missing transverse energy (GeV);Events",500,0, 1000));
  addHistogram(new TH1F("METsig",";Significance of the Missing transverse energy (GeV);Events",100,0, 200));
  addHistogram(new TH1F("METoverPt",";MET/p^{boson}_{T};Events",100,0, 100));
  addHistogram(new TH1F("METoverPt_zoom",";MET/p^{boson}_{T};Events",40,0, 2));
  addHistogram(new TH1F("METperp","MET_{#perp}", 100, -200, 200)); 
  addHistogram(new TH1F("METpar","MET_{#parallel}+q_{T}", 100, -200, 200)); 

  addHistogram(new TH1F("MET_HT300",";Missing transverse energy (GeV);Events (HT>300GeV)",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt0-300",";Missing transverse energy (GeV);Events (Pt0-300)",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt300-400",";Missing transverse energy (GeV);Events (Pt300-400)",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt400-600",";Missing transverse energy (GeV);Events (Pt400-600)",nMETAxis-1,METaxis));
  addHistogram(new TH1F("MET_Pt600-Inf",";Missing transverse energy (GeV);Events (Pt600-Inf)",nMETAxis-1,METaxis));



  //TProfile for closure
  addHistogram(new TProfile("METvsBosonPt",";boson p_{T} (GeV);MET profile (GeV)", nzptAxis-1, zptaxis, 0, 500));
  addHistogram(new TProfile("METvsMT",";m_{T} (GeV);MET profile (GeV)", nmTAxis-1, mTaxis, 0, 500));
  addHistogram(new TProfile("METvsDPhiMETBos",";#Delta #phi(boson,E_{T}^{miss});MET profile (GeV)", 20, 0, 4, 0, 500));
  addHistogram(new TProfile("METvsDPhiMETJet",";min(#Delta#phi(jet,E_{T}^{miss}));MET profile (GeV)",20, 0, 4, 0, 500));
  addHistogram(new TProfile("METvsJetPt",";leading jet p_{T} (GeV);MET profile (GeV)", 40, 0, 800, 0, 500));
  addHistogram(new TProfile("METvsNJets",";# jets;MET profile (GeV)", 10, 0, 10, 0, 500));
  addHistogram(new TProfile("METvsBosonEta",";boson #eta;MET profile (GeV)", 40, -4, 4, 0, 500));
  addHistogram(new TProfile("METvsHT",";HT;MET profile (GeV)", 15, 0, 1500, 0, 500));
  addHistogram(new TProfile("HTvsBosonEta",";boson #eta;HT profile (GeV)", 40, -4, 4, 0, 500));

  return true;  
} 

bool SmartSelectionMonitor_hzz::fillHistoForAllCategories(TString name, double variable, evt currentEvt, TString tag, double weight){
  fillHisto(name, tag+currentEvt.s_lepCat+currentEvt.s_jetCat, variable, weight,true); //The 'true' means that we normalize by the bin width!
  fillHisto(name, tag+currentEvt.s_lepCat, variable, weight, true);
  fillHisto(name, tag, variable, weight, true);
  return true;
}

bool SmartSelectionMonitor_hzz::fillHistoForAllCategories(TString name, double variable, photon_evt currentEvt, TString tag, double weight){
  fillHisto(name, tag+currentEvt.s_jetCat, variable, weight,true); //The 'true' means that we normalize by the bin width!
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
  if (tag.Contains("inZpeak")){
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
bool SmartSelectionMonitor_hzz::fillAnalysisHistos_InstrMET(photon_evt currentEvt, TString tag, double weight){
  std::map<std::string, double> data;
  data["mT"] = currentEvt.MT;
  data["mT_unif"] = currentEvt.MT;
  data["M_Z"] = currentEvt.M;
  data["pT_Z"] = currentEvt.pT;
  data["pT_Z_unif"] = currentEvt.pT;
  data["eta_Z"] = currentEvt.eta;
  data["MET"] = currentEvt.MET;
  data["MET_unif"] = currentEvt.MET;
  data["nJets"] = currentEvt.nJets;
  for(std::map<std::string,double>::iterator it = data.begin() ; it != data.end() ; it++) fillHistoForAllCategories(it->first, it->second, currentEvt, tag, weight);
  return true;
}

bool SmartSelectionMonitor_hzz::fillPhotonIDHistos_InstrMET(photon_evt currentEvt, TString tag, double weight){
  std::map<std::string, double> data;
  data["mT"] = currentEvt.MT;
  data["M_Z"] = currentEvt.M;
  data["pT_Z"] = currentEvt.pT;
  data["eta_Z"] = currentEvt.eta;
  data["MET"] = currentEvt.MET;
  data["nJets"] = currentEvt.nJets;
  data["nb_selPhoton"] = currentEvt.nb_selPhoton;
  data["MET_phi"] = currentEvt.METphi;
  data["DeltaPhi_MET_Jet"] = currentEvt.deltaPhi_MET_jet;
  data["DeltaPhi_MET_Phot"] = currentEvt.deltaPhi_MET_phot;
  data["METoverPt"] = currentEvt.METoPT;
  data["METoverPt_zoom"] = currentEvt.METoPT;
  data["METpar"] = currentEvt.METpar;
  data["METperp"] = currentEvt.METperp;
  data["METsig"] = currentEvt.METsig;
  data["nvtx"] = currentEvt.nVtx;
  data["HoE"] = currentEvt.HoE;
  data["HoE_zoom"] = currentEvt.HoE;
  data["SigmaIetaIeta"] = currentEvt.sigmaIEtaIEta;
  data["SigmaIetaIeta_zoom"] = currentEvt.sigmaIEtaIEta;
  data["RhoCorrPfIsoChHad"] = currentEvt.chIsoRhoCorr;
  data["RhoCorrPfIsoChHad_zoom"] = currentEvt.chIsoRhoCorr;
  data["RhoCorrPfIsoNeutralHad"] = currentEvt.neuIsoRhoCorr;
  data["RhoCorrPfIsoNeutralHad_zoom"] = currentEvt.neuIsoRhoCorr;
  data["RhoCorrPfIsoPhot"] = currentEvt.phoIsoRhoCorr;
  data["RhoCorrPfIsoPhot_zoom"] = currentEvt.phoIsoRhoCorr;
  data["R9"] = currentEvt.R9;
  data["rho"] = currentEvt.rho;
  data["pT_jet0"] = currentEvt.jet0_pT;
  data["pT_jet1"] = currentEvt.jet1_pT;
  data["pT_jet2"] = currentEvt.jet2_pT;
  data["pT_jet3"] = currentEvt.jet3_pT;
  data["selJetsHT"] = currentEvt.HT_selJets;
  data["PFJetsHT"] = currentEvt.HT_PFTightJetsAbove5GeV;
  for(std::map<std::string,double>::iterator it = data.begin() ; it != data.end() ; it++) fillHistoForAllCategories(it->first, it->second, currentEvt, tag, weight);
  return true;
}
