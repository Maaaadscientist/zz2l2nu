#define HZZ2l2nuLooper_cxx
#include "../Includes/HZZ2l2nuLooper.h"
#include "../Includes/SmartSelectionMonitor.h"
#include "../Includes/Utils.h"
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>

/*
double deltaR(TLorentzVector v1, TLorentzVector v2)
{
   double dEta = v1.Eta()-v2.Eta();
   double dPhi = v1.Phi()-v2.Phi();
   return sqrt(pow(dEta,2) + pow(dPhi,2));
}
*/

void HZZ2l2nuLooper::Loop()
{
   if (fChain == 0) return;

   SmartSelectionMonitor mon;
   mon.addHistogram(new TH1F("pile-up",";Number of PU events;Events",50,0,50));
   TH1F *h =(TH1F*) mon.addHistogram(new TH1F("eventflow",";;Events",10,0,10));
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
   mon.addHistogram(new TH1F("pT_mu",";p_{T} of muon;Events",200,0,800));
   mon.addHistogram(new TH1F("pT_e",";p_{T} of electron;Events",200,0,800));
   mon.addHistogram(new TH1F("nb_mu",";number of muons;Events",10,0,10));
   mon.addHistogram(new TH1F("nb_e",";number of electrons;Events",10,0,10));
   mon.addHistogram(new TH1F("M_Z",";M_{Z};Events",200,0,200));
   mon.addHistogram(new TH1F("MET",";Missing transverse energy (GeV);Events",200,0,800));
   TH1F *hc = (TH1F*) mon.addHistogram(new TH1F("jetCategory",";Jet Category;Events",3,0,3));
   hc->GetXaxis()->SetBinLabel(1,"= 0 jets");
   hc->GetXaxis()->SetBinLabel(2,"#geq 1 jets");
   hc->GetXaxis()->SetBinLabel(3,"vbf");
   mon.addHistogram(new TH1F("mT",";m_{T};Events",200,0,800));
   mon.addHistogram(new TH1F("nJets",";N_{jets #geq 30 GeV};Events",20,0,20)); //Jets using the same criteria as the ones for jet bin category


   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   cout << "nb of entries in the input file =" << nentries << endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(jentry % 10000 ==0) cout << jentry << " of " << nentries << endl;

      double weight = 1.;

      mon.fillHisto("eventflow","all",0,weight);

      for(int i =0 ; i < MuPt->size() ; i++) mon.fillHisto("pT_mu","all",MuPt->at(i),weight);
      for(int i =0 ; i < ElPt->size() ; i++) mon.fillHisto("pT_e","all",ElPt->at(i),weight);
      mon.fillHisto("nb_mu","all",MuPt->size(),weight);
      mon.fillHisto("nb_e","all",ElPt->size(),weight);
      mon.fillHisto("pile-up","all",EvtPuCnt,weight);

      //Leptons passing final cuts, or additional leptons veto
      vector<TLorentzVector> selElectrons;
      vector<TLorentzVector> selMuons;
      vector<TLorentzVector> extraElectrons;
      vector<TLorentzVector> extraMuons;
      vector<TLorentzVector> selPhotons;
      vector<TLorentzVector> tagJets;
      vector<TLorentzVector> selJets;
      vector<double> btags;

      for(int i = 0 ; i<ElPt->size() ; i++){ //Electrons
         bool passEta = false, passIso = false, passId = false, passPt = false, passLoosePt = false, passLooseId = false;
	 TLorentzVector currentLepton; currentLepton.SetPtEtaPhiE(ElPt->at(i),ElEta->at(i),ElPhi->at(i),ElE->at(i));
	 //Id //Very temporary!!! Used without much cross-checking.
	 passId = ElId->at(i) & (1<<17);
	 passLooseId = ElId->at(i) & (1<<16);
	 int eta = fabs(ElEtaSc->at(i));//I took the supercluster eta since it's really the geometry which is taken here.
	 passEta = (eta<=2.5 && (eta>=1.5660 || eta<=1.4442));
	 //Iso //We use ElPfIsoRho for now, we'll see after if it's mandatory to refine it. Iso is applied only for the "tight" selection, not for the extra lepton veto.
	 if(eta>=1.5660 && ElPfIsoRho->at(i)<0.0646) passIso = true;
	 if(eta<=1.4442 && ElPfIsoRho->at(i)<0.0354) passIso = true; //Numbers are taken from llvv_fwk and have not been checked.
	 passPt = (currentLepton.Pt() >=25);
	 passLoosePt = (currentLepton.Pt() >=10);
	 if(passEta && passLooseId && passLoosePt && selElectrons.size()==2) extraElectrons.push_back(currentLepton); //No iso criteria for extra leptons.
         if(passEta && passIso && passId && passPt && selElectrons.size()<2) selElectrons.push_back(currentLepton);
      }
      for(int i = 0 ; i<MuPt->size() ; i++){ //Muons
         bool passEta = false, passIso = false, passId = false, passPt = false, passLoosePt = false, passLooseId = false, passSoftId = false, passSoftPt = false;
	 TLorentzVector currentLepton; currentLepton.SetPtEtaPhiE(MuPt->at(i),MuEta->at(i),MuPhi->at(i),MuE->at(i));
	 //Id //Very temporary!!! Used without much cross-checking.
	 passId = MuIdTight->at(i) & (1<<0); //Look at the first vertex, hence the bit 0.
	 passLooseId = MuId->at(i) & (1<<0);
	 passSoftId = false; //Don't know yet how to implement this thing. FIXME
	 int eta = fabs(MuEta->at(i));
	 passEta = (eta<=2.4);
	 //Iso //We use MuPfIso for now, we'll see after if it's mandatory to refine it. Iso is applied only for the "tight" selection, not for the extra lepton veto.
	 passIso = (MuPfIso->at(i)<0.15); //Numbers are taken from llvv_fwk and have not been checked.
	 passPt = (currentLepton.Pt() >=25);
	 passLoosePt = (currentLepton.Pt() >=10);
	 passSoftPt = (currentLepton.Pt() >=3);
	 if(passEta && passLooseId && passLoosePt && selMuons.size()==2) extraMuons.push_back(currentLepton); //No iso criteria for extra leptons.
	 if(passEta && !(passLooseId && passLoosePt) && passSoftId && passSoftPt && selMuons.size()==2) extraMuons.push_back(currentLepton); //Soft leptons. Need a particular cut?
         if(passEta && passIso && passId && passPt && selMuons.size()<2) selMuons.push_back(currentLepton);
      }
      for(int i = 0 ; i<PhotPt->size() ; i++){ //Photons
         TLorentzVector currentPhoton; currentPhoton.SetPtEtaPhiE(PhotPt->at(i),PhotEta->at(i),PhotPhi->at(i),PhotE3x3->at(i)); //photon energy???
	 //FIXME Criteria need come here...
	 selPhotons.push_back(currentPhoton);
      }
      for(int i =0 ; i<JetAk04Pt->size() ; i++){ //Jets
         bool passPt = false, passSelPt = false, passEta = false, passTightEta = false, passId = false, passLeptonCleaning = false, passPhotonCleaning = false;
	 TLorentzVector currentJet; currentJet.SetPtEtaPhiE(JetAk04Pt->at(i),JetAk04Eta->at(i),JetAk04Phi->at(i),JetAk04E->at(i));
	 passPt = (currentJet.Pt() >=15);
	 passSelPt = (currentJet.Pt() >=30);
	 double eta = fabs(currentJet.Eta());
	 passEta = (eta <= 4.7);
	 passTightEta = (eta <= 2.5);
	 if(eta<2.7){
	    passId = (JetAk04Id->at(i) >= 1); //This case is simple, it corresponds exactly to what we apply for now
	 }
	 float nef = JetAk04NeutralEmFrac->at(i);
	 float nhf = JetAk04NeutralHadAndHfFrac->at(i);
	 float nnp = 20; //FIXME this variable could not be found in the tree, it has obviously to be fixed if we are to reproduce what stands in the AN.
	 if(eta<3.0 && eta >=2.7){
	    //passId = (nef > 0.01 && nhf < 0.98 && nnp > 2);
	    passId = (JetAk04Id->at(i) >= 2); //Simpler criterium, but not equivalent to what is mentionned in the AN
	 }
	 if(eta>=3.0){
	    //passId = (nef<0.90 && nnp > 10);
	    passId = (JetAk04Id->at(i) >= 3); //Simpler criterium, but not equivalent to what is mentionned in the AN
	 }
	 double minDRlj(9999.); for(int ilep=0; ilep<selMuons.size(); ilep++) minDRlj = TMath::Min( minDRlj, utils::deltaR(currentJet,selMuons[ilep]) );
	 for(int ilep=0; ilep<selElectrons.size(); ilep++) minDRlj = TMath::Min( minDRlj, utils::deltaR(currentJet,selElectrons[ilep]) );
	 passLeptonCleaning = (minDRlj>=0.4);
	 passPhotonCleaning = true; //FIXME Needs implementation of the photon cleaning
	 if(passPt && passEta && passId && passLeptonCleaning && passPhotonCleaning) tagJets.push_back(currentJet);
	 if(passSelPt && passEta && passId && passLeptonCleaning && passPhotonCleaning) selJets.push_back(currentJet);
	 if(passSelPt && passTightEta && passId && passLeptonCleaning && passPhotonCleaning) btags.push_back(JetAk04BDiscCisvV2->at(i));
      }

      //Discriminate ee and mumu
      bool isEE = (selElectrons.size()==2); //2 good electrons
      bool isMuMu = (selMuons.size()==2); //2 good muons

      //Jet category
      enum {eq0jets,geq1jets,vbf};
      int jetCat = geq1jets;
      if(selJets.size()==0) jetCat = eq0jets;
      if(selJets.size()>=2){
         float etamin=0., etamax=0;
	 if(selJets[0].Eta()>selJets[1].Eta()) {etamax = selJets[0].Eta(); etamin = selJets[1].Eta();}
	 else {etamax = selJets[1].Eta(); etamin = selJets[0].Eta();}
	 bool centralJetVeto = true;
	 if(selJets.size()>2){
	    for(int i = 2 ; i < selJets.size() ; i++){
	       if(selJets[i].Eta()>etamin && selJets[i].Eta()<etamax) centralJetVeto = false;
	    }
	 }
	 bool passDeltaEta = (fabs(selJets[0].Eta() - selJets[1].Eta())>4);
	 bool passMjj = ((selJets[0]+selJets[1]).M()>500);
	 if(centralJetVeto && passDeltaEta && passMjj) jetCat = vbf;
      }
      mon.fillHisto("jetCategory","all",jetCat,weight);
      mon.fillHisto("nJets","all",selJets.size(),weight);

      mon.fillHisto("nb_mu","sel",selMuons.size(),weight);
      mon.fillHisto("nb_e","sel",selElectrons.size(),weight);
      mon.fillHisto("nb_mu","extra",extraMuons.size(),weight);
      mon.fillHisto("nb_e","extra",extraElectrons.size(),weight);

      //Begin selection

      if(!isEE && !isMuMu) continue; //not a good lepton pair
      mon.fillHisto("eventflow","all",1,weight);

      //Take care of the problematic case of 2 good electrons and 2 good muons
      if(isEE && isMuMu){
         if(fabs((selElectrons[0]+selElectrons[1]).M()-91.1876)<fabs((selMuons[0]+selMuons[1]).M()-91.1876)) isMuMu=false; //Closest one to Z mass
	 else isEE=false;
	 }
      
      vector<TLorentzVector> selLeptons;
      if(isEE) selLeptons = selElectrons;
      if(isMuMu) selLeptons = selMuons;
      TLorentzVector boson = selLeptons[0] + selLeptons[1];
      TLorentzVector METVector; METVector.SetPxPyPzE(METPx->at(0),METPy->at(0),METPz->at(0),METE->at(0));
      double transverseMass = sqrt(pow(sqrt(pow(boson.Pt(),2)+pow(boson.M(),2))+sqrt(pow(METVector.Pt(),2)+pow(91.1876,2)),2)-pow((boson+METVector).Pt(),2)); //Pretty long formula. Please check that it's correct.
      mon.fillHisto("M_Z","all",boson.M(),weight);
      mon.fillHisto("MET","all",METVector.Pt(),weight);
      mon.fillHisto("mT","all",transverseMass,weight);

      if(fabs(boson.M()-91.1876)>15.) continue;
      mon.fillHisto("eventflow","all",2,weight);

      if(boson.Pt() < 55.) continue;
      mon.fillHisto("eventflow","all",3,weight);

      if(isEE && extraElectrons.size()>0) continue;
      if(isMuMu && extraMuons.size()>0) continue;
      mon.fillHisto("eventflow","all",4,weight);

      //b veto
      bool passBTag = true;
      for(int i =0 ; i < btags.size() ; i++){
         if (btags[i] > 0.5426) passBTag = false;
      }
      if(!passBTag) continue;
      
      mon.fillHisto("eventflow","all",5,weight);
      
      //Phi(jet,MET)
      bool passDeltaPhiJetMET = true;
      for(int i = 0 ; i < tagJets.size() ; i++){
         if (fabs(tagJets[i].Phi()-METVector.Phi())<0.5) passDeltaPhiJetMET = false;
      }
      if(!passDeltaPhiJetMET) continue;
      
      mon.fillHisto("eventflow","all",6,weight);
      
      //Phi(Z,MET)
      double deltaPhiZMet = fabs(boson.Phi()-METVector.Phi());
      if(deltaPhiZMet<0.5) continue;
      mon.fillHisto("eventflow","all",7,weight);

      mon.fillHisto("M_Z","beforeMETcut",boson.M(),weight);
      mon.fillHisto("MET","beforeMETcut",METVector.Pt(),weight);
      mon.fillHisto("mT","beforeMETcut",transverseMass,weight);
      mon.fillHisto("jetCategory","beforeMETcut",jetCat,weight);
      mon.fillHisto("nJets","beforeMETcut",selJets.size(),weight);
      
      //MET>80
      if(METVector.Pt()<80) continue;
      mon.fillHisto("eventflow","all",8,weight);
      
      //MET>125
      if(METVector.Pt()<125) continue;
      mon.fillHisto("eventflow","all",9,weight);

      mon.fillHisto("M_Z","final",boson.M(),weight);
      mon.fillHisto("MET","final",METVector.Pt(),weight);
      mon.fillHisto("mT","final",transverseMass,weight);
      mon.fillHisto("jetCategory","final",jetCat,weight);
      mon.fillHisto("nJets","final",selJets.size(),weight);

      if(jetCat==vbf){
         mon.fillHisto("M_Z","final_vbf",boson.M(),weight);
         mon.fillHisto("MET","final_vbf",METVector.Pt(),weight);
         mon.fillHisto("mT","final_vbf",transverseMass,weight);
         mon.fillHisto("nJets","final_vbf",selJets.size(),weight);
      }
      if(jetCat==geq1jets){
         mon.fillHisto("M_Z","final_geq1jets",boson.M(),weight);
         mon.fillHisto("MET","final_geq1jets",METVector.Pt(),weight);
         mon.fillHisto("mT","final_geq1jets",transverseMass,weight);
         mon.fillHisto("nJets","final_geq1jets",selJets.size(),weight);
      }
      if(jetCat==eq0jets){
         mon.fillHisto("M_Z","final_eq0jets",boson.M(),weight);
         mon.fillHisto("MET","final_eq0jets",METVector.Pt(),weight);
         mon.fillHisto("mT","final_eq0jets",transverseMass,weight);
         mon.fillHisto("nJets","final_eq0jets",selJets.size(),weight);
      }
      


      //cout << "nb of interactions=" << EvtPuCnt << endl;

   }

   TFile* outFile=TFile::Open("out.root","recreate");
   mon.Write();
   outFile->Close();

}
