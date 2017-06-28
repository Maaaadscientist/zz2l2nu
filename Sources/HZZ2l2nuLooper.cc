#define HZZ2l2nuLooper_cxx
#include "../Includes/HZZ2l2nuLooper.h"
#include "../Includes/SmartSelectionMonitor.h"
#include "../Includes/Utils.h"
#include "../Includes/ObjectSelection.h"
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>

void HZZ2l2nuLooper::Loop()
{
   if (fChain == 0) return;

   //###############################################################
   //################## DECLARATION OF HISTOGRAMS ##################
   //###############################################################

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

   //###############################################################
   //##################     EVENT LOOP STARTS     ##################
   //###############################################################

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

     //###############################################################
     //##################     OBJECT SELECTION      ##################
     //###############################################################

      vector<TLorentzVector> selElectrons; //Leptons passing final cuts
      vector<TLorentzVector> selMuons; //Muons passing final cuts
      vector<TLorentzVector> extraElectrons; //Additional electrons, used for veto
      vector<TLorentzVector> extraMuons; //Additional muons, used for veto
      vector<TLorentzVector> selPhotons; //Photons
      vector<TLorentzVector> tagJets; //Jets passing Id and cleaning, with |eta|<4.7 and pT>15GeV. Used for the DeltaPhi cut.
      vector<TLorentzVector> selJets; //Jets passing Id and cleaning, with |eta|<4.7 and pT>30GeV. Used for jet categorization.
      vector<double> btags; //B-Tag discriminant, recorded for selJets with |eta|<2.5. Used for b-tag veto.

      objectSelection::selectElectrons(selElectrons, extraElectrons, ElPt, ElEta, ElPhi, ElE, ElId, ElEtaSc, ElPfIsoRho);
      objectSelection::selectMuons(selMuons, extraMuons, MuPt, MuEta, MuPhi, MuE, MuId, MuIdTight, MuPfIso);
      objectSelection::selectPhotons(selPhotons, PhotPt, PhotEta, PhotPhi, PhotScRawE, PhotId);
      objectSelection::selectJets(tagJets, selJets, btags, JetAk04Pt, JetAk04Eta, JetAk04Phi, JetAk04E, JetAk04Id, JetAk04NeutralEmFrac, JetAk04NeutralHadAndHfFrac, JetAk04BDiscCisvV2, selMuons, selElectrons, selPhotons);

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

     //###############################################################
     //##################       ANALYSIS CUTS       ##################
     //###############################################################

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

     //###############################################################
     //##################     END OF SELECTION      ##################
     //###############################################################

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
      



   }

   //###############################################################
   //##################        END OF LOOP        ##################
   //###############################################################

   TFile* outFile=TFile::Open("out.root","recreate");
   mon.Write();
   outFile->Close();

}
