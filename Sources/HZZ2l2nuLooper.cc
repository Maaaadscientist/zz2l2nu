#define HZZ2l2nuLooper_cxx
#include "../Includes/HZZ2l2nuLooper.h"
#include "../Includes/SmartSelectionMonitor.h"
#include "../Includes/SmartSelectionMonitor_hzz.h"
#include "../Includes/Utils.h"
#include "../Includes/ObjectSelection.h"
#include <ctime>
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

   SmartSelectionMonitor_hzz mon;
   mon.declareHistos();

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   cout << "nb of entries in the input file =" << nentries << endl;

   //###############################################################
   //##################     EVENT LOOP STARTS     ##################
   //###############################################################

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      if ((jentry>maxEvents_)&&(maxEvents_>=0)) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(jentry % 10000 ==0) cout << jentry << " of " << nentries << " it is now " << std::time(0) << endl;

      evt currentEvt;

      double weight = 1.;
      double totEventWeight = 1.;

      //get the MC event weight if exists
      if (isMC_) { 
        weight = (EvtWeights->size()>0 ? EvtWeights->at(0) : 1);
        if ((sumWeightInBonzai_>0)&&(sumWeightInBaobab_>0)) totEventWeight = weight*sumWeightInBaobab_/sumWeightInBonzai_;
      }
      else {
        totEventWeight = totalEventsInBaobab_/nentries;
      }

      mon.fillHisto("totEventInBaobab","tot",EvtPuCnt,totEventWeight);
      mon.fillHisto("eventflow","tot",0,weight);

      for(int i =0 ; i < MuPt->size() ; i++) mon.fillHisto("pT_mu","tot",MuPt->at(i),weight);
      for(int i =0 ; i < ElPt->size() ; i++) mon.fillHisto("pT_e","tot",ElPt->at(i),weight);
      mon.fillHisto("nb_mu","tot",MuPt->size(),weight);
      mon.fillHisto("nb_e","tot",ElPt->size(),weight);
      mon.fillHisto("pile-up","tot",EvtPuCnt,weight);

     //###############################################################
     //##################     OBJECT SELECTION      ##################
     //###############################################################

      vector<TLorentzVector> selElectrons; //Leptons passing final cuts
      vector<TLorentzVector> selMuons; //Muons passing final cuts
      vector<TLorentzVector> extraElectrons; //Additional electrons, used for veto
      vector<TLorentzVector> extraMuons; //Additional muons, used for veto
      vector<TLorentzVector> selPhotons; //Photons
      vector<TLorentzVector> selJets; //Jets passing Id and cleaning, with |eta|<4.7 and pT>30GeV. Used for jet categorization and deltaPhi cut.
      vector<double> btags; //B-Tag discriminant, recorded for selJets with |eta|<2.5. Used for b-tag veto.

      objectSelection::selectElectrons(selElectrons, extraElectrons, ElPt, ElEta, ElPhi, ElE, ElId, ElEtaSc, ElPfIsoRho);
      objectSelection::selectMuons(selMuons, extraMuons, MuPt, MuEta, MuPhi, MuE, MuId, MuIdTight, MuIdSoft, MuPfIso);
      objectSelection::selectPhotons(selPhotons, PhotPt, PhotEta, PhotPhi, PhotId, PhotScEta, PhotHasPixelSeed, selMuons, selElectrons);
      objectSelection::selectJets(selJets, btags, JetAk04Pt, JetAk04Eta, JetAk04Phi, JetAk04E, JetAk04Id, JetAk04NeutralEmFrac, JetAk04NeutralHadAndHfFrac, JetAk04NeutMult, JetAk04BDiscCisvV2, selMuons, selElectrons, selPhotons);

      //Discriminate ee and mumu
      bool isEE = (selElectrons.size()==2); //2 good electrons
      bool isMuMu = (selMuons.size()==2); //2 good muons

      mon.fillHisto("nb_mu","sel",selMuons.size(),weight);
      mon.fillHisto("nb_e","sel",selElectrons.size(),weight);
      mon.fillHisto("nb_mu","extra",extraMuons.size(),weight);
      mon.fillHisto("nb_e","extra",extraElectrons.size(),weight);

     //###############################################################
     //##################       ANALYSIS CUTS       ##################
     //###############################################################

      if(!isEE && !isMuMu) continue; //not a good lepton pair
      mon.fillHisto("eventflow","tot",1,weight);

      if(isEE) currentEvt.s_lepCat = "_ee";
      else currentEvt.s_lepCat = "_mumu";

      //Definition of the relevant analysis variables
      vector<TLorentzVector> selLeptons;
      if(isEE) selLeptons = selElectrons;
      if(isMuMu) selLeptons = selMuons;
      TLorentzVector boson = selLeptons[0] + selLeptons[1];
      TLorentzVector METVector; METVector.SetPxPyPzE(METPx->at(0),METPy->at(0),METPz->at(0),METE->at(0));
      currentEvt.transverseMass = sqrt(pow(sqrt(pow(boson.Pt(),2)+pow(boson.M(),2))+sqrt(pow(METVector.Pt(),2)+pow(91.1876,2)),2)-pow((boson+METVector).Pt(),2)); //Pretty long formula. Please check that it's correct.
      currentEvt.MZ = boson.M();
      currentEvt.pTZ = boson.Pt();
      currentEvt.MET = METVector.Pt();

      //Jet category
      enum {eq0jets,geq1jets,vbf};
      TString v_jetCat[3] = {"_eq0jets","_geq1jets","_vbf"};
      int jetCat = geq1jets;
      if(selJets.size()==0) jetCat = eq0jets;
      if(utils::passVBFcuts(selJets, boson)) jetCat = vbf;
      currentEvt.s_jetCat = v_jetCat[jetCat];
      mon.fillHisto("jetCategory","tot",jetCat,weight);
      currentEvt.nJets = selJets.size();
      mon.fillHisto("nJets","tot",currentEvt.nJets,weight);

      mon.fillAnalysisHistos(currentEvt, "tot", weight);

      if(fabs(boson.M()-91.1876)>15.) continue;
      mon.fillHisto("eventflow","tot",2,weight);

      if(boson.Pt() < 55.) continue;
      mon.fillHisto("eventflow","tot",3,weight);

      if(extraElectrons.size()>0 || extraMuons.size()>0) continue;
      if(isEE && selMuons.size()>0) continue;
      if(isMuMu && selElectrons.size()>0) continue;
      mon.fillHisto("eventflow","tot",4,weight);

      //b veto
      bool passBTag = true;
      for(int i =0 ; i < btags.size() ; i++){
         if (btags[i] > 0.5426) passBTag = false;
      }
      if(!passBTag) continue;

      mon.fillHisto("eventflow","tot",5,weight);

      //Phi(jet,MET)
      bool passDeltaPhiJetMET = true;
      for(int i = 0 ; i < selJets.size() ; i++){
         if (fabs(utils::deltaPhi(selJets[i], METVector))<0.5) passDeltaPhiJetMET = false;
      }
      if(!passDeltaPhiJetMET) continue;

      mon.fillHisto("eventflow","tot",6,weight);

      //Phi(Z,MET)
      double deltaPhiZMet = fabs(utils::deltaPhi(boson, METVector));
      if(deltaPhiZMet<0.5) continue;
      mon.fillHisto("eventflow","tot",7,weight);

      mon.fillAnalysisHistos(currentEvt, "beforeMETcut", weight);
      mon.fillHisto("jetCategory","beforeMETcut",jetCat,weight);

      //MET>80
      if(METVector.Pt()<80) continue;
      mon.fillHisto("eventflow","tot",8,weight);

      //MET>125
      if(METVector.Pt()<125) continue;
      mon.fillHisto("eventflow","tot",9,weight);

     //###############################################################
     //##################     END OF SELECTION      ##################
     //###############################################################

      mon.fillHisto("jetCategory","final",jetCat,weight);
      mon.fillAnalysisHistos(currentEvt, "final", weight);

   }

   //###############################################################
   //##################        END OF LOOP        ##################
   //###############################################################

   TFile* outFile=TFile::Open(outputFile_,"recreate");
   mon.Write();
   outFile->Close();

}
