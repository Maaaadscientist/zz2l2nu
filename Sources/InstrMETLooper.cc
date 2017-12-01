#define InstrMETLooper_cxx
#include "../Includes/LooperMain.h"
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
#include "../Includes/PhotonEfficiencySF.h"

void LooperMain::Loop_InstrMET()
{
  std::cout<<"Starting the InstrMET Looper..."<<std::endl;
  if (fChain == 0) return;

  //###############################################################
  //################## DECLARATION OF HISTOGRAMS ##################
  //###############################################################


  SmartSelectionMonitor_hzz mon;
  mon.declareHistos_InstrMET();

  Long64_t nentries = fChain->GetEntries();

  TString fileName = fChain->GetCurrentFile()->GetName();

  bool isMC_QCD = (isMC_ && fileName.Contains("QCD"));
  bool isMC_GJet = (isMC_ && fileName.Contains("GJet"));
  bool isMC_Wlnu_inclusive = (isMC_ && fileName.Contains("_WJets_") && !fileName.Contains("HT"));
  bool isMC_Wlnu_HT100 = (isMC_ && fileName.Contains("_WJets_HT-") );
  bool isMC_WGToLNuG = (isMC_ && fileName.Contains("WGToLNuG") );
  bool isMC_ZNuNuGJets = (isMC_ && fileName.Contains("ZNuNuGJets"));
  bool isMC_ZJetsToNuNu = (isMC_ && fileName.Contains("ZJetsToNuNu"));



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
    int eventflowStep=0;

    //get the MC event weight if exists
    if (isMC_) { 
      weight = (EvtWeights->size()>0 ? EvtWeights->at(0) : 1);
      if ((sumWeightInBonzai_>0)&&(sumWeightInBaobab_>0)) totEventWeight = weight*sumWeightInBaobab_/sumWeightInBonzai_;
    }
    else {
      totEventWeight = totalEventsInBaobab_/nentries;
    }

    mon.fillHisto("totEventInBaobab","tot",EvtPuCnt,totEventWeight);
    mon.fillHisto("eventflow","tot",eventflowStep++,weight); //output of bonzais

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

    objectSelection::selectElectrons(selElectrons, extraElectrons, ElPt, ElEta, ElPhi, ElE, ElId, ElEtaSc);
    objectSelection::selectMuons(selMuons, extraMuons, MuPt, MuEta, MuPhi, MuE, MuId, MuIdTight, MuIdSoft, MuPfIso);
    objectSelection::selectPhotons(selPhotons, PhotPt, PhotEta, PhotPhi, PhotId, PhotScEta, PhotHasPixelSeed, selMuons, selElectrons);
    objectSelection::selectJets(selJets, btags, JetAk04Pt, JetAk04Eta, JetAk04Phi, JetAk04E, JetAk04Id, JetAk04NeutralEmFrac, JetAk04NeutralHadAndHfFrac, JetAk04NeutMult, JetAk04BDiscCisvV2, selMuons, selElectrons, selPhotons);

    //Ask for a prompt photon
    if(selPhotons.size() != 1) continue;

    //Apply prescales here
    //for the moment this function just return true or false, next it will return 0 or the prescale
    int triggerWeight =0;
    triggerWeight = utils::passTrigger(utils::SinglePhoton, TrigHltDiMu, TrigHltMu, TrigHltDiEl, TrigHltEl, TrigHltElMu, TrigHltPhot, TrigHltDiMu_prescale, TrigHltMu_prescale, TrigHltDiEl_prescale, TrigHltEl_prescale, TrigHltElMu_prescale, TrigHltPhot_prescale, selPhotons[0].Pt());
    if(triggerWeight==0) continue; //trigger not found
    mon.fillHisto("pT_Z","noPrescale",selPhotons[0].Pt(),weight);
    weight *= triggerWeight;
    mon.fillHisto("pT_Z","withPrescale",selPhotons[0].Pt(),weight);
    mon.fillHisto("eventflow","tot",eventflowStep++,weight); //after prescale

    //photon efficiencies
    PhotonEfficiencySF phoEff;
    if(isMC_) weight *= phoEff.getPhotonEfficiency(selPhotons[0].Pt(), selPhotons[0].Eta(), "tight",utils::CutVersion::ICHEP16Cut ).first; //NB: By definition of selPhotons, Eta is in fact the supercluster eta.
    mon.fillHisto("pT_Z","withPrescale_and_phoEff",selPhotons[0].Pt(),weight);

    mon.fillHisto("eventflow","tot",eventflowStep++,weight); //after Pho eff

    if(isMC_){
      //get the PU weights
      float weightPU = pileUpWeight(EvtPuCntTruth);
      weight *= weightPU;
    }
    mon.fillHisto("eventflow","tot",eventflowStep++,weight); //after PU reweighting

    std::vector<std::pair<int, int> > listMETFilter; //after the passMetFilter function, it contains the bin number of the cut in .first and if it passed 1 or not 0 the METfilter
    bool passMetFilter = utils::passMetFilter(TrigMET, listMETFilter, isMC_);
    //now fill the metFilter eventflow
    mon.fillHisto("metFilters","tot",26,weight); //the all bin, i.e. the last one
    for(unsigned int i =0; i < listMETFilter.size(); i++){
      if(listMETFilter[i].second ==1) mon.fillHisto("metFilters","tot",listMETFilter[i].first,weight);
    }

    if (!passMetFilter) continue;
    mon.fillHisto("eventflow","tot",eventflowStep++,weight); // after met filters

    //Resolve G+jet/QCD mixing (avoid double counting of photons)
    if (isMC_GJet || isMC_QCD ||
        isMC_Wlnu_inclusive || isMC_Wlnu_HT100 || isMC_WGToLNuG ||
        isMC_ZNuNuGJets || isMC_ZJetsToNuNu ) {
      // iF GJet sample; accept only event with prompt photons
      // if QCD sample; reject events with prompt photons in final state
      bool gPromptFound=false;
      //for(unsigned int i = 0; i < GPhotMotherId->size(); i++) std::cout<< GPhotMotherId->at(i) << std::endl;
      /*for(size_t ig=0; ig<gen.size(); ig++){
        if((abs(gen[ig].pdgId())==22) && gen[ig].isPromptFinalState())  gPromptFound=true;
      }*/
      if ( (isMC_GJet) && (!gPromptFound) ) continue; //reject event
      if ( (isMC_QCD) && gPromptFound ) continue; //reject event
      if ( ( isMC_Wlnu_inclusive || isMC_Wlnu_HT100) && gPromptFound ) continue;
      //             if ( (isMC_WGToLNuG) && (!gPromptFound) ) continue;
      // if ( (isMC_ZNuNuGJets) && (!gPromptFound) ) continue;
      if ( (isMC_ZJetsToNuNu) && gPromptFound ) continue;
    }



    for(int i =0 ; i < MuPt->size() ; i++) mon.fillHisto("pT_mu","afterWeight",MuPt->at(i),weight);
    for(int i =0 ; i < ElPt->size() ; i++) mon.fillHisto("pT_e","afterWeight",ElPt->at(i),weight);
    for(int i =0 ; i < PhotPt->size() ; i++) mon.fillHisto("pT_Z","afterWeight",PhotPt->at(i),weight);
    mon.fillHisto("nb_pho","afterWeight",PhotPt->size(),weight);
    mon.fillHisto("pile-up","afterWeight",EvtPuCnt,weight);
    mon.fillHisto("truth-pile-up","afterWeight",EvtPuCntTruth,weight);
    mon.fillHisto("reco-vtx","afterWeight",EvtVtxCnt,weight);

    mon.fillHisto("nb_mu","sel_afterWeight",selMuons.size(),weight);
    mon.fillHisto("nb_e","sel_afterWeight",selElectrons.size(),weight);
    mon.fillHisto("nb_pho","sel_afterWeight",selPhotons.size(),weight);
    mon.fillHisto("nb_mu","extra_afterWeight",extraMuons.size(),weight);
    mon.fillHisto("nb_e","extra_afterWeight",extraElectrons.size(),weight);

    //###############################################################
    //##################       ANALYSIS CUTS       ##################
    //###############################################################


    currentEvt.s_lepCat = "_gamma";

    //Definition of the relevant analysis variables
    TLorentzVector boson = selPhotons[0];
    //generate mass from the line shape:
    //Will read the mass line shape of the ee and mumu case
    TLorentzVector METVector; METVector.SetPxPyPzE(METPx->at(0),METPy->at(0),METPz->at(0),METE->at(0));
    currentEvt.transverseMass = sqrt(pow(sqrt(pow(boson.Pt(),2)+pow(boson.M(),2))+sqrt(pow(METVector.Pt(),2)+pow(91.1876,2)),2)-pow((boson+METVector).Pt(),2)); //Pretty long formula. Please check that it's correct.
    currentEvt.MZ = boson.M();
    currentEvt.pTZ = boson.Pt();
    currentEvt.MET = METVector.Pt();
    currentEvt.etaZ = boson.Eta();

    //Jet category
    enum {eq0jets,geq1jets,vbf};
    TString v_jetCat[3] = {"_eq0jets","_geq1jets","_vbf"};
    int jetCat = geq1jets;
    if(selJets.size()==0) jetCat = eq0jets;
    if(utils::passVBFcuts(selJets, boson)) jetCat = vbf;
    currentEvt.s_jetCat = v_jetCat[jetCat];
    mon.fillHisto("jetCategory","afterWeight",jetCat,weight);
    currentEvt.nJets = selJets.size();
    mon.fillHisto("nJets","afterWeight",currentEvt.nJets,weight);

    mon.fillAnalysisHistos_InstrMET(currentEvt, "afterWeight", weight);

    if(boson.Pt() < 55.) continue;
    mon.fillHisto("eventflow","tot",eventflowStep++,weight); //after pt cut

    if(selElectrons.size()+extraElectrons.size()+selMuons.size()+extraMuons.size()>0) continue;
    mon.fillHisto("eventflow","tot",eventflowStep++,weight); //after no extra leptons

    //b veto
    bool passBTag = true;
    for(int i =0 ; i < btags.size() ; i++){
      if (btags[i] > 0.5426) passBTag = false;
    }
    if(!passBTag) continue;

    mon.fillHisto("eventflow","tot",eventflowStep++,weight); //after b-tag veto

    //Phi(jet,MET)
    bool passDeltaPhiJetMET = true;
    for(int i = 0 ; i < selJets.size() ; i++){
      if (fabs(utils::deltaPhi(selJets[i], METVector))<0.5) passDeltaPhiJetMET = false;
    }
    if(!passDeltaPhiJetMET) continue;

    mon.fillHisto("eventflow","tot",eventflowStep++,weight); //after delta phi (jet, met)

    //Phi(Z,MET)
    double deltaPhiZMet = fabs(utils::deltaPhi(boson, METVector));
    if(deltaPhiZMet<0.5) continue;
    mon.fillHisto("eventflow","tot",eventflowStep++,weight); //after delta phi (Z, met)

    mon.fillAnalysisHistos(currentEvt, "beforeMETcut", weight);
    mon.fillHisto("reco-vtx","beforeMETcut",EvtVtxCnt,weight);
    mon.fillHisto("jetCategory","beforeMETcut",jetCat,weight);


    //MET>80
    if(METVector.Pt()<80) continue;
    mon.fillHisto("eventflow","tot",eventflowStep++,weight); //after MET > 80

    //MET>125
    if(METVector.Pt()<125) continue;
    mon.fillHisto("eventflow","tot",eventflowStep++,weight); //after MET > 125

    //###############################################################
    //##################     END OF SELECTION      ##################
    //###############################################################
    mon.fillHisto("reco-vtx","final",EvtVtxCnt,weight);
    mon.fillHisto("jetCategory","final",jetCat,weight);
    mon.fillAnalysisHistos_InstrMET(currentEvt, "final", weight);

  }

  //###############################################################
  //##################        END OF LOOP        ##################
  //###############################################################

  TFile* outFile=TFile::Open(outputFile_,"recreate");
  mon.Write();
  outFile->Close();

}
