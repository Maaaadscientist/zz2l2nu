#define InstrMETLooper_cxx
#include "../Loopers_Includes/LooperMain.h"
#include "../Loopers_Includes/SmartSelectionMonitor.h"
#include "../Loopers_Includes/SmartSelectionMonitor_hzz.h"
#include "../Common/Utils.h"
#include "../Common/ObjectSelection.h"
#include "../Common/PhotonEfficiencySF.h"
#include "../Common/Trigger.h"
#include "../Common/TLorentzVectorWithIndex.h"
#include <ctime>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

#define MAXIMAL_AMOUNT_OF_HISTOS false //Set this to FALSE to only have the main histos. If you want all debug histo at different level, set it to TRUE

void LooperMain::Loop_InstrMET()
{
  std::cout<<"Starting the InstrMET Looper..."<<std::endl;
  if (fChain == 0) return;

  //###############################################################
  //################## DECLARATION OF HISTOGRAMS ##################
  //###############################################################


  SmartSelectionMonitor_hzz mon;
  mon.declareHistos_InstrMET();

  Long64_t nentries = fChain->GetEntries(); //Warning: GetEntries has to be called before any other work on fChain

  TString fileName = fChain->GetCurrentFile()->GetName();

  bool isMC_QCD = (isMC_ && fileName.Contains("-QCD_"));
  bool isMC_QCD_HT = (isMC_ && fileName.Contains("-QCD_") && fileName.Contains("HT"));
  bool isMC_GJet = (isMC_ && fileName.Contains("-GJets_"));
  bool isMC_GJet_HT = (isMC_ && fileName.Contains("-GJets_HT"));
  bool isMC_Wlnu_inclusive = (isMC_ && fileName.Contains("-WJetsToLNu_") && !fileName.Contains("HT"));
  bool isMC_Wlnu_HT100 = (isMC_ && fileName.Contains("-WJetsToLNu_HT-") );
  bool isMC_WGToLNuG = (isMC_ && fileName.Contains("-WGToLNuG_") );
  bool isMC_LO_ZNuNuGJets = (isMC_ && fileName.Contains("-ZNuNuGJets_"));
  bool isMC_NLO_ZGTo2NuG_inclusive = (isMC_ && fileName.Contains("-ZGTo2NuG_") && !fileName.Contains("PtG-130"));
  bool isMC_NLO_ZGTo2NuG_Pt130 = (isMC_ && fileName.Contains("-ZGTo2NuG_PtG-130_"));
  bool isMC_ZJetsToNuNu = (isMC_ && fileName.Contains("-ZJetsToNuNu_"));

  Long64_t nbytes = 0, nb = 0;
  cout << "nb of entries in the input file " << fileName << " = " << nentries << endl;

  //Compute once weights for Instr. MET reweighting if needed
  std::vector<string> v_jetCat = {"_eq0jets","_geq1jets","_vbf"};
  string const base_path = string(getenv("HZZ2L2NU_BASE")) + "/";
  bool doClosureTest = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/please_do_closure_test_when_running_InstrMETLooper");
  if(doClosureTest) std::cout << "/!\\/!\\ CLOSURE TEST ONGOING - not wanted? Then remove 'WeightsAndDatadriven/InstrMET/please_do_closure_test_when_running_InstrMETLooper' /!\\/!\\" << std::endl;
  std::string weightFileType = (doClosureTest) ? "closureTest" : "InstrMET";
  bool weight_NVtx_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root");
  bool weight_Pt_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root");
  bool weight_Mass_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_lineshape_mass.root");
  std::map<TString, std::map<double, std::pair<double, double> > > NVtxWeight_map, PtWeight_map;
  std::map<TString, TH1D*> LineshapeMassWeight_map;
  utils::loadInstrMETWeights(weight_NVtx_exist, weight_Pt_exist, weight_Mass_exist, NVtxWeight_map, PtWeight_map, LineshapeMassWeight_map, weightFileType, base_path, v_jetCat);
  std::vector<TString> tagsR;
  tagsR.push_back("_gamma"); //_gamma, i.e. no reweighting to ee or mumu
  if(weight_NVtx_exist){
    tagsR.push_back("_ee");
    tagsR.push_back("_mumu");
    tagsR.push_back("_ll");
  }
  unsigned int tagsR_size =  tagsR.size();

  bool divideFinalHistoByBinWidth = false; //For final plots, we don't divide by the bin width to ease computations of the yields by eye.

  //###############################################################
  //##################     EVENT LOOP STARTS     ##################
  //###############################################################

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if ((jentry>maxEvents_)&&(maxEvents_>=0)) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(jentry % 10000 ==0) cout << jentry << " of " << nentries << " it is now " << std::time(0) << endl;
    photon_evt currentEvt;

    double weight = 1.;
    double totEventWeight = 1.;
    int eventflowStep=0;

    //get the MC event weight if exists
    if (isMC_) { 
      weight = (EvtWeights->size()>0 ? EvtWeights->at(1) : 1); //Value 0 is not filled properly for LO generated samples (MadgraphMLM)
      if ((sumWeightInBonzai_>0)&&(sumWeightInBaobab_>0)) totEventWeight = weight*sumWeightInBaobab_/sumWeightInBonzai_;
      if (jentry == 0){
        std::cout<< "Printing once the content of EvtWeights for event " << jentry << ":" << std::endl;
        for(unsigned int i = 0; i < EvtWeights->size(); i++ ) std::cout<< i << " " << EvtWeights->at(i) << std::endl;
      }
    }
    else {
      totEventWeight = totalEventsInBaobab_/nentries;
    }


    mon.fillHisto("totEventInBaobab","tot",EvtPuCnt,totEventWeight);
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //output of bonzais
    eventflowStep++;

    //Cleaning of low stats MC spikes that are gathering in some specific spot (in MET phi, pt, MET delta phi(MET, spike)...).
    bool isPathologicEvent=false;
    if(isMC_) isPathologicEvent = objectSelection::cleanPathologicEventsInPhotons(fileName, EvtRunNum, EvtLumiNum, EvtNum);
    if(isPathologicEvent) continue;

    // Remove events with 0 vtx
    if(EvtVtxCnt == 0 ) continue;

    //###############################################################
    //##################     OBJECT SELECTION      ##################
    //###############################################################

    vector<TLorentzVectorWithIndex> selElectrons; //Leptons passing final cuts
    vector<TLorentzVectorWithIndex> selMuons; //Muons passing final cuts
    vector<TLorentzVectorWithIndex> extraElectrons; //Additional electrons, used for veto
    vector<TLorentzVectorWithIndex> extraMuons; //Additional muons, used for veto
    vector<TLorentzVectorWithIndex> selPhotons; //Photons
    vector<TLorentzVectorWithIndex> selJets, selCentralJets; //Jets passing Id and cleaning, with |eta|<4.7 and pT>30GeV. Used for jet categorization and deltaPhi cut.
    vector<double> btags; //B-Tag discriminant, recorded for selJets with |eta|<2.5. Used for b-tag veto.

    objectSelection::selectElectrons(selElectrons, extraElectrons, ElPt, ElEta, ElPhi, ElE, ElId, ElEtaSc);
    objectSelection::selectMuons(selMuons, extraMuons, MuPt, MuEta, MuPhi, MuE, MuId, MuIdTight, MuIdSoft, MuPfIso);
    objectSelection::selectPhotons(selPhotons, PhotPt, PhotEta, PhotPhi, PhotId, PhotScEta, PhotHasPixelSeed, PhotSigmaIetaIeta, PhotSigmaIphiIphi, selMuons, selElectrons);
    objectSelection::selectJets(selJets, selCentralJets, btags, JetAk04Pt, JetAk04Eta, JetAk04Phi, JetAk04E, JetAk04Id, JetAk04NeutralEmFrac, JetAk04NeutralHadAndHfFrac, JetAk04NeutMult, JetAk04BDiscCisvV2, selMuons, selElectrons, selPhotons);

    //Ask for a prompt photon
    if(selPhotons.size() != 1) continue;

    //Check trigger and find prescale
    int triggerWeight =0;
    int triggerType;
    if(isMC_) triggerType = trigger::MC_Photon;
    else triggerType = trigger::SinglePhoton;

    triggerWeight = trigger::passTrigger(triggerType, TrigHltDiMu, TrigHltMu, TrigHltDiEl, TrigHltEl, TrigHltElMu, TrigHltPhot, TrigHltDiMu_prescale, TrigHltMu_prescale, TrigHltDiEl_prescale, TrigHltEl_prescale, TrigHltElMu_prescale, TrigHltPhot_prescale, selPhotons[0].Pt());
    if(triggerWeight==0) continue; //trigger not found

    if(MAXIMAL_AMOUNT_OF_HISTOS) mon.fillHisto("pT_Boson","noPrescale",selPhotons[0].Pt(),weight);
    weight *= triggerWeight;
    if(MAXIMAL_AMOUNT_OF_HISTOS) mon.fillHisto("pT_Boson","withPrescale",selPhotons[0].Pt(),weight);
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after prescale
    eventflowStep++;

    //photon efficiencies
    PhotonEfficiencySF phoEff;
    if(isMC_) weight *= phoEff.getPhotonEfficiency(selPhotons[0].Pt(), PhotScEta->at(selPhotons[0].GetIndex()), "tight",utils::CutVersion::Moriond17Cut ).first; 
    if(MAXIMAL_AMOUNT_OF_HISTOS) mon.fillHisto("pT_Boson","withPrescale_and_phoEff",selPhotons[0].Pt(),weight);

    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after Photon Efficiency
    eventflowStep++;

    float weightPU =1.;
    if(isMC_){
      //get the PU weights
      weightPU = pileUpWeight(EvtPuCntTruth); //on full 2016 data
    }
    weight *= weightPU;
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after PU reweighting
    eventflowStep++;

    std::vector<std::pair<int, int> > listMETFilter; //after the passMetFilter function, it contains the bin number of the cut in .first and if it passed 1 or not 0 the METfilter
    bool passMetFilter = utils::passMetFilter(TrigMET, listMETFilter, isMC_);
    //now fill the metFilter eventflow
    mon.fillHisto("metFilters","tot",26,weight); //the all bin, i.e. the last one
    for(unsigned int i =0; i < listMETFilter.size(); i++){
      if(listMETFilter[i].second ==1 && MAXIMAL_AMOUNT_OF_HISTOS) mon.fillHisto("metFilters","tot",listMETFilter[i].first,weight);
    }

    if (!passMetFilter) continue;
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after MET filters
    eventflowStep++;


    //Resolve G+jet/QCD mixing (avoid double counting of photons)
    bool passVetoQCDevent = true;
    if(isMC_QCD){
      for(unsigned int i = 0; i < GPhotPrompt->size(); i++){
        if(GPhotPrompt->at(i) && GPhotPt->at(i) >25){ //Gjets generated prompt photon above 25 GeV. QCD above 10 GeV, so the double counting occurs above 25.
          passVetoQCDevent = false;
        }
      }
    }
    if(!passVetoQCDevent) continue;

    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after avoiding G+jets and QCD mixing
    eventflowStep++;

    //LO to NLO k-factor for ZNuNuGamma (ref: fig 16 (bottom right) of http://link.springer.com/article/10.1007%2FJHEP02%282016%29057)
    double kFactor_ZNuNuGWeight = 1.;
    if(isMC_LO_ZNuNuGJets){
      //reconstruct the gen transverse energy
      std::vector<TLorentzVector> genNeutrinosFromZ;
      TLorentzVector tmpVector;
      for (unsigned int i =0; i < GLepBarePt->size(); i++){
        if(fabs(GLepBareId->at(i))==12 || fabs(GLepBareId->at(i))==14 || fabs(GLepBareId->at(i))==16){
          if(fabs(GLepBareMomId->at(i))==23 /*&& genParticle.mother()->status()==62*/){ //after testing, the status is not needed at all.
            tmpVector.SetPtEtaPhiE(GLepBarePt->at(i), GLepBareEta->at(i), GLepBarePhi->at(i), GLepBareE->at(i));
            genNeutrinosFromZ.push_back(tmpVector);//neutrino originating directly from Z boson
          }
        }
      }
      if(genNeutrinosFromZ.size() < 2) continue;
      TLorentzVector genZnunuBoson;
      genZnunuBoson = genNeutrinosFromZ[0] + genNeutrinosFromZ[1]; //Z from neutrinos at gen lvl

      //Apply LO to NLO k-factor for ZNuNuGamma (ref: fig 16 (bottom right) of http://link.springer.com/article/10.1007%2FJHEP02%282016%29057)
      if(      genZnunuBoson.Pt() > 960 ) kFactor_ZNuNuGWeight = 2.05;
      else if( genZnunuBoson.Pt() > 920 ) kFactor_ZNuNuGWeight = 2.10;
      else if( genZnunuBoson.Pt() > 880 ) kFactor_ZNuNuGWeight = 2.13;
      else if( genZnunuBoson.Pt() > 800 ) kFactor_ZNuNuGWeight = 2.16;
      else if( genZnunuBoson.Pt() > 440 ) kFactor_ZNuNuGWeight = 2.20;
      else if( genZnunuBoson.Pt() > 400 ) kFactor_ZNuNuGWeight = 2.16;
      else if( genZnunuBoson.Pt() > 360 ) kFactor_ZNuNuGWeight = 2.13;
      else if( genZnunuBoson.Pt() > 320 ) kFactor_ZNuNuGWeight = 2.07;
      else if( genZnunuBoson.Pt() > 280 ) kFactor_ZNuNuGWeight = 2.03;
      else if( genZnunuBoson.Pt() > 240 ) kFactor_ZNuNuGWeight = 1.96;
      else if( genZnunuBoson.Pt() > 200 ) kFactor_ZNuNuGWeight = 1.90;
      else if( genZnunuBoson.Pt() > 160 ) kFactor_ZNuNuGWeight = 1.75;
      else if( genZnunuBoson.Pt() > 120 ) kFactor_ZNuNuGWeight = 1.50;
      else if( genZnunuBoson.Pt() > 100 ) kFactor_ZNuNuGWeight = 1.32;
      else kFactor_ZNuNuGWeight = 1.;
    }
    weight *= kFactor_ZNuNuGWeight;

    if( isMC_NLO_ZGTo2NuG_inclusive && selPhotons[0].Pt() >= 130) continue;
    if( isMC_NLO_ZGTo2NuG_Pt130 && selPhotons[0].Pt() < 130) continue;

    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after LO-to-NLO k-factor for ZnunuGamma
    eventflowStep++;

    //Avoid double couting for W+jets
    if (isMC_Wlnu_inclusive || isMC_Wlnu_HT100){ //Avoid double counting and make our W#rightarrow l#nu exclusif of the dataset with a cut on HT...
      bool isHT100 = false;

      //Let's create our own HT variable
      double vHT =0;
      TLorentzVector genJet_uncleaned;
      for(size_t ig=0; ig<GJetAk04Pt->size(); ig++){
        genJet_uncleaned.SetPtEtaPhiE(GJetAk04Pt->at(ig), GJetAk04Eta->at(ig), GJetAk04Phi->at(ig), GJetAk04E->at(ig));
        //cross-clean with selected leptons and photons
        double minDRmj(9999.); for(size_t ilepM=0; ilepM<selMuons.size();     ilepM++)  minDRmj = TMath::Min( minDRmj, utils::deltaR(genJet_uncleaned,selMuons[ilepM]) );
        double minDRej(9999.); for(size_t ilepE=0; ilepE<selElectrons.size(); ilepE++)  minDRej = TMath::Min( minDRej, utils::deltaR(genJet_uncleaned,selElectrons[ilepE]) );
        //double minDRgj(9999.); for(size_t ipho=0;  ipho <selPhotons.size();   ipho++)   minDRgj = TMath::Min( minDRgj, utils::deltaR(genJet_uncleaned,selPhotons[ipho]) );
        //if(minDRmj<0.4 || minDRej<0.4 || minDRgj<0.4) continue;
        if(minDRmj<0.4 || minDRej<0.4) continue;

        vHT += GJetAk04Pt->at(ig);
      }
      if(vHT >100) isHT100 = true;
      if(MAXIMAL_AMOUNT_OF_HISTOS){
        if(isMC_Wlnu_inclusive) mon.fillHisto("custom_HT","forWlnu_inclusive",vHT,weight);
        if(isMC_Wlnu_HT100) mon.fillHisto("custom_HT","forWlnu_HT100",vHT,weight);
      }
      if(isMC_Wlnu_inclusive && isHT100) continue; //reject event
      if(isMC_Wlnu_HT100 && !isHT100) continue; //reject event

    }

    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after avoiding double counting for W+jets
    eventflowStep++;

    //###############################################################
    //##################       ANALYSIS CUTS       ##################
    //###############################################################


    //Definition of the relevant analysis variables and storage in the currentEvt
    TLorentzVector boson = selPhotons[0];
    TLorentzVector METVector; METVector.SetPtEtaPhiE(METPtType1XY->at(0),0.,METPhiType1XY->at(0),METPtType1XY->at(0));

    //Jet category
    enum {eq0jets,geq1jets,vbf};
    int jetCat = geq1jets;
    if(selJets.size()==0) jetCat = eq0jets;
    if(utils::passVBFcuts(selJets, boson)) jetCat = vbf;


    currentEvt.Fill_photonEvt(v_jetCat[jetCat], tagsR[0], boson, METVector, selJets, EvtRunNum, EvtVtxCnt, EvtFastJetRho, METsig->at(0), PhotHoE->at(selPhotons[0].GetIndex()), PhotSigmaIetaIeta->at(selPhotons[0].GetIndex()), utils::photon_rhoCorrectedIso(PhotPfIsoChHad->at(selPhotons[0].GetIndex()), EvtFastJetRho, PhotScEta->at(selPhotons[0].GetIndex()), "chIso"), utils::photon_rhoCorrectedIso(PhotPfIsoNeutralHad->at(selPhotons[0].GetIndex()), EvtFastJetRho, PhotScEta->at(selPhotons[0].GetIndex()), "nhIso"), utils::photon_rhoCorrectedIso(PhotPfIsoPhot->at(selPhotons[0].GetIndex()), EvtFastJetRho, PhotScEta->at(selPhotons[0].GetIndex()), "gIso"), PhotR9->at(selPhotons[0].GetIndex())); 

    //PUPPI variables
    TLorentzVector PUPPIMETVector; PUPPIMETVector.SetPtEtaPhiE(METPtType1XY->at(2),0.,METPhiType1XY->at(2),METPtType1XY->at(2));
    double transverseMass_PUPPI = sqrt(pow(sqrt(pow(boson.Pt(),2)+pow(boson.M(),2))+sqrt(pow(PUPPIMETVector.Pt(),2)+pow(91.1876,2)),2)-pow((boson+PUPPIMETVector).Pt(),2));

    // compute the parallele and the orthogonal MET
    double PUPPIMETorth = 0;
    double PUPPIMETpar = 0;
    TVector3 bosonDir = boson.Vect();
    if (bosonDir.Perp()>0){
      TVector3 bosonDir2D(bosonDir.x()/bosonDir.Perp(), bosonDir.y()/bosonDir.Perp(), 0);
      TVector3 bosonPerp2D(-bosonDir2D.y(),bosonDir2D.x(),0);
      TVector3 PUPPIMETVector3D = PUPPIMETVector.Vect();
      PUPPIMETpar = - PUPPIMETVector3D.Dot(bosonDir2D);
      PUPPIMETorth = PUPPIMETVector3D.Dot(bosonPerp2D);
    }

    mon.fillHisto("jetCategory","afterWeight",jetCat,weight);
    mon.fillAnalysisHistos(currentEvt, "afterWeight", weight);

    //std::cout<<"Event info: " << EvtRunNum<<":"<<EvtLumiNum<<":"<<EvtNum << "; boson pt = "<<boson.Pt()<<"; weight = "<<weight<<"; triggerPrescale = "<<triggerWeight<<"; met = "<<currentEvt.MET<<"; mt = "<<currentEvt.MT<<"; njets = "<<currentEvt.nJets<<"; vtx = "<<EvtVtxCnt<<"; rho = "<<EvtFastJetRho<<"; puWeight = "<<weightPU<<std::endl;

    //More MET variables
    // -- PUPPI MET
    mon.fillHisto("MET", "afterWeight_PUPPI"+currentEvt.s_jetCat, PUPPIMETVector.Pt(), weight, true);
    mon.fillHisto("MET_phi", "afterWeight_PUPPI"+currentEvt.s_jetCat, PUPPIMETVector.Phi(), weight, true);
    mon.fillHisto("DeltaPhi_MET_Boson", "afterWeight_PUPPI"+currentEvt.s_jetCat, fabs(utils::deltaPhi(boson, PUPPIMETVector)), weight);
    double minDeltaPhiJetMET_PUPPI = 4.;
    for(int i = 0 ; i < selJets.size() ; i++){
      if (fabs(utils::deltaPhi(selJets[i], PUPPIMETVector)) < minDeltaPhiJetMET_PUPPI) minDeltaPhiJetMET_PUPPI = fabs(utils::deltaPhi(selJets[i], PUPPIMETVector));
    }
    mon.fillHisto("DeltaPhi_MET_Jet","afterWeight_PUPPI"+currentEvt.s_jetCat,minDeltaPhiJetMET_PUPPI,weight);
    mon.fillHisto("mT",  "afterWeight_PUPPI"+currentEvt.s_jetCat, transverseMass_PUPPI,       weight, true);
    // -- MET significance (PUPPI MET)
    mon.fillHisto("METsigx2", "afterWeight_PUPPI"+currentEvt.s_jetCat, METsigx2->at(2), weight);
    mon.fillHisto("METsigxy", "afterWeight_PUPPI"+currentEvt.s_jetCat, METsigxy->at(2), weight);
    mon.fillHisto("METsigy2", "afterWeight_PUPPI"+currentEvt.s_jetCat, METsigy2->at(2), weight);
    mon.fillHisto("METsig", "afterWeight_PUPPI"+currentEvt.s_jetCat, METsig->at(2), weight);
    //MET/pt (PF MET and PUPPI MET)
    mon.fillHisto("METoverPt", "afterWeight_PUPPI"+currentEvt.s_jetCat, PUPPIMETVector.Pt()/(1.*boson.Pt()), weight);
    mon.fillHisto("METoverPt_zoom", "afterWeight_PUPPI"+currentEvt.s_jetCat, PUPPIMETVector.Pt()/(1.*boson.Pt()), weight);

    //parallele and the orthogonal PUPPI MET
    mon.fillHisto("METperp", "afterWeight_PUPPI"+currentEvt.s_jetCat, PUPPIMETorth, weight);
    mon.fillHisto("METpar", "afterWeight_PUPPI"+currentEvt.s_jetCat, PUPPIMETpar, weight);

    if(boson.Pt() < 55.) continue;
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after pt cut
    eventflowStep++;

    //LO-to-NLO k-factor for GJets_HT
    if(isMC_GJet_HT) weight *= max(1., 1.716910-0.001221*boson.Pt()); //We use the same weights than the ones used by JME-17-001. However when the weight becomes lower than one (at 587.15 GeV) we keep the weight =1. This looks like the weights we found when comparing our LO samples to our NLO samples.

    //Phi(Z,MET)
    if(currentEvt.deltaPhi_MET_Boson<0.5) continue;
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after delta phi (Z, met)
    eventflowStep++;

    //No Extra Lepton
    if(selElectrons.size()+extraElectrons.size()+selMuons.size()+extraMuons.size()>0) continue;
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after no extra leptons
    eventflowStep++;

      //b veto
      bool passBTag = true;
      for(int i =0 ; i < btags.size() ; i++){
        if (btags[i] > 0.5426) passBTag = false;
      }
      if(!passBTag) continue;

      //mon.fillHisto("eventflow","tot"+tagsR[c],eventflowStep,weight); //after b-tag veto
      for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight);
      eventflowStep++;

      //Phi(jet,MET)
      bool passDeltaPhiJetMET = true;
      for(int i = 0 ; i < selJets.size() ; i++){
        if (fabs(utils::deltaPhi(selJets[i], METVector))<0.5) passDeltaPhiJetMET = false;
      }
      if(!passDeltaPhiJetMET) continue;

      //mon.fillHisto("eventflow","tot"+tagsR[c],eventflowStep,weight); //after delta phi (jet, met)
      for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight);
      eventflowStep++;

    mon.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweighting", weight);

    // -- Histograms used to compute weights for the Instr. MET estimation : NVtx part --
    if(METVector.Pt()<125){
      mon.fillHisto("reco-vtx_MET125",    "InstrMET_reweighting"+currentEvt.s_jetCat+currentEvt.s_lepCat, EvtVtxCnt, weight, true);
      mon.fillHisto("reco-vtx_MET125",    "InstrMET_reweighting"+currentEvt.s_lepCat, EvtVtxCnt, weight, true); //for all jet cats
    }
      
    //Apply NVtx reweighting if file exist!
    //Starting from here, plots won't be "gamma" anymore but "eeR" or "mumuR". R for Reweighted.   
    double weightBeforeLoop = weight;
    double MTBeforeLoop = currentEvt.MT;
    double MBeforeLoop = currentEvt.M_Boson;
    TLorentzVector bosonBeforeLoop = boson;

    for(unsigned int c = 0; c < tagsR_size; c++){
      weight = weightBeforeLoop;
      currentEvt.MT = MTBeforeLoop;
      currentEvt.M_Boson = MBeforeLoop;
      boson = bosonBeforeLoop;


      if(c > 0){ //c=0 corresponds to no reweighting
        std::map<double, std::pair<double,double> >::iterator itlow;
        itlow = NVtxWeight_map[tagsR[c]].upper_bound(EvtVtxCnt); //look at which bin in the map currentEvt.rho corresponds
        if(itlow == NVtxWeight_map[tagsR[c]].begin()) throw std::out_of_range("You are trying to access your NVtx reweighting map outside of bin boundaries)");
        itlow--;

        weight *= itlow->second.first; //don't apply for first element of the map which is the normal one without reweighting.
      }

      mon.fillHisto("eventflow","tot"+tagsR[c],eventflowStep,weight); //after ee or mumu reweighting
      eventflowStep++;

      mon.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweightingAfter"+tagsR[c], weight);

      // -- Histograms used to compute weights for the Instr. MET estimation : Pt part --
      if(METVector.Pt()<125){
        mon.fillHisto("pT_Boson_MET125", "InstrMET_reweightingAfter"+tagsR[c]+currentEvt.s_jetCat, boson.Pt(), weight, true); // all jet cats
        mon.fillHisto("pT_Boson_MET125", "InstrMET_reweightingAfter"+tagsR[c], boson.Pt(), weight, true); // all jet cats
      }

      //Apply pt weight on top of NVtxWeight... so if i>0:
      if(c > 0 && weight_Pt_exist){
        std::map<double, std::pair<double,double> >::iterator itlow;
        itlow = PtWeight_map[tagsR[c]+currentEvt.s_jetCat].upper_bound(currentEvt.pT_Boson); //look at which bin in the map currentEvt.pT corresponds
        if(itlow == PtWeight_map[tagsR[c]+currentEvt.s_jetCat].begin()) throw std::out_of_range("You are trying to access your Pt reweighting map outside of bin boundaries)");
        itlow--;
        weight *= itlow->second.first; //don't apply for first element of the map which is the normal one without reweighting.
      }

      mon.fillHisto("eventflow","tot"+tagsR[c],eventflowStep,weight); //after Pt reweighting
      eventflowStep++;

      mon.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweightingAfter"+tagsR[c]+"AfterPtR", weight);

      //Apply mass on the photon:
      if(c > 0 && weight_Mass_exist){
        utils::giveMassToPhoton(boson, LineshapeMassWeight_map[tagsR[c]]);
        currentEvt.MT = sqrt(pow(sqrt(pow(boson.Pt(),2)+pow(boson.M(),2))+sqrt(pow(METVector.Pt(),2)+pow(91.1876,2)),2)-pow((boson+METVector).Pt(),2));
        currentEvt.M_Boson = boson.M();
      }

      mon.fillInstrMETControlRegionHisto(currentEvt, "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[c], weight);
      mon.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweightingAfter"+tagsR[c]+"AfterPtR_andMassivePhoton", weight);

      mon.fillAnalysisHistos(currentEvt, "beforeMETcut_After"+tagsR[c], weight);

      //MET>80
      if(METVector.Pt()<80) continue;
      mon.fillHisto("eventflow","tot"+tagsR[c],eventflowStep,weight); //after MET > 80
      eventflowStep++;

      //MET>125
      if(METVector.Pt()<125) continue;
      mon.fillHisto("eventflow","tot"+tagsR[c],eventflowStep,weight); //after MET > 125
      eventflowStep++;

      //###############################################################
      //##################     END OF SELECTION      ##################
      //###############################################################
      mon.fillAnalysisHistos(currentEvt, "final"+tagsR[c], weight, divideFinalHistoByBinWidth);
      //THE RESULT:
      mon.fillHisto("mT_final"+currentEvt.s_jetCat, currentEvt.s_lepCat, currentEvt.MT, weight, divideFinalHistoByBinWidth);

    }

  }

  //###############################################################
  //##################        END OF LOOP        ##################
  //###############################################################

  TFile* outFile=TFile::Open(outputFile_,"recreate");
  mon.Write();
  outFile->Close();

}
