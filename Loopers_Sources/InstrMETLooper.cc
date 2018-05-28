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
  string base_path = std::string(getenv("CMSSW_BASE")) + "/src/shears/HZZ2l2nu/";
  bool doClosureTest = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/please_do_closure_test_when_running_InstrMETLooper");
  if(doClosureTest) std::cout << "/!\\/!\\ CLOSURE TEST ONGOING - not wanted? Then remove 'WeightsAndDatadriven/InstrMET/please_do_closure_test_when_running_InstrMETLooper' /!\\/!\\" << std::endl;
  std::string weightFileType = (doClosureTest) ? "closureTest" : "InstrMET";
  bool weight_NVtx_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root");
  bool weight_Pt_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root");
  bool weight_Mass_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_lineshape_mass.root");
  std::map<TString, std::map<double, double> > NVtxWeight_map, PtWeight_map;
  std::map<TString, TH1D*> LineshapeMassWeight_map;
  if(weight_NVtx_exist){
    std::cout << "NVtx weight file has been found! Some histo (called 'After_eeR' and 'After_mumuR') will have the NVtx reweighting applied :)" << std::endl;
    NVtxWeight_map["_eeR"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root", "WeightHisto__ee_AllBins");
    NVtxWeight_map["_mumuR"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root", "WeightHisto__mumu_AllBins");
    NVtxWeight_map["_llR"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root", "WeightHisto__ll_AllBins");
    if(weight_Pt_exist){
      std::cout << "Pt weight file has also been found! Some histo (called 'AfterPtR') will have both reweighting applied :)" << std::endl;
      PtWeight_map["_eeR"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto__ee_AllBins");
      PtWeight_map["_mumuR"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto__mumu_AllBins");
      PtWeight_map["_llR"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto__ll_AllBins");
      for(unsigned int i =0; i < v_jetCat.size(); i++){
        PtWeight_map["_eeR"+v_jetCat[i]] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto"+v_jetCat[i]+"_ee_AllBins");
        PtWeight_map["_mumuR"+v_jetCat[i]] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto"+v_jetCat[i]+"_mumu_AllBins");
        PtWeight_map["_llR"+v_jetCat[i]] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto"+v_jetCat[i]+"_ll_AllBins");
      }
    }
    if(weight_Mass_exist){
      std::cout << "Lineshape mass file has been found! Some histo (named 'andMassivePhoton') will have the lineshape applied :)" << std::endl;
      TFile *f_weight_lineshape = TFile::Open((TString) base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_lineshape_mass.root");
      LineshapeMassWeight_map["_eeR"] = (TH1D*) ((TH1D*) f_weight_lineshape->Get("WeightHisto__ee_AllBins"))->Clone();
      LineshapeMassWeight_map["_mumuR"] = (TH1D*) ((TH1D*) f_weight_lineshape->Get("WeightHisto__mumu_AllBins"))->Clone();
      LineshapeMassWeight_map["_llR"] = (TH1D*) ((TH1D*) f_weight_lineshape->Get("WeightHisto__ll_AllBins"))->Clone();
      //f_weight_lineshape->Close();
    }
  }
  std::vector<TString> tagsR;
  tagsR.push_back("_gamma"); //_gamma, i.e. no reweighting to ee or mumu
  if(weight_NVtx_exist){
    tagsR.push_back("_eeR");
    tagsR.push_back("_mumuR");
    tagsR.push_back("_llR");
  }
  unsigned int tagsR_size =  tagsR.size();

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
      weight = (EvtWeights->size()>0 ? EvtWeights->at(0) : 1);
      if ((sumWeightInBonzai_>0)&&(sumWeightInBaobab_>0)) totEventWeight = weight*sumWeightInBaobab_/sumWeightInBonzai_;
    }
    else {
      totEventWeight = totalEventsInBaobab_/nentries;
    }

    //Cleaning of low stats MC spikes that are gathering in some specific spot (in MET phi, pt, MET delta phi(MET, spike)...).
    bool isPathologicEvent=false;
    if(isMC_) isPathologicEvent = objectSelection::cleanPathologicEventsInPhotons(fileName, EvtRunNum, EvtLumiNum, EvtNum);
    if(isPathologicEvent) continue;

    // Remove events with 0 vtx
    if(EvtVtxCnt == 0 ) continue;

    mon.fillHisto("totEventInBaobab","tot",EvtPuCnt,totEventWeight);
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //output of bonzais
    eventflowStep++;

    //###############################################################
    //##################     OBJECT SELECTION      ##################
    //###############################################################

    vector<TLorentzVectorWithIndex> selElectrons; //Leptons passing final cuts
    vector<TLorentzVectorWithIndex> selMuons; //Muons passing final cuts
    vector<TLorentzVectorWithIndex> extraElectrons; //Additional electrons, used for veto
    vector<TLorentzVectorWithIndex> extraMuons; //Additional muons, used for veto
    vector<TLorentzVectorWithIndex> selPhotons; //Photons
    vector<TLorentzVectorWithIndex> selJets; //Jets passing Id and cleaning, with |eta|<4.7 and pT>30GeV. Used for jet categorization and deltaPhi cut.
    vector<double> btags; //B-Tag discriminant, recorded for selJets with |eta|<2.5. Used for b-tag veto.

    objectSelection::selectElectrons(selElectrons, extraElectrons, ElPt, ElEta, ElPhi, ElE, ElId, ElEtaSc);
    objectSelection::selectMuons(selMuons, extraMuons, MuPt, MuEta, MuPhi, MuE, MuId, MuIdTight, MuIdSoft, MuPfIso);
    objectSelection::selectPhotons(selPhotons, PhotPt, PhotEta, PhotPhi, PhotId, PhotScEta, PhotHasPixelSeed, PhotSigmaIetaIeta, selMuons, selElectrons);
    objectSelection::selectJets(selJets, btags, JetAk04Pt, JetAk04Eta, JetAk04Phi, JetAk04E, JetAk04Id, JetAk04NeutralEmFrac, JetAk04NeutralHadAndHfFrac, JetAk04NeutMult, JetAk04BDiscCisvV2, selMuons, selElectrons, selPhotons);

    //Ask for a prompt photon
    if(selPhotons.size() != 1) continue;

    //Check trigger and find prescale
    int triggerWeight =0;
    int triggerType;
    if(isMC_) triggerType = trigger::MC_Photon;
    else triggerType = trigger::SinglePhoton;

    triggerWeight = trigger::passTrigger(triggerType, TrigHltDiMu, TrigHltMu, TrigHltDiEl, TrigHltEl, TrigHltElMu, TrigHltPhot, TrigHltDiMu_prescale, TrigHltMu_prescale, TrigHltDiEl_prescale, TrigHltEl_prescale, TrigHltElMu_prescale, TrigHltPhot_prescale, selPhotons[0].Pt());
    if(triggerWeight==0) continue; //trigger not found

    if(MAXIMAL_AMOUNT_OF_HISTOS) mon.fillHisto("pT_Z","noPrescale",selPhotons[0].Pt(),weight);
    weight *= triggerWeight;
    if(MAXIMAL_AMOUNT_OF_HISTOS) mon.fillHisto("pT_Z","withPrescale",selPhotons[0].Pt(),weight);
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after prescale
    eventflowStep++;

    //photon efficiencies
    PhotonEfficiencySF phoEff;
    if(isMC_) weight *= phoEff.getPhotonEfficiency(selPhotons[0].Pt(), PhotScEta->at(selPhotons[0].GetIndex()), "tight",utils::CutVersion::Moriond17Cut ).first; 
    if(MAXIMAL_AMOUNT_OF_HISTOS) mon.fillHisto("pT_Z","withPrescale_and_phoEff",selPhotons[0].Pt(),weight);

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
    TLorentzVector PUPPIMETVector; PUPPIMETVector.SetPtEtaPhiE(METPtType1XY->at(2),0.,METPhiType1XY->at(2),METPtType1XY->at(2));
    double transverseMass_PUPPI = sqrt(pow(sqrt(pow(boson.Pt(),2)+pow(boson.M(),2))+sqrt(pow(PUPPIMETVector.Pt(),2)+pow(91.1876,2)),2)-pow((boson+PUPPIMETVector).Pt(),2));
    currentEvt.MT = sqrt(pow(sqrt(pow(boson.Pt(),2)+pow(boson.M(),2))+sqrt(pow(METVector.Pt(),2)+pow(91.1876,2)),2)-pow((boson+METVector).Pt(),2));
    currentEvt.M = boson.M();
    currentEvt.pT = boson.Pt();
    currentEvt.MET = METVector.Pt();
    currentEvt.eta = boson.Eta();
    currentEvt.s_lepCat = tagsR[0];
    double minDeltaPhiJetMET = 4.;
    for(int i = 0 ; i < selJets.size() ; i++){
      if (fabs(utils::deltaPhi(selJets[i], METVector)) < minDeltaPhiJetMET) minDeltaPhiJetMET = fabs(utils::deltaPhi(selJets[i], METVector));
    }
    // compute the parallele and the orthogonal MET
    double METorth = 0, PUPPIMETorth = 0;
    double METpar = 0, PUPPIMETpar = 0;
    TVector3 bosonDir = boson.Vect();
    if (bosonDir.Perp()>0){
      TVector3 bosonDir2D(bosonDir.x()/bosonDir.Perp(), bosonDir.y()/bosonDir.Perp(), 0);
      TVector3 bosonPerp2D(-bosonDir2D.y(),bosonDir2D.x(),0);
      TVector3 METVector3D = METVector.Vect();
      TVector3 PUPPIMETVector3D = PUPPIMETVector.Vect();
      METpar = - METVector3D.Dot(bosonDir2D);
      METorth = METVector3D.Dot(bosonPerp2D);
      PUPPIMETpar = - PUPPIMETVector3D.Dot(bosonDir2D);
      PUPPIMETorth = PUPPIMETVector3D.Dot(bosonPerp2D);
    }
    currentEvt.METphi = METVector.Phi();
    currentEvt.deltaPhi_MET_jet = minDeltaPhiJetMET;
    currentEvt.deltaPhi_MET_phot = fabs(utils::deltaPhi(boson, METVector));
    currentEvt.METoPT = METVector.Pt()/(1.*boson.Pt());
    currentEvt.METpar = METpar;
    currentEvt.METperp = METorth;
    currentEvt.METsig = METsig->at(0);
    currentEvt.nVtx = EvtVtxCnt;
    currentEvt.HoE = PhotHoE->at(selPhotons[0].GetIndex());
    currentEvt.sigmaIEtaIEta = PhotSigmaIetaIeta->at(selPhotons[0].GetIndex());
    currentEvt.chIsoRhoCorr =  utils::photon_rhoCorrectedIso(PhotPfIsoChHad->at(selPhotons[0].GetIndex()), EvtFastJetRho, PhotScEta->at(selPhotons[0].GetIndex()), "chIso");
    currentEvt.neuIsoRhoCorr = utils::photon_rhoCorrectedIso(PhotPfIsoNeutralHad->at(selPhotons[0].GetIndex()), EvtFastJetRho, PhotScEta->at(selPhotons[0].GetIndex()), "nhIso");
    currentEvt.phoIsoRhoCorr = utils::photon_rhoCorrectedIso(PhotPfIsoPhot->at(selPhotons[0].GetIndex()), EvtFastJetRho, PhotScEta->at(selPhotons[0].GetIndex()), "gIso");
    currentEvt.R9 = PhotR9->at(selPhotons[0].GetIndex());
    currentEvt.rho = EvtFastJetRho;
    currentEvt.jet0_pT = ((selJets.size() > 0) ? selJets[0].Pt() : 0);
    currentEvt.jet1_pT = ((selJets.size() > 1) ? selJets[1].Pt() : 0);
    currentEvt.jet2_pT = ((selJets.size() > 2) ? selJets[2].Pt() : 0);
    currentEvt.jet3_pT = ((selJets.size() > 3) ? selJets[3].Pt() : 0);
    double allRecoJets_HT = 0.;
    double allSelJets_HT = 0.;
    for( int i =0; i< JetAk04Pt->size(); i++) allRecoJets_HT += JetAk04Pt->at(i);
    for (unsigned int i =0; i < selJets.size(); i++) allSelJets_HT += selJets[i].Pt();
    double HT_PFTightJetsAbove5GeV  = 0;
    for(unsigned int i =0 ; i<JetAk04Pt->size() ; i++){
      bool passSelPt = false, passEta = false, passTightEta = false, passId = false, passLeptonCleaning = false, passPhotonCleaning = false;
      TLorentzVectorWithIndex currentJet = TLorentzVectorWithIndex::PtEtaPhiEIndex(JetAk04Pt->at(i),JetAk04Eta->at(i),JetAk04Phi->at(i),JetAk04E->at(i), i);
      passSelPt = (currentJet.Pt() >=5);
      double eta = fabs(currentJet.Eta());
      passEta = (eta <= 4.7);
      passTightEta = (eta <= 2.5);
      if(eta<2.7){
        passId = (JetAk04Id->at(i) >= 1); //This case is simple, it corresponds exactly to what we apply for now
        }
      float nef = JetAk04NeutralEmFrac->at(i);
      float nhf = JetAk04NeutralHadAndHfFrac->at(i);
      float nnp = JetAk04NeutMult->at(i);
      if(eta<3.0 && eta >=2.7){
        passId = (nef > 0.01 && nhf < 0.98 && nnp > 2);
        //passId = (JetAk04Id->at(i) >= 2); //Simpler criterium, but not equivalent to what is mentionned in the AN. Was applied before having access to nnp.
        }
      if(eta>=3.0){
        passId = (nef<0.90 && nnp > 10);
        //passId = (JetAk04Id->at(i) >= 3); //Simpler criterium, but not equivalent to what is mentionned in the AN. Was applied before having access to nnp.
        }
      double minDRlj(9999.); for(unsigned int ilep=0; ilep<selMuons.size(); ilep++) minDRlj = TMath::Min( minDRlj, utils::deltaR(currentJet,selMuons[ilep]) );
      for(unsigned int ilep=0; ilep<selElectrons.size(); ilep++) minDRlj = TMath::Min( minDRlj, utils::deltaR(currentJet,selElectrons[ilep]) );
      passLeptonCleaning = (minDRlj>=0.4);
      double minDRgj(9999.); for(unsigned int ilep=0; ilep<selPhotons.size(); ilep++) minDRgj = TMath::Min( minDRgj, utils::deltaR(currentJet,selPhotons[ilep]) );
      passPhotonCleaning = (minDRgj>=0.4);
      if(passSelPt && passEta && passId && passLeptonCleaning && passPhotonCleaning) HT_PFTightJetsAbove5GeV += currentJet.Pt();
    }
    currentEvt.HT_selJets = allSelJets_HT;
    currentEvt.HT_PFTightJetsAbove5GeV = HT_PFTightJetsAbove5GeV;
    currentEvt.nb_selPhoton = selPhotons.size();
    double deltaPhiZMet = fabs(utils::deltaPhi(boson, METVector));



    //Jet category
    enum {eq0jets,geq1jets,vbf};
    int jetCat = geq1jets;
    if(selJets.size()==0) jetCat = eq0jets;
    if(utils::passVBFcuts(selJets, boson)) jetCat = vbf;
    currentEvt.s_jetCat = v_jetCat[jetCat];
    currentEvt.nJets = selJets.size();

    mon.fillHisto("jetCategory","afterWeight",jetCat,weight);
    mon.fillAnalysisHistos_InstrMET(currentEvt, "afterWeight", weight);

    //std::cout<<"Event info: " << EvtRunNum<<":"<<EvtLumiNum<<":"<<EvtNum << "; boson pt = "<<boson.Pt()<<"; weight = "<<weight<<"; triggerPrescale = "<<triggerWeight<<"; met = "<<currentEvt.MET<<"; mt = "<<currentEvt.MT<<"; njets = "<<currentEvt.nJets<<"; vtx = "<<EvtVtxCnt<<"; rho = "<<EvtFastJetRho<<"; puWeight = "<<weightPU<<std::endl;

    //More MET variables
    // -- PUPPI MET
    mon.fillHisto("MET", "afterWeight_PUPPI"+currentEvt.s_jetCat, PUPPIMETVector.Pt(), weight, true);
    mon.fillHisto("MET_phi", "afterWeight_PUPPI"+currentEvt.s_jetCat, PUPPIMETVector.Phi(), weight, true);
    mon.fillHisto("DeltaPhi_MET_Phot", "afterWeight_PUPPI"+currentEvt.s_jetCat, fabs(utils::deltaPhi(boson, PUPPIMETVector)), weight);
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

    //Addition of control plots
    //R9
    mon.fillHisto("R9", "afterWeight_photonStudy"+currentEvt.s_jetCat, PhotR9->at(selPhotons[0].GetIndex()), weight);
    //Sigma ieta ieta
    mon.fillHisto("SigmaIetaIeta", "afterWeight_photonStudy"+currentEvt.s_jetCat, PhotSigmaIetaIeta->at(selPhotons[0].GetIndex()), weight);
    //Isolation
    mon.fillHisto("PfIsoChHad", "afterWeight_photonStudy"+currentEvt.s_jetCat, PhotPfIsoChHad->at(selPhotons[0].GetIndex()), weight, true);
    mon.fillHisto("PfIsoChHad_showQCD", "afterWeight_photonStudy"+currentEvt.s_jetCat, PhotPfIsoChHad->at(selPhotons[0].GetIndex()), weight, true);
    mon.fillHisto("PfIsoNeutralHad", "afterWeight_photonStudy"+currentEvt.s_jetCat, PhotPfIsoNeutralHad->at(selPhotons[0].GetIndex()), weight);
    mon.fillHisto("PfIsoPhot", "afterWeight_photonStudy"+currentEvt.s_jetCat, PhotPfIsoPhot->at(selPhotons[0].GetIndex()), weight);
    mon.fillHisto("IsoEcal", "afterWeight_photonStudy"+currentEvt.s_jetCat, PhotIsoEcal->at(selPhotons[0].GetIndex()), weight);
    mon.fillHisto("IsoHcal", "afterWeight_photonStudy"+currentEvt.s_jetCat, PhotIsoHcal->at(selPhotons[0].GetIndex()), weight);
    mon.fillHisto("IsoTk", "afterWeight_photonStudy"+currentEvt.s_jetCat, PhotIsoTk->at(selPhotons[0].GetIndex()), weight);
    mon.fillHisto("PfIsoPuChHad", "afterWeight_photonStudy"+currentEvt.s_jetCat, PhotPfIsoPuChHad->at(selPhotons[0].GetIndex()), weight);
    mon.fillHisto("PfIsoEcalClus", "afterWeight_photonStudy"+currentEvt.s_jetCat, PhotPfIsoEcalClus->at(selPhotons[0].GetIndex()), weight);
    mon.fillHisto("PfIsoHcalClus", "afterWeight_photonStudy"+currentEvt.s_jetCat, PhotPfIsoHcalClus->at(selPhotons[0].GetIndex()), weight);
    //HoverE
    mon.fillHisto("HoE", "afterWeight_photonStudy"+currentEvt.s_jetCat, PhotHoE->at(selPhotons[0].GetIndex()), weight);


    //Reco HT
    mon.fillHisto("custom_HT","afterWeight_allRecoJets_HT"+currentEvt.s_jetCat,allRecoJets_HT,weight);
    mon.fillHisto("custom_HT","afterWeight_allSelJets_HT"+currentEvt.s_jetCat,allSelJets_HT,weight);

    
    if(boson.Pt() < 55.) continue;
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after pt cut
    eventflowStep++;

    //LO-to-NLO k-factor for GJets_HT
    if(isMC_GJet_HT) weight *= max(1., 1.716910-0.001221*boson.Pt()); //We use the same weights than the ones used by JME-17-001. However when the weight becomes lower than one (at 587.15 GeV) we keep the weight =1. This looks like the weights we found when comparing our LO samples to our NLO samples.





    //Phi(Z,MET)
    if(deltaPhiZMet<0.5) continue;
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after delta phi (Z, met)
    eventflowStep++;

    //No Extra Lepton
    if(selElectrons.size()+extraElectrons.size()+selMuons.size()+extraMuons.size()>0) continue;
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after no extra leptons
    eventflowStep++;

    mon.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweighting", weight);

    // -- Histograms used to compute weights for the Instr. MET estimation : Nvtx part --
    mon.fillHisto("pT_Z",        "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, boson.Pt(), weight);
    mon.fillHisto("reco-vtx",    "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, EvtVtxCnt,  weight);
    mon.fillHisto("reco-vtx",    "InstrMET_reweighting__"+currentEvt.s_lepCat, EvtVtxCnt,  weight); //for all jet cats
    mon.fillHisto("zpt_vs_nvtx", "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, boson.Pt(), EvtVtxCnt, weight);
    mon.fillHisto("zpt_vs_nvtx_ee"+currentEvt.s_jetCat, "InstrMET_reweighting_"+currentEvt.s_lepCat, boson.Pt(), EvtVtxCnt, weight);
    mon.fillHisto("zpt_vs_nvtx_mumu"+currentEvt.s_jetCat, "InstrMET_reweighting_"+currentEvt.s_lepCat, boson.Pt(), EvtVtxCnt, weight);

    //if(currentEvt.s_jetCat == "_eq0jets"){
    //  mon.fillHisto("zpt_vs_nvtx_ee_eq0jets", "InstrMET_reweighting_"+currentEvt.s_lepCat, boson.Pt(), EvtVtxCnt, weight);
    //  mon.fillHisto("zpt_vs_nvtx_mumu_eq0jets", "InstrMET_reweighting_"+currentEvt.s_lepCat, boson.Pt(), EvtVtxCnt, weight);
    //}
    //else if(currentEvt.s_jetCat == "_geq1jets"){
    //  mon.fillHisto("zpt_vs_nvtx_ee_geq1jets", "InstrMET_reweighting_"+currentEvt.s_lepCat, boson.Pt(), EvtVtxCnt, weight);
    //  mon.fillHisto("zpt_vs_nvtx_mumu_geq1jets", "InstrMET_reweighting_"+currentEvt.s_lepCat, boson.Pt(), EvtVtxCnt, weight);
    //}
    //else if(currentEvt.s_jetCat == "_vbf"){
    //   mon.fillHisto("zpt_vs_nvtx_ee_vbf", "InstrMET_reweighting_"+currentEvt.s_lepCat, boson.Pt(), EvtVtxCnt, weight);
    //   mon.fillHisto("zpt_vs_nvtx_mumu_vbf", "InstrMET_reweighting_"+currentEvt.s_lepCat, boson.Pt(), EvtVtxCnt, weight);
    //}

    //Apply Nvtx reweighting if file exist!
    //Starting from here, plots won't be "gamma" anymore but "eeR" or "mumuR". R for Reweighted.   
    double weightBeforeLoop = weight;
    double MTBeforeLoop = currentEvt.MT;
    double MBeforeLoop = currentEvt.M;
    for(unsigned int i = 0; i < tagsR_size; i++){
      weight = weightBeforeLoop;
      currentEvt.MT = MTBeforeLoop;
      currentEvt.M = MBeforeLoop;


      if(i > 0){ //i=0 corresponds to no reweighting
        std::map<double,double>::iterator itlow;
        itlow = NVtxWeight_map[tagsR[i]].upper_bound(currentEvt.nVtx); //look at which bin in the map currentEvt.nVtx corresponds
        if(itlow == NVtxWeight_map[tagsR[i]].begin()) throw std::out_of_range("You are trying to access your NVtx reweighting map outside of bin boundaries)");
        itlow--;
        
        weight *= itlow->second ; //don't apply for first element of the map which is the normal one without reweighting.
      }

      mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after ee or mumu reweighting
      eventflowStep++;
 
      mon.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweightingAfter"+tagsR[i], weight);

      // -- Histograms used to compute weights for the Instr. MET estimation : Pt part --
      mon.fillHisto("pT_Z",        "InstrMET_reweightingAfter"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, boson.Pt(), weight);
      mon.fillHisto("pT_Z",        "InstrMET_reweightingAfter"+tagsR[i]+"__"+currentEvt.s_lepCat, boson.Pt(), weight); // all jets
      mon.fillHisto("reco-vtx",    "InstrMET_reweightingAfter"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, EvtVtxCnt,  weight);
      mon.fillHisto("reco-vtx",    "InstrMET_reweightingAfter"+tagsR[i]+"__"+currentEvt.s_lepCat, EvtVtxCnt,  weight); //for all jet cats

      //Apply pt weight on top of NVtxWeight... so if i>0:
      if(i > 0 && weight_Pt_exist){
        std::map<double,double>::iterator itlow;
        itlow = PtWeight_map[tagsR[i]+currentEvt.s_jetCat].upper_bound(currentEvt.pT); //look at which bin in the map currentEvt.pT corresponds
        if(itlow == PtWeight_map[tagsR[i]+currentEvt.s_jetCat].begin()) throw std::out_of_range("You are trying to access your Pt reweighting map outside of bin boundaries)");
        itlow--;

        weight *= itlow->second ; //don't apply for first element of the map which is the normal one without reweighting.
      }
      
      mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after Pt reweighting
      eventflowStep++;
 
      mon.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweightingAfter"+tagsR[i]+"AfterPtR", weight);

      //Apply mass on the photon:
      if(i > 0 && weight_Mass_exist){
        utils::giveMassToPhoton(boson, LineshapeMassWeight_map[tagsR[i]] );
        currentEvt.MT = sqrt(pow(sqrt(pow(boson.Pt(),2)+pow(boson.M(),2))+sqrt(pow(METVector.Pt(),2)+pow(91.1876,2)),2)-pow((boson+METVector).Pt(),2));
        currentEvt.M = boson.M();
      }
     
      //Plots for closure test
      mon.fillHisto("pT_Z",        "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, boson.Pt(), weight);
      mon.fillHisto("pT_Z",        "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, boson.Pt(), weight); // all jets
      mon.fillHisto("reco-vtx",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, EvtVtxCnt,  weight);
      mon.fillHisto("reco-vtx",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, EvtVtxCnt,  weight); //for all jet cats
      mon.fillHisto("M_Z",        "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, boson.M(), weight);
      mon.fillHisto("M_Z",        "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, boson.M(), weight); // all jets
      mon.fillHisto("MET",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.MET,  weight);
      mon.fillHisto("MET",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.MET,  weight); //for all jet cats
      mon.fillHisto("mT",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.MT,  weight);
      mon.fillHisto("mT",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.MT,  weight); //for all jet cats
      mon.fillHisto("DeltaPhi_MET_Phot",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.deltaPhi_MET_phot,  weight);
      mon.fillHisto("DeltaPhi_MET_Phot",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.deltaPhi_MET_phot, weight); //for all jet cats
      mon.fillHisto("DeltaPhi_MET_Jet",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.deltaPhi_MET_jet,  weight);
      mon.fillHisto("DeltaPhi_MET_Jet",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.deltaPhi_MET_jet, weight); //for all jet cats
      mon.fillHisto("METoverPt_zoom",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.METoPT,  weight);
      mon.fillHisto("METoverPt_zoom",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.METoPT, weight); //for all jet cats
      mon.fillHisto("eta_Z",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.eta,  weight);
      mon.fillHisto("eta_Z",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.eta, weight); //for all jet cats
      mon.fillHisto("pT_jet0",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.jet0_pT,  weight);
      mon.fillHisto("pT_jet0",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.jet0_pT, weight); //for all jet cats
      mon.fillHisto("nJets",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.nJets,  weight);
      mon.fillHisto("nJets",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.nJets, weight); //for all jet cats
      mon.fillHisto("custom_HT",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, allSelJets_HT,  weight);
      mon.fillHisto("custom_HT",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, allSelJets_HT, weight); //for all jet cats
      if(allSelJets_HT > 300){
        mon.fillHisto("MET_HT300",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.MET,  weight);
        mon.fillHisto("MET_HT300",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.MET, weight); //for all jet cats
      }
      if(boson.Pt() < 300){
        mon.fillHisto("MET_Pt0-300",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.MET,  weight);
        mon.fillHisto("MET_Pt0-300",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.MET, weight); //for all jet cats
      }
      else if(boson.Pt() < 400){
        mon.fillHisto("MET_Pt300-400",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.MET,  weight);
        mon.fillHisto("MET_Pt300-400",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.MET, weight); //for all jet cats
      }
      else if(boson.Pt() < 600){
        mon.fillHisto("MET_Pt400-600",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.MET,  weight);
        mon.fillHisto("MET_Pt400-600",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.MET, weight); //for all jet cats
      }
      else{
        mon.fillHisto("MET_Pt600-Inf",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.MET,  weight);
        mon.fillHisto("MET_Pt600-Inf",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.MET, weight); //for all jet cats
      }

      //TProfile for closure test
      mon.fillProfile("METvsBosonPt",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, boson.Pt(), currentEvt.MET,  weight);
      mon.fillProfile("METvsBosonPt",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, boson.Pt(), currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("METvsMT",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.MT, currentEvt.MET,  weight);
      mon.fillProfile("METvsMT",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.MT, currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("METvsDPhiMETBos",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.deltaPhi_MET_phot, currentEvt.MET, weight);
      mon.fillProfile("METvsDPhiMETBos",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.deltaPhi_MET_phot, currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("METvsDPhiMETJet",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.deltaPhi_MET_jet, currentEvt.MET, weight);
      mon.fillProfile("METvsDPhiMETJet",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.deltaPhi_MET_jet, currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("METvsJetPt",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.jet0_pT, currentEvt.MET, weight);
      mon.fillProfile("METvsJetPt",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.jet0_pT, currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("METvsNJets",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.nJets, currentEvt.MET, weight);
      mon.fillProfile("METvsNJets",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.nJets, currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("METvsBosonEta",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.eta, currentEvt.MET, weight);
      mon.fillProfile("METvsBosonEta",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, currentEvt.eta, currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("METvsHT",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, allSelJets_HT, currentEvt.MET, weight);
      mon.fillProfile("METvsHT",    "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[i]+"__"+currentEvt.s_lepCat, allSelJets_HT, currentEvt.MET, weight); //for all jet cats

      mon.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweightingAfter"+tagsR[i]+"AfterPtR_andMassivePhoton", weight);

      //FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME 
      //FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME
      // This is where I should add the first reweighting. Then I could add a second set of histo for the second reweighting. So in one run I can see:
      // 1) Important plots without reweighting
      // 2) Important plots with NVtx reweighting
      // 3) Important plots with full reweighting
      // L'ideal serait aussi d'ajouter ca dans l'eventflow, pour voir comment ca affecte le bazar
      //FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME 
      //FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME

      //b veto
      bool passBTag = true;
      for(int i =0 ; i < btags.size() ; i++){
        if (btags[i] > 0.5426) passBTag = false;
      }
      if(!passBTag) continue;

      mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after b-tag veto
      eventflowStep++;

      //Phi(jet,MET)
      bool passDeltaPhiJetMET = true;
      for(int i = 0 ; i < selJets.size() ; i++){
        if (fabs(utils::deltaPhi(selJets[i], METVector))<0.5) passDeltaPhiJetMET = false;
      }
      if(!passDeltaPhiJetMET) continue;

      mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after delta phi (jet, met)
      eventflowStep++;

      mon.fillAnalysisHistos_InstrMET(currentEvt, "beforeMETcut_After"+tagsR[i], weight);

      //MET>80
      if(METVector.Pt()<80) continue;
      mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after MET > 80
      eventflowStep++;

      //MET>125
      if(METVector.Pt()<125) continue;
      mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after MET > 125
      eventflowStep++;

      //###############################################################
      //##################     END OF SELECTION      ##################
      //###############################################################
      mon.fillAnalysisHistos_InstrMET(currentEvt, "final", weight);

    }

  }

  //###############################################################
  //##################        END OF LOOP        ##################
  //###############################################################

  TFile* outFile=TFile::Open(outputFile_,"recreate");
  mon.Write();
  outFile->Close();

}
