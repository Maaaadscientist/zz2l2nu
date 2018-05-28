#define HZZ2l2nuLooper_cxx
#include "../Loopers_Includes/LooperMain.h"
#include "../Loopers_Includes/SmartSelectionMonitor.h"
#include "../Loopers_Includes/SmartSelectionMonitor_hzz.h"
#include "../Common/Utils.h"
#include "../Common/ObjectSelection.h"
#include "../Common/TLorentzVectorWithIndex.h"
#include "../Common/Trigger.h"
#include "../Common/PhotonEfficiencySF.h"
#include "../Common/LeptonsEfficiencySF.h"
#include <ctime>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

void LooperMain::Loop()
{
  if (fChain == 0) return;

  //Get file info
  Long64_t nentries = fChain->GetEntries(); //Warning: GetEntries has to be called before any other work on fChain
  TString fileName = fChain->GetCurrentFile()->GetName();
  bool isMC_Wlnu_inclusive = (isMC_ && fileName.Contains("-WJetsToLNu_") && !fileName.Contains("HT"));
  bool isMC_Wlnu_HT100 = (isMC_ && fileName.Contains("-WJetsToLNu_HT-") );

  //###############################################################
  //################## DECLARATION OF HISTOGRAMS ##################
  //###############################################################

  SmartSelectionMonitor_hzz mon;
  mon.declareHistos();

  Long64_t nbytes = 0, nb = 0;
  cout << "nb of entries in the input file =" << nentries << endl;

  std::vector<string> v_jetCat = {"_eq0jets","_geq1jets","_vbf"};
  // ***--- Instr. MET building ---*** \\
  //Compute once weights for Instr. MET reweighting if needed
  string base_path = std::string(getenv("CMSSW_BASE")) + "/src/shears/HZZ2l2nu/";
  std::string weightFileType = "InstrMET";
  bool weight_NVtx_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root");
  bool weight_Pt_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root");
  bool weight_Mass_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_lineshape_mass.root");
  std::map<TString, std::map<double, double> > NVtxWeight_map, PtWeight_map;
  std::map<TString, TH1D*> LineshapeMassWeight_map;
  std::vector<TString> tagsR;
  tagsR.push_back("_ee");
  tagsR.push_back("_mumu");
  tagsR.push_back("_ll");
  if(isPhotonDatadriven_ && (!weight_NVtx_exist || !weight_Pt_exist || !weight_Mass_exist) ) throw std::logic_error("You tried to run datadriven method without having weights for Instr.MET. This is bad :-) Please compute weights first!");
  if(isPhotonDatadriven_){
    if(weight_NVtx_exist){
      std::cout << "NVtx weight file has been found! Some histo (called 'After_eeR' and 'After_mumuR') will have the NVtx reweighting applied :)" << std::endl;
      NVtxWeight_map["_ee"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root", "WeightHisto__ee_AllBins");
      NVtxWeight_map["_mumu"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root", "WeightHisto__mumu_AllBins");
      NVtxWeight_map["_ll"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root", "WeightHisto__ll_AllBins");
      if(weight_Pt_exist){
        std::cout << "Pt weight file has also been found! Some histo (called 'AfterPtR') will have both reweighting applied :)" << std::endl;
        PtWeight_map["_ee"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto__ee_AllBins");
        PtWeight_map["_mumu"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto__mumu_AllBins");
        PtWeight_map["_ll"] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto__ll_AllBins");
        for(unsigned int i =0; i < v_jetCat.size(); i++){
          PtWeight_map["_ee"+v_jetCat[i]] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto"+v_jetCat[i]+"_ee_AllBins");
          PtWeight_map["_mumu"+v_jetCat[i]] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto"+v_jetCat[i]+"_mumu_AllBins");
          PtWeight_map["_ll"+v_jetCat[i]] = utils::TH1toMap(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root", "WeightHisto"+v_jetCat[i]+"_ll_AllBins");
        }
      }
      if(weight_Mass_exist){
        std::cout << "Lineshape mass file has been found! Some histo (named 'andMassivePhoton') will have the lineshape applied :)" << std::endl;
        gRandom = new TRandom3(0); //New seed for the mass generation. If one wants to always have the same mass and mT distribution, please remove this line
        TFile *f_weight_lineshape = TFile::Open((TString) base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_lineshape_mass.root");
        LineshapeMassWeight_map["_ee"] = (TH1D*) ((TH1D*) f_weight_lineshape->Get("WeightHisto__ee_AllBins"))->Clone();
        LineshapeMassWeight_map["_mumu"] = (TH1D*) ((TH1D*) f_weight_lineshape->Get("WeightHisto__mumu_AllBins"))->Clone();
        LineshapeMassWeight_map["_ll"] = (TH1D*) ((TH1D*) f_weight_lineshape->Get("WeightHisto__ll_AllBins"))->Clone();
      }
    }
  }
  unsigned int tagsR_size =  tagsR.size();  
  // ***--- End Instr. MET ---*** \\



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
      //get the MC event weight if exists
      weight = (EvtWeights->size()>0 ? EvtWeights->at(0) : 1);
      if ((sumWeightInBonzai_>0)&&(sumWeightInBaobab_>0)) totEventWeight = weight*sumWeightInBaobab_/sumWeightInBonzai_;
      //get the PU weights
      float weightPU = pileUpWeight(EvtPuCntTruth); 
      weight = weight*weightPU;
    }
    else {
      totEventWeight = totalEventsInBaobab_/nentries;
    }

    // Remove events with 0 vtx
    if(EvtVtxCnt == 0 ) continue;

    mon.fillHisto("totEventInBaobab","tot",EvtPuCnt,totEventWeight);
    mon.fillHisto("eventflow","tot",0,weight);

    for(int i =0 ; i < MuPt->size() ; i++) mon.fillHisto("pT_mu","tot",MuPt->at(i),weight);
    for(int i =0 ; i < ElPt->size() ; i++) mon.fillHisto("pT_e","tot",ElPt->at(i),weight);
    mon.fillHisto("nb_mu","tot",MuPt->size(),weight);
    mon.fillHisto("nb_e","tot",ElPt->size(),weight);
    mon.fillHisto("pile-up","tot",EvtPuCnt,weight);
    mon.fillHisto("truth-pile-up","tot",EvtPuCntTruth,weight);
    mon.fillHisto("reco-vtx","tot",EvtVtxCnt,weight);


    //###############################################################
    //##################     OBJECTS CORRECTIONS   ##################
    //###############################################################
    // muon momentum correction (Rochester)
    vector<float> *correctedMuPt = computeCorrectedMuPt(isMC_);


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
    objectSelection::selectMuons(selMuons, extraMuons, correctedMuPt, MuEta, MuPhi, MuE, MuId, MuIdTight, MuIdSoft, MuPfIso);
    objectSelection::selectPhotons(selPhotons, PhotPt, PhotEta, PhotPhi, PhotId, PhotScEta, PhotHasPixelSeed, PhotSigmaIetaIeta, selMuons, selElectrons);
    objectSelection::selectJets(selJets, btags, JetAk04Pt, JetAk04Eta, JetAk04Phi, JetAk04E, JetAk04Id, JetAk04NeutralEmFrac, JetAk04NeutralHadAndHfFrac, JetAk04NeutMult, JetAk04BDiscCisvV2, selMuons, selElectrons, selPhotons);

    //Discriminate ee and mumu
    bool isEE = (selElectrons.size()==2 && !isPhotonDatadriven_); //2 good electrons
    bool isMuMu = (selMuons.size()==2 && !isPhotonDatadriven_); //2 good muons
    bool isGamma = (selPhotons.size() == 1 && isPhotonDatadriven_); //1 good photon

    mon.fillHisto("nb_mu","sel",selMuons.size(),weight);
    mon.fillHisto("nb_e","sel",selElectrons.size(),weight);
    mon.fillHisto("nb_mu","extra",extraMuons.size(),weight);
    mon.fillHisto("nb_e","extra",extraElectrons.size(),weight);

    //###############################################################
    //##################       ANALYSIS CUTS       ##################
    //###############################################################


    if(!isEE && !isMuMu && !isGamma) continue; //not a good lepton pair or photon (if datadriven)
    mon.fillHisto("eventflow","tot",1,weight);

    if(isEE) currentEvt.s_lepCat = "_ee";
    else if(isMuMu) currentEvt.s_lepCat = "_mumu";

    //compute and apply the efficiency SFs
    if (isMC_){
      if(!isPhotonDatadriven_){ //for leptons
        float weightLeptonsSF= (isEE ? trigAndIDsfs::diElectronEventSFs(utils::CutVersion::CutSet::Moriond17Cut, selElectrons[0].Pt(), ElEtaSc->at(selElectrons[0].GetIndex()), selElectrons[1].Pt(), ElEtaSc->at(selElectrons[1].GetIndex())) : trigAndIDsfs::diMuonEventSFs( utils::CutVersion::CutSet::Moriond17Cut, MuPt->at(selMuons[0].GetIndex()), selMuons[0].Eta(), MuPt->at(selMuons[1].GetIndex()), selMuons[1].Eta()));
        weight*=weightLeptonsSF;
      }
      else{ //for photons
        PhotonEfficiencySF phoEff;
        weight *= phoEff.getPhotonEfficiency(selPhotons[0].Pt(), PhotScEta->at(selPhotons[0].GetIndex()), "tight",utils::CutVersion::Moriond17Cut ).first;
      }
    }

    //trigger weights for photon data
    if(isPhotonDatadriven_){
      int triggerWeight =0;
      int triggerType;
      if(isMC_) triggerType = trigger::MC_Photon;
      else triggerType = trigger::SinglePhoton;

      triggerWeight = trigger::passTrigger(triggerType, TrigHltDiMu, TrigHltMu, TrigHltDiEl, TrigHltEl, TrigHltElMu, TrigHltPhot, TrigHltDiMu_prescale, TrigHltMu_prescale, TrigHltDiEl_prescale, TrigHltEl_prescale, TrigHltElMu_prescale, TrigHltPhot_prescale, selPhotons[0].Pt());
      if(triggerWeight==0) continue; //trigger not found
      weight *= triggerWeight;
    }

    //MET filters
    std::vector<std::pair<int, int> > listMETFilter; //after the passMetFilter function, it contains the bin number of the cut in .first and if it passed 1 or not 0 the METfilter
    bool passMetFilter = utils::passMetFilter(TrigMET, listMETFilter, isMC_);
    //now fill the metFilter eventflow
    mon.fillHisto("metFilters","tot",26,weight); //the all bin, i.e. the last one
    for(unsigned int i =0; i < listMETFilter.size(); i++){
      if(listMETFilter[i].second ==1) mon.fillHisto("metFilters","tot",listMETFilter[i].first,weight);
    }
    if (!passMetFilter) continue;

    //Avoid double couting for W+jets
    //For some reasons we just have the inclusive sample for the Dilepton region while we have both HT and inclusive samples for the photon region. Hence this cleaning only applies to the photon region.
    if(isPhotonDatadriven_){
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
        if(isMC_Wlnu_inclusive) mon.fillHisto("custom_HT","forWlnu_inclusive",vHT,weight);
        if(isMC_Wlnu_HT100) mon.fillHisto("custom_HT","forWlnu_HT100",vHT,weight);
        if(isMC_Wlnu_inclusive && isHT100) continue; //reject event
        if(isMC_Wlnu_HT100 && !isHT100) continue; //reject event

      }
    }

    //Definition of the relevant analysis variables
    vector<TLorentzVectorWithIndex> selLeptons;
    if(isEE) selLeptons = selElectrons;
    if(isMuMu) selLeptons = selMuons;
    TLorentzVector boson = (isPhotonDatadriven_) ? selPhotons[0] : selLeptons[0] + selLeptons[1];

    //Loop on lepton type. This is important also to apply Instr.MET if needed:
    TLorentzVector METVector; METVector.SetPtEtaPhiE(METPtType1XY->at(0),0.,METPhiType1XY->at(0),METPtType1XY->at(0));
    currentEvt.transverseMass = sqrt(pow(sqrt(pow(boson.Pt(),2)+pow(boson.M(),2))+sqrt(pow(METVector.Pt(),2)+pow(91.1876,2)),2)-pow((boson+METVector).Pt(),2));
    double weightBeforeLoop = weight;
    double MTBeforeLoop = currentEvt.transverseMass;
    double MZBeforeLoop = boson.M();
    TLorentzVector bosonBeforeLoop = boson;
    for(unsigned int c = 0; c < tagsR_size; c++){
      weight = weightBeforeLoop;
      currentEvt.transverseMass = MTBeforeLoop;
      currentEvt.MZ = MZBeforeLoop;
      boson = bosonBeforeLoop;
      
      if(!isPhotonDatadriven_){
        if(tagsR[c] == "_ee" && !isEE) continue;
        else if(tagsR[c] == "_mumu" && !isMuMu) continue;
      }
      else{
        //Apply photon reweighting
        //1. #Vtx
        std::map<double,double>::iterator itlow;
        itlow = NVtxWeight_map[tagsR[c]].upper_bound(EvtVtxCnt); //look at which bin in the map currentEvt.nVtx corresponds
        if(itlow == NVtxWeight_map[tagsR[c]].begin()) throw std::out_of_range("You are trying to access your NVtx reweighting map outside of bin boundaries");
        itlow--;
        weight *= itlow->second ;

        //2. Pt
        std::map<double,double>::iterator itlow2;
        itlow2 = PtWeight_map[tagsR[c]+currentEvt.s_jetCat].upper_bound(boson.Pt()); //look at which bin in the map currentEvt.pT corresponds
        if(itlow2 == PtWeight_map[tagsR[c]+currentEvt.s_jetCat].begin()) throw std::out_of_range("You are trying to access your Pt reweighting map outside of bin boundaries");
        itlow2--;
        weight *= itlow2->second ; 

        //3. Mass lineshape
        utils::giveMassToPhoton(boson, LineshapeMassWeight_map[tagsR[c]]);
      }
      currentEvt.s_lepCat = tagsR[c]; //Warning, starting from here ALL plots have to have the currentEvt.s_lepCat in their name, otherwise the reweighting will go crazy

      currentEvt.transverseMass = sqrt(pow(sqrt(pow(boson.Pt(),2)+pow(boson.M(),2))+sqrt(pow(METVector.Pt(),2)+pow(91.1876,2)),2)-pow((boson+METVector).Pt(),2));
      currentEvt.MZ = boson.M();
      currentEvt.pTZ = boson.Pt();
      currentEvt.MET = METVector.Pt();
      currentEvt.METphi = METVector.Phi();
      currentEvt.runNumber= EvtRunNum; 
      if(!isPhotonDatadriven_){
        if (selLeptons[0].Pt() > selLeptons[1].Pt()) {
          currentEvt.lep1pT = selLeptons[0].Pt();
          currentEvt.lep1eta = selLeptons[0].Eta();
          currentEvt.lep2pT = selLeptons[1].Pt();
          currentEvt.lep2eta = selLeptons[1].Eta();
        }
        else {
          currentEvt.lep1pT = selLeptons[1].Pt();
          currentEvt.lep1eta = selLeptons[1].Eta();
          currentEvt.lep2pT = selLeptons[0].Pt();
          currentEvt.lep2eta = selLeptons[0].Eta();
        }
      }
      currentEvt.nVtx = EvtVtxCnt;


      //Jet category
      enum {eq0jets,geq1jets,vbf};
      int jetCat = geq1jets;
      if(selJets.size()==0) jetCat = eq0jets;
      else if(utils::passVBFcuts(selJets, boson)) jetCat = vbf;
      currentEvt.s_jetCat = v_jetCat[jetCat];
      mon.fillHisto("jetCategory","tot"+currentEvt.s_lepCat,jetCat,weight);
      currentEvt.nJets = selJets.size();
      mon.fillHisto("nJets","tot"+currentEvt.s_lepCat,currentEvt.nJets,weight);

      if(currentEvt.s_lepCat == "_ll") mon.fillAnalysisHistos(currentEvt, "tot", weight);

      if(fabs(boson.M()-91.1876)>15.) continue; //We keep this cut for the photon... because we've already given a mass to the photon :)
      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",2,weight);
      mon.fillAnalysisHistos(currentEvt, "inZpeak"+currentEvt.s_lepCat, weight);

      if(boson.Pt() < 55.) continue;
      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",3,weight);

      //Phi(Z,MET)
      double deltaPhiZMet = fabs(utils::deltaPhi(boson, METVector));
      if(deltaPhiZMet<0.5) continue;
      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",4,weight);

      if(extraElectrons.size()>0 || extraMuons.size()>0) continue;
      if(isPhotonDatadriven_ && (selMuons.size()>0 || selElectrons.size()>0) ) continue;
      if(isEE && selMuons.size()>0) continue;
      if(isMuMu && selElectrons.size()>0) continue;
      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",5,weight);

      // -- Histograms used to compute weights for the Instr. MET estimation --
      mon.fillHisto("reco-vtx","InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat,EvtVtxCnt,weight);
      mon.fillHisto("reco-vtx","InstrMET_reweighting__"+currentEvt.s_lepCat,EvtVtxCnt,weight); //for all jet cats
      mon.fillHisto("pT_Z","InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat,currentEvt.pTZ,weight);
      mon.fillHisto("pT_Z","InstrMET_reweighting__"+currentEvt.s_lepCat,currentEvt.pTZ,weight);
      mon.fillHisto("M_Z","InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat,currentEvt.MZ,weight);
      mon.fillHisto("M_Z","InstrMET_reweighting__"+currentEvt.s_lepCat,currentEvt.MZ,weight);
      mon.fillHisto("MET","InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat,currentEvt.MET,weight);
      mon.fillHisto("MET","InstrMET_reweighting__"+currentEvt.s_lepCat,currentEvt.MET,weight);
      mon.fillHisto("mT","InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat,currentEvt.transverseMass,weight);
      mon.fillHisto("mT","InstrMET_reweighting__"+currentEvt.s_lepCat,currentEvt.transverseMass,weight);
      double minDeltaPhiJetMET = 4.;
      for(int i = 0 ; i < selJets.size() ; i++){
        if (fabs(utils::deltaPhi(selJets[i], METVector)) < minDeltaPhiJetMET) minDeltaPhiJetMET = fabs(utils::deltaPhi(selJets[i], METVector));
      }
      mon.fillHisto("DeltaPhi_MET_Phot",    "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, fabs(utils::deltaPhi(boson, METVector)),  weight);
      mon.fillHisto("DeltaPhi_MET_Phot",    "InstrMET_reweighting__"+currentEvt.s_lepCat, fabs(utils::deltaPhi(boson, METVector)), weight); //for all jet cats
      mon.fillHisto("DeltaPhi_MET_Jet",    "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, minDeltaPhiJetMET,  weight);
      mon.fillHisto("DeltaPhi_MET_Jet",    "InstrMET_reweighting__"+currentEvt.s_lepCat, minDeltaPhiJetMET, weight); //for all jet cats
      mon.fillHisto("METoverPt_zoom",    "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, METVector.Pt()/(1.*boson.Pt()),  weight);
      mon.fillHisto("METoverPt_zoom",    "InstrMET_reweighting__"+currentEvt.s_lepCat, METVector.Pt()/(1.*boson.Pt()), weight); //for all jet cats
      mon.fillHisto("eta_Z",    "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, boson.Eta(),  weight);
      mon.fillHisto("eta_Z",    "InstrMET_reweighting__"+currentEvt.s_lepCat, boson.Eta(), weight); //for all jet cats
      mon.fillHisto("pT_jet0",    "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, ((selJets.size() > 0) ? selJets[0].Pt() : 0),  weight);
      mon.fillHisto("pT_jet0",    "InstrMET_reweighting__"+currentEvt.s_lepCat, ((selJets.size() > 0) ? selJets[0].Pt() : 0), weight); //for all jet cats
      mon.fillHisto("nJets",    "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.nJets,  weight);
      mon.fillHisto("nJets",    "InstrMET_reweighting__"+currentEvt.s_lepCat, currentEvt.nJets, weight); //for all jet cats
      double allSelJets_HT =0.;
      for (unsigned int i =0; i < selJets.size(); i++) allSelJets_HT += selJets[i].Pt();
      mon.fillHisto("selJetsHT",    "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, allSelJets_HT, weight);
      mon.fillHisto("selJetsHT",    "InstrMET_reweighting__"+currentEvt.s_lepCat, allSelJets_HT, weight); //for all jet cats
      if(allSelJets_HT > 300){
        mon.fillHisto("MET_HT300",  "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.MET, weight);
        mon.fillHisto("MET_HT300",  "InstrMET_reweighting__"+currentEvt.s_lepCat, currentEvt.MET, weight); //for all jet cats
      }
      if(boson.Pt() < 300){
        mon.fillHisto("MET_Pt0-300",      "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.MET,  weight);
        mon.fillHisto("MET_Pt0-300",      "InstrMET_reweighting__"+currentEvt.s_lepCat, currentEvt.MET, weight); //for all jet cats
      }
      else if(boson.Pt() < 400){
        mon.fillHisto("MET_Pt300-400",    "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.MET,  weight);
        mon.fillHisto("MET_Pt300-400",    "InstrMET_reweighting__"+currentEvt.s_lepCat, currentEvt.MET, weight); //for all jet cats
      }
      else if(boson.Pt() < 600){
        mon.fillHisto("MET_Pt400-600",    "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.MET,  weight);
        mon.fillHisto("MET_Pt400-600",    "InstrMET_reweighting__"+currentEvt.s_lepCat, currentEvt.MET, weight); //for all jet cats
      }
      else{
        mon.fillHisto("MET_Pt600-Inf",    "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.MET,  weight);
        mon.fillHisto("MET_Pt600-Inf",    "InstrMET_reweighting__"+currentEvt.s_lepCat, currentEvt.MET, weight); //for all jet cats
      }

      //TProfile for closure test
      mon.fillProfile("METvsBosonPt",     "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, boson.Pt(), currentEvt.MET, weight);
      mon.fillProfile("METvsBosonPt",     "InstrMET_reweighting__"+currentEvt.s_lepCat, boson.Pt(), currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("METvsMT",          "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.transverseMass, currentEvt.MET, weight);
      mon.fillProfile("METvsMT",          "InstrMET_reweighting__"+currentEvt.s_lepCat, currentEvt.transverseMass, currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("METvsDPhiMETBos",  "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, fabs(utils::deltaPhi(boson, METVector)), currentEvt.MET, weight);
      mon.fillProfile("METvsDPhiMETBos",  "InstrMET_reweighting__"+currentEvt.s_lepCat, fabs(utils::deltaPhi(boson, METVector)), currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("METvsDPhiMETJet",  "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, minDeltaPhiJetMET, currentEvt.MET, weight);
      mon.fillProfile("METvsDPhiMETJet",  "InstrMET_reweighting__"+currentEvt.s_lepCat, minDeltaPhiJetMET, currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("METvsJetPt",       "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, ((selJets.size() > 0) ? selJets[0].Pt() : 0), currentEvt.MET, weight);
      mon.fillProfile("METvsJetPt",       "InstrMET_reweighting__"+currentEvt.s_lepCat, ((selJets.size() > 0) ? selJets[0].Pt() : 0), currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("METvsNJets",       "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, currentEvt.nJets, currentEvt.MET, weight);
      mon.fillProfile("METvsNJets",       "InstrMET_reweighting__"+currentEvt.s_lepCat, currentEvt.nJets, currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("METvsBosonEta",    "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, boson.Eta(), currentEvt.MET, weight);
      mon.fillProfile("METvsBosonEta",    "InstrMET_reweighting__"+currentEvt.s_lepCat, boson.Eta(), currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("METvsHT",          "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, allSelJets_HT, currentEvt.MET, weight);
      mon.fillProfile("METvsHT",          "InstrMET_reweighting__"+currentEvt.s_lepCat, allSelJets_HT, currentEvt.MET, weight); //for all jet cats
      mon.fillProfile("HTvsBosonEta",     "InstrMET_reweighting"+currentEvt.s_jetCat+"_"+currentEvt.s_lepCat, boson.Eta(), allSelJets_HT, weight);
      mon.fillProfile("HTvsBosonEta",     "InstrMET_reweighting__"+currentEvt.s_lepCat, boson.Eta(), allSelJets_HT, weight); //for all jet cats
     

      //b veto
      bool passBTag = true;
      for(int i =0 ; i < btags.size() ; i++){
        if (btags[i] > 0.5426) passBTag = false;
      }
      if(!passBTag) continue;

      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",6,weight);

      //Phi(jet,MET)
      bool passDeltaPhiJetMET = true;
      for(int i = 0 ; i < selJets.size() ; i++){
        if (fabs(utils::deltaPhi(selJets[i], METVector))<0.5) passDeltaPhiJetMET = false;
      }
      if(!passDeltaPhiJetMET) continue;

      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",7,weight);

      mon.fillAnalysisHistos(currentEvt, "beforeMETcut"+currentEvt.s_lepCat, weight);
      mon.fillHisto("reco-vtx","beforeMETcut"+currentEvt.s_lepCat,EvtVtxCnt,weight);
      mon.fillHisto("jetCategory","beforeMETcut"+currentEvt.s_lepCat,jetCat,weight);


      //MET>80
      if(METVector.Pt()<80) continue;
      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",8,weight);

      //MET>125
      if(METVector.Pt()<125) continue;
      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",9,weight);

      //###############################################################
      //##################     END OF SELECTION      ##################
      //###############################################################
      mon.fillHisto("reco-vtx","final"+currentEvt.s_lepCat,EvtVtxCnt,weight);
      mon.fillHisto("jetCategory","final"+currentEvt.s_lepCat,jetCat,weight);
      mon.fillAnalysisHistos(currentEvt, "final"+currentEvt.s_lepCat, weight);

    }
  }

  //###############################################################
  //##################        END OF LOOP        ##################
  //###############################################################

  TFile* outFile=TFile::Open(outputFile_,"recreate");
  mon.Write();
  outFile->Close();

}

