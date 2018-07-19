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
#include "../Common/EwkCorrections.h"
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
  bool isMC_NLO_ZGTo2NuG_inclusive = (isMC_ && fileName.Contains("-ZGTo2NuG_") && !fileName.Contains("PtG-130"));
  bool isMC_NLO_ZGTo2NuG_Pt130 = (isMC_ && fileName.Contains("-ZGTo2NuG_PtG-130_"));

  //###############################################################
  //################## DECLARATION OF HISTOGRAMS ##################
  //###############################################################

  SmartSelectionMonitor_hzz mon;
  mon.declareHistos();

  Long64_t nbytes = 0, nb = 0;
  cout << "nb of entries in the input file =" << nentries << endl;

  cout << "fileName is " << fileName << endl;

  bool applyElectroweakCorrections = (fileName.Contains("ZZTo2L2Nu") || fileName.Contains("WZTo3LNu")  && !(fileName.Contains("GluGlu") || fileName.Contains("VBF")));
  if(applyElectroweakCorrections) cout << "Will apply electroweak corrections." << endl;
  else cout << "Will NOT apply electroweak corrections." << endl;

  // Table for electroweak corrections.
  vector<vector<float>> ewkTable;
  if(applyElectroweakCorrections) ewkTable = EwkCorrections::readFile_and_loadEwkTable(fileName);

  enum {ee, mumu, ll, lepCat_size};
  enum {eq0jets,geq1jets,vbf, jetCat_size};
  std::vector<string> v_jetCat = {"_eq0jets","_geq1jets","_vbf"};
  std::vector<TString> tagsR = {"_ee", "_mumu", "_ll"};
  unsigned int tagsR_size =  tagsR.size();  
  // ***--- Instr. MET building ---*** \\
  //Compute once weights for Instr. MET reweighting if needed
  string base_path = std::string(getenv("CMSSW_BASE")) + "/src/shears/HZZ2l2nu/";
  string weightFileType = "InstrMET";
  bool weight_NVtx_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root");
  bool weight_Pt_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root");
  bool weight_Mass_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_lineshape_mass.root");
  std::map<TString, std::map<double, std::pair<double, double> > > NVtxWeight_map, PtWeight_map;
  std::map<TString, TH1D*> LineshapeMassWeight_map;
  TH1 *h_mT = (TH1*) mon.getHisto("mT", "toGetBins", false); int h_mT_size = h_mT->GetNbinsX();
  TH1 *h_Vtx = (TH1*) mon.getHisto("reco-vtx", "toGetBins", false); int h_Vtx_size = h_Vtx->GetNbinsX();
  TH1 *h_pT = (TH1*) mon.getHisto("pT_Boson", "toGetBins", false); int h_pT_size = h_pT->GetNbinsX();
  std::vector<std::vector<std::vector<std::vector<std::vector<std::pair<double, double> > > > > > mT_InstrMET_map(lepCat_size, std::vector<std::vector<std::vector<std::vector<std::pair<double, double> > > > >(jetCat_size, std::vector<std::vector<std::vector<std::pair<double, double> > > >(h_mT_size+1, std::vector<std::vector<std::pair<double, double> > >(h_Vtx_size+1, std::vector<std::pair<double, double> >(h_pT_size+1/*, std::pair<double, double>*/)))));
  std::vector<std::vector<std::vector<std::vector<std::pair<double, double> > > > > photon_reweighting(lepCat_size, std::vector<std::vector<std::vector<std::pair<double, double> > > >(jetCat_size, std::vector<std::vector<std::pair<double, double> > >(h_Vtx_size+1, std::vector<std::pair<double, double> >(h_pT_size+1/*, std::pair<double, double>*/))));

  if(isPhotonDatadriven_ && (!weight_NVtx_exist || !weight_Pt_exist || !weight_Mass_exist) ) throw std::logic_error("You tried to run datadriven method without having weights for Instr.MET. This is bad :-) Please compute weights first!");
  if(isPhotonDatadriven_){
    utils::loadInstrMETWeights(weight_NVtx_exist, weight_Pt_exist, weight_Mass_exist, NVtxWeight_map, PtWeight_map, LineshapeMassWeight_map, weightFileType, base_path, v_jetCat);
    for(unsigned int lepCat = 0; lepCat < tagsR.size(); lepCat++){
      for(unsigned int jetCat = 0; jetCat < v_jetCat.size(); jetCat++){
        for(unsigned int Vtx = 1; Vtx <= h_Vtx_size; Vtx++){
          //1. #Vtx
          std::map<double, std::pair<double, double> >::iterator Vtx_low;
          Vtx_low = NVtxWeight_map[tagsR[lepCat]].upper_bound(h_Vtx->GetBinCenter(Vtx)); //look at which bin in the map this nVtx corresponds
          if(Vtx_low == NVtxWeight_map[tagsR[lepCat]].begin()) continue;
          Vtx_low--;
          for(unsigned int pT = 1; pT <= h_pT_size; pT++){
            //2. Pt
            std::map<double, std::pair<double, double> >::iterator pT_low;
            pT_low = PtWeight_map[tagsR[lepCat]+v_jetCat[jetCat]].upper_bound(h_pT->GetBinCenter(pT)); //look at which bin in the map this pT corresponds
            if(pT_low == PtWeight_map[tagsR[lepCat]+v_jetCat[jetCat]].begin()) continue;
            pT_low--;
            photon_reweighting[lepCat][jetCat][Vtx][pT].first = Vtx_low->second.first * pT_low->second.first;
            photon_reweighting[lepCat][jetCat][Vtx][pT].second = sqrt(Vtx_low->second.second*Vtx_low->second.second*pT_low->second.first*pT_low->second.first + pT_low->second.second*pT_low->second.second*Vtx_low->second.first*Vtx_low->second.first);
          }
        }
      }
    }
  }
  // ***--- End Instr. MET ---*** \\

  bool divideFinalHistoByBinWidth = false; //For final plots, we don't divide by the bin width to ease computations of the yields by eye.

  //###############################################################
  //##################     EVENT LOOP STARTS     ##################
  //###############################################################

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if ((jentry>maxEvents_)&&(maxEvents_>=0)) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    std::time_t currentTime = std::time(nullptr);
    if(jentry % 10000 ==0) cout << jentry << " of " << nentries << ". It is now " << std::asctime(std::localtime(&currentTime));

    evt currentEvt;

    double weight = 1.;
    double totEventWeight = 1.;
    //get the MC event weight if exists
    if (isMC_) {
      //get the MC event weight if exists
      weight *= (EvtWeights->size()>0 ? EvtWeights->at(0) : 1);
      if ((sumWeightInBonzai_>0)&&(sumWeightInBaobab_>0)) totEventWeight = weight*sumWeightInBaobab_/sumWeightInBonzai_;
      //get the PU weights
      float weightPU = pileUpWeight(EvtPuCntTruth); 
      weight *= weightPU;
    }
    else {
      totEventWeight = totalEventsInBaobab_/nentries;
    }

    mon.fillHisto("totEventInBaobab","tot",EvtPuCnt,totEventWeight);
    mon.fillHisto("eventflow","tot",0,weight);

    // Remove events with 0 vtx
    if(EvtVtxCnt == 0 ) continue;

    for(int i =0 ; i < MuPt->size() ; i++) mon.fillHisto("pT_mu","tot",MuPt->at(i),weight);
    for(int i =0 ; i < ElPt->size() ; i++) mon.fillHisto("pT_e","tot",ElPt->at(i),weight);
    mon.fillHisto("nb_mu","tot",MuPt->size(),weight);
    mon.fillHisto("nb_e","tot",ElPt->size(),weight);
    mon.fillHisto("pile-up","tot",EvtPuCnt,weight);
    mon.fillHisto("truth-pile-up","tot",EvtPuCntTruth,weight);
    mon.fillHisto("reco-vtx","tot",EvtVtxCnt,weight);


    //###############################################################
    //##################     OBJECT CORRECTIONS    ##################
    //###############################################################
    // muon momentum correction (Rochester)
    vector<float> *correctedMuPt = computeCorrectedMuPt(isMC_);

    // electroweak corrections
    map<string,pair<TLorentzVector,TLorentzVector>> genLevelLeptons;
    if(applyElectroweakCorrections) genLevelLeptons = EwkCorrections::reconstructGenLevelBosons(GLepBarePt, GLepBareEta, GLepBarePhi, GLepBareE, GLepBareId, GLepBareSt, GLepBareMomId);
    double ewkCorrections_error = 0.;
    double ewkCorrections_factor = 1.;
    if(applyElectroweakCorrections) ewkCorrections_factor = EwkCorrections::getEwkCorrections(fileName, genLevelLeptons, ewkTable, ewkCorrections_error, GPdfx1, GPdfx2, GPdfId1, GPdfId2);
    if(syst_=="ewk_up") weight *= (ewkCorrections_factor + ewkCorrections_error);
    else if (syst_=="ewk_down") weight *= (ewkCorrections_factor - ewkCorrections_error);
    else weight *= ewkCorrections_factor;

    // Theory uncertainties
    if(syst_ !="") weight *= utils::getTheoryUncertainties(EvtWeights, syst_);

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
    objectSelection::selectPhotons(selPhotons, PhotPt, PhotEta, PhotPhi, PhotId, PhotScEta, PhotHasPixelSeed, PhotSigmaIetaIeta, PhotSigmaIphiIphi, selMuons, selElectrons);
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
      int triggerWeight =0, triggerType;
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

    //Avoid double counting for W+jets
    //For some reasons we just have the inclusive sample for the Dilepton region while we have both HT and inclusive samples for the photon region. Hence this cleaning only applies to the photon region.
    if(isPhotonDatadriven_){
      if (isMC_Wlnu_inclusive || isMC_Wlnu_HT100){ //Avoid double counting and make our W#rightarrow l#nu exclusif of the dataset with a cut on HT...
        bool isHT100 = false;
        //Let's create our own gen HT variable
        double vHT =0;
        TLorentzVector genJet_uncleaned;
        for(size_t ig=0; ig<GJetAk04Pt->size(); ig++){
          genJet_uncleaned.SetPtEtaPhiE(GJetAk04Pt->at(ig), GJetAk04Eta->at(ig), GJetAk04Phi->at(ig), GJetAk04E->at(ig));
          double minDRmj(9999.); for(size_t ilepM=0; ilepM<selMuons.size();     ilepM++)  minDRmj = TMath::Min( minDRmj, utils::deltaR(genJet_uncleaned,selMuons[ilepM]) );
          double minDRej(9999.); for(size_t ilepE=0; ilepE<selElectrons.size(); ilepE++)  minDRej = TMath::Min( minDRej, utils::deltaR(genJet_uncleaned,selElectrons[ilepE]) );
          if(minDRmj<0.4 || minDRej<0.4) continue;
          vHT += GJetAk04Pt->at(ig);
        }
        if(vHT >100) isHT100 = true;
        if(isMC_Wlnu_inclusive) mon.fillHisto("custom_HT","forWlnu_inclusive",vHT,weight);
        if(isMC_Wlnu_HT100) mon.fillHisto("custom_HT","forWlnu_HT100",vHT,weight);
        if(isMC_Wlnu_inclusive && isHT100) continue; //reject event
        if(isMC_Wlnu_HT100 && !isHT100) continue; //reject event
      }

      //Avoid double counting for NLO ZvvG:
      if( isMC_NLO_ZGTo2NuG_inclusive && selPhotons[0].Pt() >= 130) continue;
      if( isMC_NLO_ZGTo2NuG_Pt130 && selPhotons[0].Pt() < 130) continue;

    }

    //Definition of the relevant analysis variables
    vector<TLorentzVectorWithIndex> selLeptons;
    if(isEE) selLeptons = selElectrons;
    if(isMuMu) selLeptons = selMuons;
    TLorentzVector boson = (isPhotonDatadriven_) ? selPhotons[0] : selLeptons[0] + selLeptons[1];
    TLorentzVector METVector; METVector.SetPtEtaPhiE(METPtType1XY->at(0),0.,METPhiType1XY->at(0),METPtType1XY->at(0));
    int jetCat = geq1jets;
    if(selJets.size()==0) jetCat = eq0jets;
    else if(utils::passVBFcuts(selJets, boson)) jetCat = vbf;
    currentEvt.s_jetCat = v_jetCat[jetCat];

    //Loop on lepton type. This is important also to apply Instr.MET if needed:
    double weightBeforeLoop = weight;
    TLorentzVector bosonBeforeLoop = boson;
    double photon_reweighting_tot = 1.;
    for(unsigned int c = 0; c < tagsR_size; c++){
      weight = weightBeforeLoop;
      boson = bosonBeforeLoop;
      photon_reweighting_tot = 1.;
      int lepCat = c;

      if(!isPhotonDatadriven_){
        if(tagsR[c] == "_ee" && !isEE) continue;
        else if(tagsR[c] == "_mumu" && !isMuMu) continue;
      }
      else{
        //Apply photon reweighting
        //1. #Vtx
        std::map<double, std::pair<double, double> >::iterator itlow;
        itlow = NVtxWeight_map[tagsR[c]].upper_bound(EvtVtxCnt); //look at which bin in the map currentEvt.nVtx corresponds
        if(itlow == NVtxWeight_map[tagsR[c]].begin()) throw std::out_of_range("You are trying to access your NVtx reweighting map outside of bin boundaries");
        itlow--;
        weight *= itlow->second.first; //(itlow->second.first = reweighting value; itlow->second.second = reweighting error)
        //2. Pt
        std::map<double, std::pair<double, double> >::iterator itlow2;
        itlow2 = PtWeight_map[tagsR[c]+v_jetCat[jetCat]].upper_bound(boson.Pt()); //look at which bin in the map currentEvt.pT corresponds
        if(itlow2 == PtWeight_map[tagsR[c]+v_jetCat[jetCat]].begin()) throw std::out_of_range("You are trying to access your Pt reweighting map outside of bin boundaries");
        itlow2--;
        weight *= itlow2->second.first; 
        //3. Mass lineshape
        utils::giveMassToPhoton(boson, LineshapeMassWeight_map[tagsR[c]]);

        photon_reweighting_tot = 1.*itlow->second.first*itlow2->second.first;
        if(photon_reweighting_tot == 0) continue;
      }
      //Jet category

      //Warning, starting from here ALL plots have to have the currentEvt.s_lepCat in their name, otherwise the reweighting will go crazy
      currentEvt.Fill_evt(v_jetCat[jetCat], tagsR[c], boson, METVector, selJets, EvtRunNum, EvtVtxCnt, EvtFastJetRho, METsig->at(0), selLeptons);

      mon.fillHisto("jetCategory","tot"+currentEvt.s_lepCat,jetCat,weight);
      mon.fillHisto("nJets","tot"+currentEvt.s_lepCat,currentEvt.nJets,weight);

      if(currentEvt.s_lepCat == "_ll") mon.fillAnalysisHistos(currentEvt, "tot", weight);

      if(fabs(boson.M()-91.1876)>15.) continue; //We keep this cut for the photon... because we've already given a mass to the photon :)
      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",2,weight);
      mon.fillAnalysisHistos(currentEvt, "inZpeak", weight);

      if(boson.Pt() < 55.) continue;
      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",3,weight);

      //Phi(Z,MET)
      if(currentEvt.deltaPhi_MET_Boson<0.5) continue;
      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",4,weight);

      if(extraElectrons.size()>0 || extraMuons.size()>0) continue;
      if(isPhotonDatadriven_ && (selMuons.size()>0 || selElectrons.size()>0) ) continue;
      if(isEE && selMuons.size()>0) continue;
      if(isMuMu && selElectrons.size()>0) continue;
      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",5,weight);

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

      mon.fillInstrMETControlRegionHisto(currentEvt, "InstrMET_reweighting", weight);
      mon.fillAnalysisHistos(currentEvt, "beforeMETcut", weight);
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
      mon.fillHisto("jetCategory","final"+currentEvt.s_lepCat,jetCat,weight, divideFinalHistoByBinWidth);
      mon.fillAnalysisHistos(currentEvt, "final", weight, divideFinalHistoByBinWidth);


      //Prepare the correct computation of the stat uncertainty for the mT plots with Instr.MET.
      if(isPhotonDatadriven_){
        int mT = min(h_mT_size, h_mT->FindBin(currentEvt.MT));
        int Vtx = min(h_Vtx_size, h_Vtx->FindBin(currentEvt.nVtx));
        int pT = min(h_pT_size , h_pT->FindBin(currentEvt.pT_Boson));
        mT_InstrMET_map[lepCat][jetCat][mT][Vtx][pT].first += 1.*weight/photon_reweighting_tot; //Fill with the weight before photon reweighting
        mT_InstrMET_map[lepCat][jetCat][mT][Vtx][pT].second += 1.*weight*weight/(photon_reweighting_tot*photon_reweighting_tot); //Errors due to weights without photon reweighting
      }
    }
  }

  //###############################################################
  //##################        END OF LOOP        ##################
  //###############################################################

  TFile* outFile=TFile::Open(outputFile_,"recreate");

  //Construct mT plots with correct stat uncertainty for the Instr.MET part.
  if(isPhotonDatadriven_){
    double content =0, error2 = 0, N_events = 0, N_error2 = 0, reweighting = 0, reweighting_error = 0;
    for(unsigned int lepCat = 0; lepCat < tagsR.size(); lepCat++){
      for(unsigned int jetCat = 0; jetCat < v_jetCat.size(); jetCat++){
        for(unsigned int mT = 1; mT <= h_mT_size; mT++){
          content = 0;
          error2 = 0;
          for(unsigned int Vtx = 1; Vtx <= h_Vtx_size; Vtx++){
            for(unsigned int pT = 1; pT <= h_pT_size; pT++){
              N_events = mT_InstrMET_map[lepCat][jetCat][mT][Vtx][pT].first; //sum of weights
              N_error2 = mT_InstrMET_map[lepCat][jetCat][mT][Vtx][pT].second; //sum of weights*weights
              reweighting = photon_reweighting[lepCat][jetCat][Vtx][pT].first;
              reweighting_error = photon_reweighting[lepCat][jetCat][Vtx][pT].second;
              content += N_events*reweighting;
              error2 += reweighting*reweighting*N_error2 + N_events*N_events*reweighting_error*reweighting_error;
            }
          }
          mon.setBinContentAndError("mT", "final"+v_jetCat[jetCat]+tagsR[lepCat], mT, content, sqrt(error2), divideFinalHistoByBinWidth); //erase the mT final plot
        }
      }
    }
  }

  mon.WriteForSysts(syst_,keepAllControlPlots_);
  outFile->Close();

}
