#define HZZ2l2nuLooper_cxx
#include "../Loopers_Includes/LooperMain.h"
#include "../Loopers_Includes/SmartSelectionMonitor.h"
#include "../Loopers_Includes/SmartSelectionMonitor_hzz.h"
#include "../Loopers_Includes/BTagger.h"
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
#include <algorithm>


bool LooperMain::passTrigger(int triggerType){
  bool passEventSelection = false;
  int triggerWeight=0;
  triggerWeight = trigger::passTrigger(triggerType, TrigHltDiMu, TrigHltMu, TrigHltDiEl, TrigHltEl, TrigHltElMu, TrigHltPhot, TrigHltDiMu_prescale, TrigHltMu_prescale, TrigHltDiEl_prescale, TrigHltEl_prescale, TrigHltElMu_prescale, TrigHltPhot_prescale);
  if(triggerWeight > 0) passEventSelection  = true;
  return passEventSelection;
}

void LooperMain::Loop_NRB()
{
  if (fChain == 0) return;

  //Get file info
  Long64_t nentries = fChain->GetEntries(); //Warning: GetEntries has to be called before any other work on fChain
  TString fileName = fChain->GetCurrentFile()->GetName();
  bool isMC_Wlnu_inclusive = (isMC_ && fileName.Contains("-WJetsToLNu_") && !fileName.Contains("HT"));
  bool isMC_Wlnu_HT100 = (isMC_ && fileName.Contains("-WJetsToLNu_HT-") );

  // For Baobab analysis
  int  triggerType;
  if      (!isMC_ && fileName.Contains("-DoubleMuon-"))     triggerType = trigger::DoubleMu;
  else if (!isMC_ && fileName.Contains("-DoubleEG-"))       triggerType = trigger::DoubleE ;
  else if (!isMC_ && fileName.Contains("-SingleMuon-"))     triggerType = trigger::SingleMu;
  else if (!isMC_ && fileName.Contains("-SingleElectron-")) triggerType = trigger::SingleE;
  else if (!isMC_ && fileName.Contains("-MuonEG-"))         triggerType = trigger::EMu ;
  else if ( isMC_)                                          triggerType = trigger::MC_DiLepton;
  //###############################################################
  //################## DECLARATION OF HISTOGRAMS ##################
  //###############################################################

  SmartSelectionMonitor_hzz mon;
  mon.declareHistos_NRB();

  Long64_t nbytes = 0, nb = 0;
  cout << "nb of entries in the input file =" << nentries << endl;

  cout << "fileName is " << fileName << endl;

  bool applyElectroweakCorrections =fileName.Contains("-ZZTo2L2Nu") || fileName.Contains("-WZTo3LNu") ;
  if (fileName.Contains("GluGluH")||fileName.Contains("VBF")) applyElectroweakCorrections = false;
  if(applyElectroweakCorrections) cout << "Will apply electroweak corrections." << endl;
  else cout << "Will NOT apply electroweak corrections." << endl;

  // Table for electroweak corrections.
  vector<vector<float>> ewkTable;
  if(applyElectroweakCorrections) ewkTable = EwkCorrections::readFile_and_loadEwkTable(fileName);
  // Table for btagging efficiencies
  utils::tables btagEffTable;
  if(isMC_) btagEffTable = btagger::loadEffTables();
  // CSV file for btagging SF
  BTagCalibrationReader _btag_calibration_reader = btagger::loadCalibrationReader();
  std::vector<string> v_jetCat = {"_eq0jets","_geq1jets","_vbf"};
  std::vector<string> tagsR = {"_ee", "_mumu", "_ll","_emu"};
  unsigned int tagsR_size =  tagsR.size();  

  enum {ee, mumu, ll, emu,lepCat_size};
  enum {eq0jets, geq1jets, vbf, jetCat_size};
 
  // ***--- Instr. MET building ---*** \\
  //Compute once weights for Instr. MET reweighting if needed
  string base_path = std::string(getenv("CMSSW_BASE")) + "/src/shears/HZZ2l2nu/";
  string weightFileType = "InstrMET";
  bool weight_NVtx_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root");
  bool weight_Pt_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root");
  bool weight_Mass_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_lineshape_mass.root");
  std::map<TString, std::map<double, double> > NVtxWeight_map, PtWeight_map;
  std::map<TString, TH1D*> LineshapeMassWeight_map;
 
  // ***--- End Instr. MET ---*** 
  std::vector<double> optim_Cuts1_met;
  optim_Cuts1_met.push_back(0); //add a bin in the shapes with a MET cut of 0
  for(double met=50;met<140;met+=5) {  optim_Cuts1_met    .push_back(met);  }
  TH2F *h_2D=(TH2F *) mon.addHistogram( new TH2F ("mt_shapes_NRBctrl",";cut Index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );

  h_2D->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
  h_2D->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
  h_2D->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
  h_2D->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
  h_2D->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
  h_2D->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");
  
  //Definition of the final histos (and in particular of the mT binning
  bool divideFinalHistoByBinWidth = false; //For final plots, we don't divide by the bin width to ease computations of the yields by eye.
  std::vector<TH1*> h_mT(jetCat_size); std::vector<int> h_mT_size(jetCat_size);
  h_mT[eq0jets] = (TH1*) mon.getHisto("mT_final_eq0jets", "ee", divideFinalHistoByBinWidth); h_mT_size[eq0jets] = h_mT[eq0jets]->GetNbinsX();
  h_mT[geq1jets] = (TH1*) mon.getHisto("mT_final_geq1jets", "ee", divideFinalHistoByBinWidth); h_mT_size[geq1jets] = h_mT[geq1jets]->GetNbinsX();
  h_mT[vbf] = (TH1*) mon.getHisto("mT_final_vbf", "ee", divideFinalHistoByBinWidth); h_mT_size[vbf] = h_mT[vbf]->GetNbinsX();
  mon.getHisto("mT_final_eq0jets", "mumu", divideFinalHistoByBinWidth); //The .substr(1) removes the annoying _ in the tagsR definition.
  mon.getHisto("mT_final_geq1jets", "mumu", divideFinalHistoByBinWidth);
  mon.getHisto("mT_final_vbf", "mumu", divideFinalHistoByBinWidth);
  int h_mT_maxSize = std::max({h_mT_size[eq0jets], h_mT_size[geq1jets], h_mT_size[vbf]}); 

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
    // Remove events with 0 vtx
    if(EvtVtxCnt == 0 ) continue;

    mon.fillHisto("totEventInBaobab","tot",EvtPuCnt,totEventWeight);
    if(runOnBaobabs_){
      if (!LooperMain::passTrigger(triggerType)) continue;
    }

    mon.fillHisto("eventflow","tot",0,weight);



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



    //###############################################################
    //##################     OBJECT SELECTION      ##################
    //###############################################################

    vector<TLorentzVectorWithIndex> selElectrons; //Leptons passing final cuts
    vector<TLorentzVectorWithIndex> selMuons; //Muons passing final cuts
    vector<TLorentzVectorWithIndex> extraElectrons; //Additional electrons, used for veto
    vector<TLorentzVectorWithIndex> extraMuons; //Additional muons, used for veto
    vector<TLorentzVectorWithIndex> selPhotons; //Photons
    vector<TLorentzVectorWithIndex> selJets; //Jets passing Id and cleaning, with |eta|<4.7 and pT>30GeV. Used for jet categorization and deltaPhi cut.
    vector<TLorentzVectorWithIndex> selCentralJets; //Same as the previous one, but with tracker acceptance (|eta| <= 2.5). Used to compute btag efficiency and weights. 
    vector<double> btags; //B-Tag discriminant, recorded for selCentralJets. Used for b-tag veto, efficiency and weights.

    //vector<TLorentzVectorWithPdgId> stackedLeptons;
    //vector<TLorentzVectorWithPdgId> selLeptons;
    //vector<TLorentzVectorWithPdgId> extraLeptons;
    //objectSelection::stackLeptons(stackedLeptons,ElPt, ElEta,  ElPhi,  ElE, MuPt,  MuEta,  MuPhi,  MuE);
    //objectSelection::selectLeptons(stackedLeptons,selLeptons,extraLeptons,ElPt,ElEta, ElPhi, ElE, ElId, ElEtaSc,MuPt, MuEta, MuPhi, MuE, MuId, MuIdTight, MuIdSoft, MuPfIso);
    objectSelection::selectElectrons(selElectrons, extraElectrons, ElPt, ElEta, ElPhi, ElE, ElId, ElEtaSc);
    objectSelection::selectMuons(selMuons, extraMuons, correctedMuPt, MuEta, MuPhi, MuE, MuId, MuIdTight, MuIdSoft, MuPfIso);
    //objectSelection::selectPhotons(selPhotons, PhotPt, PhotEta, PhotPhi, PhotId, PhotScEta, PhotHasPixelSeed, PhotSigmaIetaIeta, selMuons,selElectrons);
    objectSelection::selectJets(selJets, selCentralJets, btags, JetAk04Pt, JetAk04Eta, JetAk04Phi, JetAk04E, JetAk04Id, JetAk04NeutralEmFrac, JetAk04NeutralHadAndHfFrac, JetAk04NeutMult, JetAk04BDiscCisvV2, selMuons, selElectrons, selPhotons);

    //Discriminate ee and mumu
    //int lids = 0;
    //if (selLeptons.size()==2) lids = selLeptons[0].PdgId() * selLeptons[1].PdgId();
    /*bool isEE =   (selLeptons.size()==2 && lids == 121 && !isPhotonDatadriven_); //2 good electrons
    bool isMuMu = (selLeptons.size()==2 && lids == 169 && !isPhotonDatadriven_); //2 good muons
    bool isEMu  = (selLeptons.size()==2 && lids == 143 && !isPhotonDatadriven_);// 1 electron plus 1 muon
    */
    bool isEE = (selElectrons.size()==2 && !isPhotonDatadriven_); //2 good electrons
    bool isMuMu = (selMuons.size()==2 && !isPhotonDatadriven_); //2 good muons
    bool isEMu = (selMuons.size()==1 && selElectrons.size()==1 && !isPhotonDatadriven_);
    bool isGamma = (selPhotons.size() == 1 && isPhotonDatadriven_); //1 good photon


    //###############################################################
    //##################       ANALYSIS CUTS       ##################
    //###############################################################

    
    //if(!isEE && !isMuMu && isEMu&& !isGamma) continue; //not a good lepton pair or photon (if datadriven)
    for(int i =0 ; i < MuPt->size() ; i++) mon.fillHisto("pT_mu","tot",MuPt->at(i),weight);
    for(int i =0 ; i < ElPt->size() ; i++) mon.fillHisto("pT_e","tot",ElPt->at(i),weight);
    mon.fillHisto("nb_mu","tot",MuPt->size(),weight);
    mon.fillHisto("nb_e","tot",ElPt->size(),weight);  

    /*int sel_e,sel_mu,extra_e,extra_mu =0;
    for (int i=0; i<selLeptons.size();i++)
    { 
      if (selLeptons[i].PdgId()== 11) sel_e++;
      if (selLeptons[i].PdgId()== 13) sel_mu++;
    }
    for (int j=0; j<extraLeptons.size();j++)
    {
      if (extraLeptons[j].PdgId()== 11) extra_e++;
      if (extraLeptons[j].PdgId()== 13) extra_mu++;
    }
    */
    mon.fillHisto("nb_mu","sel",selMuons.size(),weight);
    mon.fillHisto("nb_e","sel",selElectrons.size(),weight);
    mon.fillHisto("nb_mu","extra",extraMuons.size(),weight);
    mon.fillHisto("nb_e","extra",extraElectrons.size(),weight);
    if(!isEE && !isMuMu && !isEMu && !isGamma) continue;
    //mon.fillHisto("pT_l1","tot",selLeptons[0].Pt(),weight);
    //mon.fillHisto("pT_l2","tot",selLeptons[1].Pt(),weight); 
    //mon.fillHisto("eta_l1","tot",selLeptons[0].Eta(),weight); 
    //mon.fillHisto("eta_l2","tot",selLeptons[1].Eta(),weight);  
    
    
    if(isEE) currentEvt.s_lepCat = "_ee";
    else if(isMuMu) currentEvt.s_lepCat = "_mumu";
    else if (isEMu) currentEvt.s_lepCat = "_emu";
    if(isMC_&& (fileName.Contains("DY")||fileName.Contains("ZZTo2L")||fileName.Contains("ZZToTauTau"))){
      int GLepId = 1;
      for(int i=0 ; i< GLepBareId->size();i++){
        if(GLepBareMomId->at(i) == 23) GLepId *= fabs(GLepBareId->at(i));
      }
      if (fileName.Contains("DYJetsToTauTau")  &&   GLepId %5 != 0  ) continue ;
      if (fileName.Contains("ZZToTauTau2Nu") && GLepId % 5 != 0 ) continue;
      if (fileName.Contains("ZZToTauTau2Q") &&  GLepId % 5 != 0 ) continue;
      if (fileName.Contains("DYJetsToLL")  &&       GLepId %5 == 0  ) continue ;
      if (fileName.Contains("ZZTo2L2Nu")&&      GLepId % 5 == 0 ) continue;
      if (fileName.Contains("ZZTo2L2Q")&&       GLepId % 5 == 0 ) continue;    
      //cout << " genLep Id product: "<< GLepId << endl;
    }
    //compute and apply the efficiency SFs
    if (isMC_){
    //for leptons
      double weightLeptonsSF = 1.0;
      /*if(isEE) weightLeptonsSF = trigAndIDsfs::diElectronEventSFs(utils::CutVersion::CutSet::Moriond17Cut, selLeptons[0].Pt(), ElEtaSc->at(selLeptons[0].GetIndex()), selLeptons[1].Pt(), ElEtaSc->at(selLeptons[1].GetIndex())) ;
      if(isMuMu) weightLeptonsSF = trigAndIDsfs::diMuonEventSFs( utils::CutVersion::CutSet::Moriond17Cut, MuPt->at(selLeptons[0].GetIndex()), selLeptons[0].Eta(), MuPt->at(selLeptons[1].GetIndex()), selLeptons[1].Eta()) ;
      if(isEMu)  weightLeptonsSF = trigAndIDsfs::EMuEventSFs(utils::CutVersion::CutSet::Moriond17Cut, MuPt->at(selLeptons[1].GetIndex()), selLeptons[1].Eta(), selLeptons[0].Pt(), ElEtaSc->at(selLeptons[0].GetIndex())) ; */
      if(isEE) weightLeptonsSF = trigAndIDsfs::diElectronEventSFs(utils::CutVersion::CutSet::Moriond17Cut, selElectrons[0].Pt(), ElEtaSc->at(selElectrons[0].GetIndex()), selElectrons[1].Pt(), ElEtaSc->at(selElectrons[1].GetIndex())) ;
      if(isMuMu) weightLeptonsSF = trigAndIDsfs::diMuonEventSFs( utils::CutVersion::CutSet::Moriond17Cut, MuPt->at(selMuons[0].GetIndex()), selMuons[0].Eta(), MuPt->at(selMuons[1].GetIndex()), selMuons[1].Eta()) ;
      if(isEMu)  weightLeptonsSF = trigAndIDsfs::EMuEventSFs(utils::CutVersion::CutSet::Moriond17Cut, MuPt->at(selMuons[0].GetIndex()), selMuons[0].Eta(), selElectrons[0].Pt(), ElEtaSc->at(selElectrons[0].GetIndex())) ;
      /*{
        if(fabs(selLeptons[0].PdgId())==11) weightLeptonsSF = trigAndIDsfs::EMuEventSFs(utils::CutVersion::CutSet::Moriond17Cut, MuPt->at(selLeptons[1].GetIndex()), selLeptons[1].Eta(), selLeptons[0].Pt(), ElEtaSc->at(selLeptons[0].GetIndex())) ;
        else if(fabs(selLeptons[0].PdgId())==13) weightLeptonsSF = trigAndIDsfs::EMuEventSFs(utils::CutVersion::CutSet::Moriond17Cut, MuPt->at(selLeptons[0].GetIndex()), selLeptons[0].Eta(), selLeptons[1].Pt(), ElEtaSc->at(selLeptons[1].GetIndex())) ;
      }*/
      weight*=weightLeptonsSF;
    }

    //MET filters
    std::vector<std::pair<int, int> > listMETFilter; //after the passMetFilter function, it contains the bin number of the cut in .first and if it passed 1 or not 0 the METfilter
    bool passMetFilter = utils::passMetFilter(TrigMET, listMETFilter, isMC_);
    mon.fillHisto("metFilters","tot",26,weight); //the all bin, i.e. the last one
    for(unsigned int i =0; i < listMETFilter.size(); i++){
      if(listMETFilter[i].second ==1) mon.fillHisto("metFilters","tot",listMETFilter[i].first,weight);
    }
    //if (!passMetFilter) continue;

    //Avoid double couting for W+jets
    //For some reasons we just have the inclusive sample for the Dilepton region while we have both HT and inclusive samples for the photon region. Hence this cleaning only applies to the photon region.

    //Definition of the relevant analysis variables
    vector<TLorentzVectorWithIndex> selLeptons;
    if(isEE) selLeptons = selElectrons;
    if(isMuMu) selLeptons = selMuons;
    if(isEMu)  
      {
        selLeptons.push_back(selElectrons[0]);
        selLeptons.push_back(selMuons[0]);
      }
    TLorentzVector boson = (isPhotonDatadriven_) ? selPhotons[0] : selLeptons[0] + selLeptons[1];


    /*TLorentzVector boson ;
    if (isEE ) boson = (isPhotonDatadriven_) ? selPhotons[0] : selElectrons[0] + selElectrons[1];
    else if (isMuMu) boson = (isPhotonDatadriven_) ? selPhotons[0] : selMuons[0] + selMuons[1];
    else if (isEMu) boson = (isPhotonDatadriven_) ? selPhotons[0] : selElectrons[0] + selMuons[0]; */
    TLorentzVector METVector; METVector.SetPtEtaPhiE(METPtType1XY->at(0),0.,METPhiType1XY->at(0),METPtType1XY->at(0));

    //Loop on lepton type. This is important also to apply Instr.MET if needed:
    double weightBeforeLoop = weight;
    TLorentzVector bosonBeforeLoop = boson;
    for(unsigned int c = 0; c < tagsR_size; c++){
      weight = weightBeforeLoop;
      boson = bosonBeforeLoop;

      if(!isPhotonDatadriven_){
        if(tagsR[c] == "_ee" && !isEE) continue;
        else if(tagsR[c] == "_mumu" && !isMuMu) continue;
        else if(tagsR[c] == "_emu" && !isEMu) continue;
        else if(tagsR[c] == "_ll" && !(isMuMu || isEE)) continue;
      }

      //Jet category
      //enum {eq0jets,geq1jets,vbf};
      int jetCat = geq1jets;
      if(selJets.size()==0) jetCat = eq0jets;
      else if(utils::passVBFcuts(selJets, boson)) jetCat = vbf;
      currentEvt.s_jetCat = v_jetCat[jetCat];

      //Warning, starting from here ALL plots have to have the currentEvt.s_lepCat in their name, otherwise the reweighting will go crazy
      currentEvt.Fill_evt(v_jetCat[jetCat], tagsR[c], boson, METVector, selJets, EvtRunNum, EvtVtxCnt, EvtFastJetRho, METsig->at(0), selLeptons);

      //mon.fillHisto("jetCategory","tot"+currentEvt.s_lepCat,jetCat,weight);
      //mon.fillHisto("nJets","tot"+currentEvt.s_lepCat,currentEvt.nJets,weight);

      // Compute the btagging efficiency
      if(isMC_) btagger::fill_eff(selCentralJets, btags, JetAk04HadFlav, weight, mon);
      // Apply the btag weights
      if(isMC_) weight *= btagger::apply_sf(selCentralJets, btags, JetAk04HadFlav, btagEffTable, _btag_calibration_reader, syst_);
      
      mon.fillAnalysisHistos(currentEvt, "tot", weight);
      //b veto
      bool passBtag = true;
      for(int i =0 ; i < btags.size() ; i++){
        if (btags[i] > 0.5426) passBtag = false;
      }
      //Phi(jet,MET)
      bool passDeltaPhiJetMET = true;
      for(int i = 0 ; i < selJets.size() ; i++){
        if (fabs(utils::deltaPhi(selJets[i], METVector))<0.5) passDeltaPhiJetMET = false;
      }

      //DPhi
      bool passDphi(currentEvt.deltaPhi_MET_Boson>0.5);
      //boson
      bool passMass(fabs(currentEvt.M_Boson-91)<15);
      bool isZ_SB ( (currentEvt.M_Boson>40  && currentEvt.M_Boson<70) || (currentEvt.M_Boson>110 && currentEvt.M_Boson<200) );
      bool isZ_upSB ( (currentEvt.M_Boson>110 && currentEvt.M_Boson<200) );
      bool passQt (currentEvt.pT_Boson > 55.);
      bool  passThirdLeptonveto ( (isEE && selMuons.size()== 0 ||isMuMu && selElectrons.size() ==0 || isEMu) && extraElectrons.size()==0 && extraMuons.size()== 0);
      //if (isPhotonDatadriven_) passThirdLeptonveto = (selLeptons.size()==0 && extraLeptons.size()==0 );
      
      TString tags = "tot"+currentEvt.s_lepCat; 

      if(currentEvt.M_Boson>40 && currentEvt.M_Boson<200 && passQt && passThirdLeptonveto  && passDeltaPhiJetMET && passDphi){
        if(passBtag)
        {
           if(METVector.Pt()>50 )mon.fillHisto("zmass_bveto50" , tags,currentEvt.M_Boson,weight);
           if(METVector.Pt()>80 )mon.fillHisto("zmass_bveto80" , tags,currentEvt.M_Boson,weight);
           if(METVector.Pt()>125)mon.fillHisto("zmass_bveto125", tags,currentEvt.M_Boson,weight);
           if(passMass)
           {
              mon.fillHisto( "met_Inbveto",tags,METVector.Pt(),weight);
              if(METVector.Pt()>50 )mon.fillHisto("mt_Inbveto50" , tags,currentEvt.MT,weight);
              if(METVector.Pt()>80 )mon.fillHisto("mt_Inbveto80" , tags,currentEvt.MT,weight);
              if(METVector.Pt()>125){
                mon.fillHisto("mt_Inbveto125", tags,currentEvt.MT,weight);
                mon.fillHisto("mT_final"+currentEvt.s_jetCat, tagsR[c].substr(1), currentEvt.MT, weight, divideFinalHistoByBinWidth);
              }
           }
           else if(isZ_SB)
           {
              mon.fillHisto( "met_Outbveto",tags,METVector.Pt(),weight);
              if(METVector.Pt()>50 )mon.fillHisto("mt_Outbveto50" , tags,currentEvt.MT,weight);
              if(METVector.Pt()>80 )mon.fillHisto("mt_Outbveto80" , tags,currentEvt.MT,weight);
              if(METVector.Pt()>125)mon.fillHisto("mt_Outbveto125", tags,currentEvt.MT,weight);
           }
        }
        else
        {
          if(METVector.Pt()>50 )mon.fillHisto("zmass_btag50" , tags,currentEvt.M_Boson,weight);
          if(METVector.Pt()>80 )mon.fillHisto("zmass_btag80" , tags,currentEvt.M_Boson,weight);
          if(METVector.Pt()>125)mon.fillHisto("zmass_btag125", tags,currentEvt.M_Boson,weight);
          if(passMass)
          {
            mon.fillHisto( "met_Inbtag",tags,METVector.Pt(),weight);
            if(METVector.Pt()>50 )mon.fillHisto("mt_Inbtag50" , tags,currentEvt.MT,weight);
            if(METVector.Pt()>80 )mon.fillHisto("mt_Inbtag80" , tags,currentEvt.MT,weight);
            if(METVector.Pt()>125)mon.fillHisto("mt_Inbtag125", tags,currentEvt.MT,weight);
          }
          else if(isZ_SB)
          {
            mon.fillHisto( "met_Outbtag",tags,METVector.Pt(),weight);
            if(METVector.Pt()>50 )mon.fillHisto("mt_Outbtag50" , tags,currentEvt.MT,weight);
            if(METVector.Pt()>80 )mon.fillHisto("mt_Outbtag80" , tags,currentEvt.MT,weight);
            if(METVector.Pt()>125)mon.fillHisto("mt_Outbtag125", tags,currentEvt.MT,weight);
          }
        }

      }
      
      if(currentEvt.M_Boson>40 && currentEvt.M_Boson<200 && passQt && passThirdLeptonveto  && passDeltaPhiJetMET && passDphi)
      {
        for(unsigned int Index=0;Index<optim_Cuts1_met.size();Index++)
        {
          if(METVector.Pt()>optim_Cuts1_met[Index])
          {
            if(passBtag && passMass)mon.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 0.5,weight);
            if(passBtag && isZ_SB)mon.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 1.5,weight);
            if(passBtag && isZ_upSB)mon.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 2.5,weight);
            if(!passBtag && passMass)mon.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 3.5,weight);
            if(!passBtag && isZ_SB)mon.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 4.5,weight);
            if(!passBtag && isZ_upSB)mon.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 5.5,weight);
          }
        }
      }
      mon.fillHisto("eventflow","tot",1,weight);
      mon.fillHisto("eventflow",tags,1,weight);
      if(!passMass) continue;
      mon.fillHisto("eventflow","tot",2,weight);
      mon.fillHisto("eventflow",tags,2,weight);
      if(!passQt) continue;
      mon.fillHisto("eventflow","tot",3,weight);
      mon.fillHisto("eventflow",tags,3,weight);
       //Phi(Z,MET)
      if(!passThirdLeptonveto) continue;
      mon.fillHisto("eventflow","tot",4,weight);
      mon.fillHisto("eventflow",tags,4,weight);
      if(!passBtag) continue;
      mon.fillHisto("eventflow","tot",5,weight);
      mon.fillHisto("eventflow",tags,5,weight);
      if(!passDeltaPhiJetMET) continue;
      mon.fillHisto("eventflow","tot",6,weight);
      mon.fillHisto("eventflow",tags,6,weight);
      if(!passDphi) continue;
      mon.fillHisto("eventflow","tot",7,weight);
      mon.fillHisto("eventflow",tags,7,weight);

      mon.fillAnalysisHistos(currentEvt, "beforeMETcut", weight);
      //mon.fillHisto("reco-vtx","beforeMETcut"+currentEvt.s_lepCat,EvtVtxCnt,weight);
      //mon.fillHisto("jetCategory","beforeMETcut"+currentEvt.s_lepCat,jetCat,weight);

      //MET>80
      if(METVector.Pt()<80) continue;
      mon.fillHisto("eventflow","tot",8,weight);
      mon.fillHisto("eventflow",tags,8,weight);
      //MET>125
      if(isMC_&& isEMu &&fileName.Contains("WJetsToLNu")){
        int ElIndex = selElectrons.at(0).GetIndex();
        int MuIndex = selMuons.at(0).GetIndex();
        ofstream Evt;
        Evt.open("WJetsAnalysis.txt", ofstream::out | ofstream::app);
        Evt<< "************Event-begin************"<<endl;
        Evt<<"Electron:"<<endl;
        Evt<<"Mother Particle PdgId:"<<GLepBareMomId->at(ElIndex)<<endl;
        Evt<<"LorentzVector:"<<selElectrons.at(0).Px()<<"\t"<<selElectrons.at(0).Py()<<"\t"<<selElectrons.at(0).Pz()<<"\t"<<selElectrons.at(0).E()<<endl;
        Evt<<"Muon:"<<endl;
        Evt<<"Mother Particle PdgId:"<<GLepBareMomId->at(MuIndex)<<endl;
        Evt<<"LorentzVector:"<<selMuons.at(0).Px()<<"\t"<<selMuons.at(0).Py()<<"\t"<<selMuons.at(0).Pz()<<"\t"<<selMuons.at(0).E()<<endl;
        Evt<<"MET:"<<METVector.Pt()<<endl;
        Evt<< "*************Event-end*************"<<endl;
        //cout<< Evt.rdbuf()<<endl;
        Evt.close();
      }
      if(METVector.Pt()<125) continue;
      mon.fillHisto("eventflow","tot",9,weight);
      mon.fillHisto("eventflow",tags,9,weight);
      //###############################################################
      //##################     END OF SELECTION      ##################
      //###############################################################
      //mon.fillHisto("reco-vtx","final"+currentEvt.s_lepCat,EvtVtxCnt,weight);
      //mon.fillHisto("jetCategory","final"+currentEvt.s_lepCat,jetCat,weight);
      //mon.fillAnalysisHistos(currentEvt, "final", weight);

    }
  }

  //###############################################################
  //##################        END OF LOOP        ##################
  //###############################################################

  TFile* outFile=TFile::Open(outputFile_,"recreate");
  mon.WriteForSysts(syst_,keepAllControlPlots_);
  outFile->Close();

}
