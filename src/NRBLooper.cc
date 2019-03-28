#define HZZ2l2nuLooper_cxx

#include <ctime>
#include <TH1.h>

#include <BTagWeight.h>
#include <ElectronBuilder.h>
#include <EWCorrectionWeight.h>
#include <GenJetBuilder.h>
#include <JetBuilder.h>
#include <LeptonsEfficiencySF.h>
#include <LooperMain.h>
#include <MuonBuilder.h>
#include <ObjectSelection.h>
#include <PhotonBuilder.h>
#include <PhotonEfficiencySF.h>
#include <PileUpWeight.h>
#include <SmartSelectionMonitor.h>
#include <SmartSelectionMonitor_hzz.h>
#include <Trigger.h>
#include <Utils.h>

#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <algorithm>


bool LooperMain::passTrigger(int triggerType){
  bool passEventSelection = false;
  int triggerWeight=0;
  triggerWeight = trigger::passTrigger(triggerType, *TrigHltDiMu, *TrigHltMu, *TrigHltDiEl, *TrigHltEl, *TrigHltElMu, *TrigHltPhot, TrigHltDiMu_prescale, TrigHltMu_prescale, TrigHltDiEl_prescale, TrigHltEl_prescale, TrigHltElMu_prescale, TrigHltPhot_prescale);
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

  ElectronBuilder electronBuilder{fReader, options_};
  MuonBuilder muonBuilder{fReader, options_, randomGenerator_};

  PhotonBuilder photonBuilder{fReader, options_};
  photonBuilder.EnableCleaning({&muonBuilder, &electronBuilder});

  GenJetBuilder genJetBuilder{fReader, options_};
  JetBuilder jetBuilder{fReader, options_, randomGenerator_};
  jetBuilder.EnableCleaning({&muonBuilder, &electronBuilder, &photonBuilder});
  jetBuilder.SetGenJetBuilder(&genJetBuilder);

  EWCorrectionWeight ewCorrectionWeight(fReader, options_);
  BTagWeight bTagWeight(options_);
  PileUpWeight pileUpWeight;

  SmartSelectionMonitor_hzz mon;
  mon.declareHistos_NRB();

  cout << "nb of entries in the input file =" << nentries << endl;

  cout << "fileName is " << fileName << endl;

  std::vector<string> v_jetCat = {"_eq0jets","_geq1jets","_vbf"};
  std::vector<string> tagsR = {"_ee", "_mumu", "_ll","_emu"};
  unsigned int tagsR_size =  tagsR.size();  

  enum {ee, mumu, ll, emu,lepCat_size};
  enum {eq0jets, geq1jets, vbf, jetCat_size};
 
  // ***--- Instr. MET building ---*** \\
  //Compute once weights for Instr. MET reweighting if needed
  string const base_path = string(getenv("HZZ2L2NU_BASE")) + "/";
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
    fReader.SetEntry(jentry);

    std::time_t currentTime = std::time(nullptr);
    if(jentry % 10000 ==0) cout << jentry << " of " << nentries << ". It is now " << std::asctime(std::localtime(&currentTime));
    

    evt currentEvt;
    
    double weight = 1.;
    double totEventWeight = 1.;
    //get the MC event weight if exists
    if (isMC_) {
      //get the MC event weight if exists
      weight *= (EvtWeights.GetSize()>0 ? EvtWeights[0] : 1);
      if ((sumWeightInBonzai_>0)&&(sumWeightInBaobab_>0)) totEventWeight = weight*sumWeightInBaobab_/sumWeightInBonzai_;
      //get the PU weights
      float weightPU = pileUpWeight(*EvtPuCntTruth); 
      weight *= weightPU;
    }
    else {
      totEventWeight = totalEventsInBaobab_/nentries;
    }
    // Remove events with 0 vtx
    if(*EvtVtxCnt == 0 ) continue;

    mon.fillHisto("totEventInBaobab","tot",*EvtPuCnt,totEventWeight);
    if(runOnBaobabs_){
      if (!LooperMain::passTrigger(triggerType)) continue;
    }

    mon.fillHisto("eventflow","tot",0,weight);



    //###############################################################
    //##################     OBJECT CORRECTIONS    ##################
    //###############################################################
    // electroweak corrections
    weight *= ewCorrectionWeight();



    //###############################################################
    //##################     OBJECT SELECTION      ##################
    //###############################################################

    auto const &tightElectrons = electronBuilder.GetTight();
    auto const &looseElectrons = electronBuilder.GetLoose();

    auto const &tightMuons = muonBuilder.GetTight();
    auto const &looseMuons = muonBuilder.GetLoose();

    auto const &photons = photonBuilder.Get();
    auto const &jets = jetBuilder.Get();

    //vector<TLorentzVectorWithPdgId> stackedLeptons;
    //vector<TLorentzVectorWithPdgId> selLeptons;
    //vector<TLorentzVectorWithPdgId> extraLeptons;
    //objectSelection::stackLeptons(stackedLeptons,ElPt, ElEta,  ElPhi,  ElE, MuPt,  MuEta,  MuPhi,  MuE);
    //objectSelection::selectLeptons(stackedLeptons,selLeptons,extraLeptons,ElPt,ElEta, ElPhi, ElE, ElId, ElEtaSc,MuPt, MuEta, MuPhi, MuE, MuId, MuIdTight, MuIdSoft, MuPfIso);
    //objectSelection::selectPhotons(photons, PhotPt, PhotEta, PhotPhi, PhotId, PhotScEta, PhotHasPixelSeed, PhotSigmaIetaIeta, selMuons, tightElectrons);

    //Discriminate ee and mumu
    //int lids = 0;
    //if (selLeptons.size()==2) lids = selLeptons[0].PdgId() * selLeptons[1].PdgId();
    /*bool isEE =   (selLeptons.size()==2 && lids == 121 && !isPhotonDatadriven_); //2 good electrons
    bool isMuMu = (selLeptons.size()==2 && lids == 169 && !isPhotonDatadriven_); //2 good muons
    bool isEMu  = (selLeptons.size()==2 && lids == 143 && !isPhotonDatadriven_);// 1 electron plus 1 muon
    */
    bool isEE = (tightElectrons.size() >= 2 and not isPhotonDatadriven_); //2 good electrons
    bool isMuMu = (tightMuons.size() >= 2 and not isPhotonDatadriven_); //2 good muons
    bool isEMu = (tightMuons.size() == 1 and tightElectrons.size() == 1 and
      not isPhotonDatadriven_);
    bool isGamma = (photons.size() == 1 && isPhotonDatadriven_); //1 good photon


    //###############################################################
    //##################       ANALYSIS CUTS       ##################
    //###############################################################

    
    //if(!isEE && !isMuMu && isEMu&& !isGamma) continue; //not a good lepton pair or photon (if datadriven)
    for(int i =0 ; i < MuPt.GetSize() ; i++) mon.fillHisto("pT_mu","tot",MuPt[i],weight);
    for(int i =0 ; i < ElPt.GetSize() ; i++) mon.fillHisto("pT_e","tot",ElPt[i],weight);
    mon.fillHisto("nb_mu","tot",MuPt.GetSize(),weight);
    mon.fillHisto("nb_e","tot",ElPt.GetSize(),weight);  

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
    
    mon.fillHisto("nb_mu", "sel",
                  std::min<int>(tightMuons.size(), 2), weight);
    mon.fillHisto("nb_e", "sel",
                  std::min<int>(tightElectrons.size(), 2), weight);
    mon.fillHisto("nb_mu", "extra",
                  looseMuons.size() - std::min<int>(tightMuons.size(), 2),
                  weight);
    mon.fillHisto("nb_e", "extra",
                  looseElectrons.size() -
                    std::min<int>(tightElectrons.size(), 2),
                  weight);

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
      for(int i=0 ; i< GLepBareId.GetSize();i++){
        if(GLepBareMomId[i] == 23) GLepId *= fabs(GLepBareId[i]);
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
    if (isMC_) {
      if (isEMu)
        weight *= trigAndIDsfs::EMuEventSFs(
          utils::CutVersion::CutSet::Moriond17Cut,
          tightMuons[0].uncorrP4.Pt(), tightMuons[0].uncorrP4.Eta(),
          tightElectrons[0].p4.Pt(), tightElectrons[0].etaSc);
      else if (isMuMu)
        weight *= trigAndIDsfs::diMuonEventSFs(
          utils::CutVersion::CutSet::Moriond17Cut,
          tightMuons[0].uncorrP4.Pt(), tightMuons[0].uncorrP4.Eta(),
          tightMuons[1].uncorrP4.Pt(), tightMuons[1].uncorrP4.Eta());
      else if (isEE)
        weight *= trigAndIDsfs::diElectronEventSFs(
          utils::CutVersion::CutSet::Moriond17Cut,
          tightElectrons[0].p4.Pt(), tightElectrons[0].etaSc,
          tightElectrons[1].p4.Pt(), tightElectrons[1].etaSc);
    }

    //MET filters
    std::vector<std::pair<int, int> > listMETFilter; //after the passMetFilter function, it contains the bin number of the cut in .first and if it passed 1 or not 0 the METfilter
    bool passMetFilter = utils::passMetFilter(*TrigMET, listMETFilter, isMC_);
    mon.fillHisto("metFilters","tot",26,weight); //the all bin, i.e. the last one
    for(unsigned int i =0; i < listMETFilter.size(); i++){
      if(listMETFilter[i].second ==1) mon.fillHisto("metFilters","tot",listMETFilter[i].first,weight);
    }
    //if (!passMetFilter) continue;

    //Avoid double couting for W+jets
    //For some reasons we just have the inclusive sample for the Dilepton region while we have both HT and inclusive samples for the photon region. Hence this cleaning only applies to the photon region.
    
    //Definition of the relevant analysis variables
    std::vector<Lepton> tightLeptons;

    if (isMuMu or isEMu) {
      for (auto const &mu : tightMuons)
        tightLeptons.emplace_back(mu);
    } else if (isEE or isEMu) {
      for (auto const &e : tightElectrons)
        tightLeptons.emplace_back(e);
    }

    if (isEMu)
      std::sort(tightLeptons.begin(), tightLeptons.end(), PtOrdered);

    TLorentzVector boson = (isPhotonDatadriven_) ? photons[0].p4 :
      tightLeptons[0].p4 + tightLeptons[1].p4;

    TLorentzVector METVector; METVector.SetPtEtaPhiE(METPtType1XY[0],0.,METPhiType1XY[0],METPtType1XY[0]);

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

      if (jets.size() == 0)
        jetCat = eq0jets;
      else if (utils::PassVbfCuts(jets, boson))
        jetCat = vbf;

      currentEvt.s_jetCat = v_jetCat[jetCat];

      //Warning, starting from here ALL plots have to have the currentEvt.s_lepCat in their name, otherwise the reweighting will go crazy
      currentEvt.Fill_evt(
        v_jetCat[jetCat], tagsR[c], boson, METVector, jets, *EvtRunNum,
        *EvtVtxCnt, *EvtFastJetRho, METsig[0], tightLeptons);

      //mon.fillHisto("jetCategory","tot"+currentEvt.s_lepCat,jetCat,weight);
      //mon.fillHisto("nJets","tot"+currentEvt.s_lepCat,currentEvt.nJets,weight);

      // Apply the btag weights
      if (isMC_)
        weight *= bTagWeight(jets);
      
      mon.fillAnalysisHistos(currentEvt, "tot", weight);

      // b veto
      bool passBTag = true;

      for (auto const &jet : jets)
        if (jet.bTagCsvV2 > 0.5426 and std::abs(jet.p4.Eta()) < 2.5) {
          passBTag = false;
          break;
        }

      if (not passBTag)
        continue;

      // Phi(jet,MET)
      bool passDeltaPhiJetMET = true;

      for (auto const &jet : jets)
        if (std::abs(utils::deltaPhi(jet.p4, METVector)) < 0.5) {
          passDeltaPhiJetMET = false;
          break;
        }

      if (not passDeltaPhiJetMET)
        continue;

      //DPhi
      bool passDphi(currentEvt.deltaPhi_MET_Boson>0.5);
      //boson
      bool passMass(fabs(currentEvt.M_Boson-91)<15);
      bool isZ_SB ( (currentEvt.M_Boson>40  && currentEvt.M_Boson<70) || (currentEvt.M_Boson>110 && currentEvt.M_Boson<200) );
      bool isZ_upSB ( (currentEvt.M_Boson>110 && currentEvt.M_Boson<200) );
      bool passQt (currentEvt.pT_Boson > 55.);

      unsigned const numExtraLeptons =
        looseMuons.size() - std::min<unsigned>(tightMuons.size(), 2);
        looseElectrons.size() - std::min<unsigned>(tightElectrons.size(), 2);
      bool passThirdLeptonveto = (isEE and tightMuons.empty() or
        isMuMu and tightElectrons.empty() or isEMu) and
        numExtraLeptons == 0;
      
      TString tags = "tot"+currentEvt.s_lepCat; 

      if(currentEvt.M_Boson>40 && currentEvt.M_Boson<200 && passQt && passThirdLeptonveto  && passDeltaPhiJetMET && passDphi){
        if(passBTag)
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
            if(passBTag && passMass)mon.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 0.5,weight);
            if(passBTag && isZ_SB)mon.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 1.5,weight);
            if(passBTag && isZ_upSB)mon.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 2.5,weight);
            if(!passBTag && passMass)mon.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 3.5,weight);
            if(!passBTag && isZ_SB)mon.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 4.5,weight);
            if(!passBTag && isZ_upSB)mon.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 5.5,weight);
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
      if(!passBTag) continue;
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
