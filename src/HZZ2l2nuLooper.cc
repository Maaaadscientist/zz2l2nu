#define HZZ2l2nuLooper_cxx

#include <limits>

#include <BTagWeight.h>
#include <ElectronBuilder.h>
#include <EWCorrectionWeight.h>
#include <LeptonsEfficiencySF.h>
#include <LooperMain.h>
#include <MuonBuilder.h>
#include <ObjectSelection.h>
#include <PhotonBuilder.h>
#include <PhotonEfficiencySF.h>
#include <PileUpWeight.h>
#include <SmartSelectionMonitor.h>
#include <Trigger.h>
#include <Utils.h>

#include <ctime>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <algorithm>




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

  ElectronBuilder electronBuilder{fReader, options_};
  MuonBuilder muonBuilder{fReader, options_, randomGenerator_};

  PhotonBuilder photonBuilder{fReader, options_};
  photonBuilder.EnableCleaning(&electronBuilder);

  EWCorrectionWeight ewCorrectionWeight(fReader, options_);
  BTagWeight bTagWeight(options_);
  PileUpWeight pileUpWeight;

  SmartSelectionMonitor_hzz mon;
  mon.declareHistos();

  cout << "nb of entries in the input file =" << nentries << endl;

  cout << "fileName is " << fileName << endl;

  enum {ee, mumu, ll, lepCat_size};
  enum {eq0jets, geq1jets, vbf, jetCat_size};
  std::vector<string> v_jetCat = {"_eq0jets","_geq1jets","_vbf"};
  std::vector<string> tagsR = {"_ee", "_mumu", "_ll"};
  unsigned int tagsR_size =  tagsR.size();  

  //Definition of the final histos (and in particular of the mT binning
  bool divideFinalHistoByBinWidth = false; //For final plots, we don't divide by the bin width to ease computations of the yields by eye.
  std::vector<TH1*> h_mT(jetCat_size); std::vector<int> h_mT_size(jetCat_size);
  h_mT[eq0jets] = (TH1*) mon.getHisto("mT_final_eq0jets", tagsR[ee].substr(1), divideFinalHistoByBinWidth); h_mT_size[eq0jets] = h_mT[eq0jets]->GetNbinsX();
  h_mT[geq1jets] = (TH1*) mon.getHisto("mT_final_geq1jets", tagsR[ee].substr(1), divideFinalHistoByBinWidth); h_mT_size[geq1jets] = h_mT[geq1jets]->GetNbinsX();
  h_mT[vbf] = (TH1*) mon.getHisto("mT_final_vbf", tagsR[ee].substr(1), divideFinalHistoByBinWidth); h_mT_size[vbf] = h_mT[vbf]->GetNbinsX();
  mon.getHisto("mT_final_eq0jets", tagsR[mumu].substr(1), divideFinalHistoByBinWidth); //The .substr(1) removes the annoying _ in the tagsR definition.
  mon.getHisto("mT_final_geq1jets", tagsR[mumu].substr(1), divideFinalHistoByBinWidth);
  mon.getHisto("mT_final_vbf", tagsR[mumu].substr(1), divideFinalHistoByBinWidth);
  int h_mT_maxSize = std::max({h_mT_size[eq0jets], h_mT_size[geq1jets], h_mT_size[vbf]}); 

  // Histograms for PDF replicas
  vector<vector<vector<TH1F*>>> pdfReplicas;
  if(syst_=="pdf_up" || syst_=="pdf_down"){
    for(unsigned int i = 0 ; i < v_jetCat.size() ; i++){
      vector<vector<TH1F*>> currentVectorOfVectors;
      for(unsigned int j = 0 ; j < tagsR.size() -1 ; j++){
        vector<TH1F*> currentVectorOfReplicas;
        for(unsigned int k = 0 ; k < 100 ; k++){
          currentVectorOfReplicas.push_back((TH1F*) mon.getHisto("mT_final"+v_jetCat.at(i), tagsR.at(j).substr(1), divideFinalHistoByBinWidth)->Clone((TString) "mT"+v_jetCat.at(i)+tagsR.at(j)+"_"+to_string(k)));
        }
        currentVectorOfVectors.push_back(currentVectorOfReplicas);
      }
      pdfReplicas.push_back(currentVectorOfVectors);
    }
  }

  // ***--- Instr. MET building ---*** \\
  //Compute once weights for Instr. MET reweighting if needed
  string const base_path = string(getenv("HZZ2L2NU_BASE")) + "/";
  string weightFileType = "InstrMET";
  bool weight_NVtx_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root");
  bool weight_Pt_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root");
  bool weight_Mass_exist = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_lineshape_mass.root");
  std::map<TString, std::map<double, std::pair<double, double> > > NVtxWeight_map, PtWeight_map;
  std::map<TString, TH1D*> LineshapeMassWeight_map;
  TH1 *h_Vtx = (TH1*) mon.getHisto("reco-vtx", "toGetBins", false); int h_Vtx_size = h_Vtx->GetNbinsX();
  TH1 *h_pT = (TH1*) mon.getHisto("pT_Boson", "toGetBins", false); int h_pT_size = h_pT->GetNbinsX();
  std::vector<std::vector<std::vector<std::vector<std::vector<std::pair<double, double> > > > > > mT_InstrMET_map(lepCat_size, std::vector<std::vector<std::vector<std::vector<std::pair<double, double> > > > >(jetCat_size, std::vector<std::vector<std::vector<std::pair<double, double> > > >(h_mT_maxSize+1, std::vector<std::vector<std::pair<double, double> > >(h_Vtx_size+1, std::vector<std::pair<double, double> >(h_pT_size+1/*, std::pair<double, double>*/)))));
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


  //###############################################################
  //##################     EVENT LOOP STARTS     ##################
  //###############################################################

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if ((jentry>maxEvents_)&&(maxEvents_>=0)) break;
    fReader.SetEntry(jentry);

    std::time_t currentTime = std::time(nullptr);
    if(jentry % 10000 ==0) cout << jentry << " of " << nentries << ". It is now " << std::asctime(std::localtime(&currentTime));

    evt currentEvt;

    double theRandomNumber = randomGenerator_.Rndm(); //Used for the uncertainty on the 3rd lepton veto.

    double weight = 1.;
    double totEventWeight = 1.;
    //get the MC event weight if exists
    if (isMC_) {
      //get the MC event weight if exists
      if(EvtWeights.GetSize()>1)
        weight *= (EvtWeights.GetSize()>0 ? EvtWeights[1] : 1); //Value 0 is not filled properly for LO generated samples (MadgraphMLM)
      if ((sumWeightInBonzai_>0)&&(sumWeightInBaobab_>0)) totEventWeight = weight*sumWeightInBaobab_/sumWeightInBonzai_;
      if (jentry == 0){
        std::cout<< "Printing once the content of EvtWeights for event " << jentry << ":" << std::endl;
        if(EvtWeights.GetSize()>1) for(unsigned int i = 0; i < EvtWeights.GetSize(); i++ ) std::cout<< i << " " << EvtWeights[i] << std::endl;
      }
      //get the PU weights
      float weightPU = pileUpWeight(*EvtPuCntTruth); 
      weight *= weightPU;
    }
    else {
      totEventWeight = totalEventsInBaobab_/nentries;
    }

    mon.fillHisto("totEventInBaobab","tot",*EvtPuCnt,totEventWeight);
    mon.fillHisto("eventflow","tot",0,weight);

    // Remove events with 0 vtx
    if(*EvtVtxCnt == 0 ) continue;

    for(int i =0 ; i < MuPt.GetSize() ; i++) mon.fillHisto("pT_mu","tot",MuPt[i],weight);
    for(int i =0 ; i < ElPt.GetSize() ; i++) mon.fillHisto("pT_e","tot",ElPt[i],weight);
    mon.fillHisto("nb_mu","tot",MuPt.GetSize(),weight);
    mon.fillHisto("nb_e","tot",ElPt.GetSize(),weight);
    mon.fillHisto("pile-up","tot",*EvtPuCnt,weight);
    mon.fillHisto("truth-pile-up","tot",*EvtPuCntTruth,weight);
    mon.fillHisto("reco-vtx","tot",*EvtVtxCnt,weight);


    //###############################################################
    //##################     OBJECT CORRECTIONS    ##################
    //###############################################################
    // electroweak corrections
    weight *= ewCorrectionWeight();

    // Theory uncertainties
    double thUncWeight = 1.;
    if(syst_ !="") thUncWeight = utils::getTheoryUncertainties(EvtWeights, syst_);
    if(thUncWeight == 0) continue; // There are some rare cases where a weight is at 0, making an indeterminate form (0/0) in the code. I figured out it was an easy (albeit a bit coward) way to avoid it without changing all the code for an effect of less than 0.01%.
    weight *= thUncWeight;

    //###############################################################
    //##################     OBJECT SELECTION      ##################
    //###############################################################

    auto const &tightElectrons = electronBuilder.GetTightElectrons();
    auto const &looseElectrons = electronBuilder.GetLooseElectrons();

    auto const &tightMuons = muonBuilder.GetTightMuons();
    auto const &looseMuons = muonBuilder.GetLooseMuons();

    auto const &photons = photonBuilder.GetPhotons();

    vector<TLorentzVectorWithIndex> selJets; //Jets passing Id and cleaning, with |eta|<4.7 and pT>30GeV. Used for jet categorization and deltaPhi cut.
    vector<TLorentzVectorWithIndex> selCentralJets; //Same as the previous one, but with tracker acceptance (|eta| <= 2.5). Used to compute btag efficiency and weights. 
    vector<double> btags; //B-Tag discriminant, recorded for selCentralJets. Used for b-tag veto, efficiency and weights.

    objectSelection::selectJets(selJets, selCentralJets, btags, JetAk04Pt, JetAk04Eta, JetAk04Phi, JetAk04E, JetAk04Id, JetAk04NeutralEmFrac, JetAk04NeutralHadAndHfFrac, JetAk04NeutMult, JetAk04BDiscCisvV2, tightMuons, tightElectrons, photons);

    //Discriminate ee and mumu
    bool isEE = (tightElectrons.size() >= 2 && !isPhotonDatadriven_); //2 good electrons
    bool isMuMu = (tightMuons.size() >= 2 && !isPhotonDatadriven_); //2 good muons
    bool isGamma = (photons.size() == 1 && isPhotonDatadriven_); //1 good photon

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

    //###############################################################
    //##################       ANALYSIS CUTS       ##################
    //###############################################################


    if(!isEE && !isMuMu && !isGamma) continue; //not a good lepton pair or photon (if datadriven)
    mon.fillHisto("eventflow","tot",1,weight);

    if(isEE) currentEvt.s_lepCat = "_ee";
    else if(isMuMu) currentEvt.s_lepCat = "_mumu";

    //compute and apply the efficiency SFs
    if (isMC_) {
      if (not isPhotonDatadriven_) {  // Leptons
        if (isEE)
          weight *= trigAndIDsfs::diElectronEventSFs(
            utils::CutVersion::CutSet::Moriond17Cut,
            tightElectrons[0].p4.Pt(), tightElectrons[0].etaSc,
            tightElectrons[1].p4.Pt(), tightElectrons[1].etaSc);
        else
          weight *= trigAndIDsfs::diMuonEventSFs(
            utils::CutVersion::CutSet::Moriond17Cut,
            tightMuons[0].uncorrP4.Pt(), tightMuons[0].uncorrP4.Eta(),
            tightMuons[1].uncorrP4.Pt(), tightMuons[1].uncorrP4.Eta());
      }
      else {  // Photons
        PhotonEfficiencySF phoEff;
        weight *= phoEff.getPhotonEfficiency(
          photons[0].p4.Pt(), photons[0].etaSc, "tight",
          utils::CutVersion::Moriond17Cut).first;
      }
    }

    //trigger weights for photon data
    if(isPhotonDatadriven_){
      int triggerWeight =0, triggerType;
      if(isMC_) triggerType = trigger::MC_Photon;
      else triggerType = trigger::SinglePhoton;

      triggerWeight = trigger::passTrigger(triggerType, *TrigHltDiMu, *TrigHltMu, *TrigHltDiEl, *TrigHltEl, *TrigHltElMu, *TrigHltPhot, TrigHltDiMu_prescale, TrigHltMu_prescale, TrigHltDiEl_prescale, TrigHltEl_prescale, TrigHltElMu_prescale, TrigHltPhot_prescale, photons[0].p4.Pt());
      if(triggerWeight==0) continue; //trigger not found
      weight *= triggerWeight;
    }

    //MET filters
    std::vector<std::pair<int, int> > listMETFilter; //after the passMetFilter function, it contains the bin number of the cut in .first and if it passed 1 or not 0 the METfilter
    bool passMetFilter = utils::passMetFilter(*TrigMET, listMETFilter, isMC_);
    //now fill the metFilter eventflow
    mon.fillHisto("metFilters","tot",26,weight); //the all bin, i.e. the last one
    for(unsigned int i =0; i < listMETFilter.size(); i++){
      if(listMETFilter[i].second ==1) mon.fillHisto("metFilters","tot",listMETFilter[i].first,weight);
    }
    // if (!passMetFilter) continue; //This cut is removed from now because of potential bugs and no expected impact on the 2016 results. We will see after if we re-apply it.

    //Avoid double counting for W+jets
    //For some reasons we just have the inclusive sample for the Dilepton region while we have both HT and inclusive samples for the photon region. Hence this cleaning only applies to the photon region.
    if(isPhotonDatadriven_){
      if (isMC_Wlnu_inclusive || isMC_Wlnu_HT100){ //Avoid double counting and make our W#rightarrow l#nu exclusif of the dataset with a cut on HT...
        bool isHT100 = false;
        //Let's create our own gen HT variable
        double vHT =0;
        TLorentzVector genJet_uncleaned;
        for(size_t ig=0; ig<GJetAk04Pt.GetSize(); ig++){
          genJet_uncleaned.SetPtEtaPhiE(GJetAk04Pt[ig], GJetAk04Eta[ig], GJetAk04Phi[ig], GJetAk04E[ig]);
          double minDRlj = std::numeric_limits<double>::infinity();

          for (auto const &l : tightMuons)
            minDRlj = std::min(minDRlj, utils::deltaR(genJet_uncleaned, l.p4));

          for (auto const &l : tightElectrons)
            minDRlj = std::min(minDRlj, utils::deltaR(genJet_uncleaned, l.p4));

          if (minDRlj < 0.4)
            continue;

          vHT += GJetAk04Pt[ig];
        }
        if(vHT >100) isHT100 = true;
        if(isMC_Wlnu_inclusive) mon.fillHisto("custom_HT","forWlnu_inclusive",vHT,weight);
        if(isMC_Wlnu_HT100) mon.fillHisto("custom_HT","forWlnu_HT100",vHT,weight);
        if(isMC_Wlnu_inclusive && isHT100) continue; //reject event
        if(isMC_Wlnu_HT100 && !isHT100) continue; //reject event
      }

      //Avoid double counting for NLO ZvvG:
      if( isMC_NLO_ZGTo2NuG_inclusive && photons[0].p4.Pt() >= 130) continue;
      if( isMC_NLO_ZGTo2NuG_Pt130 && photons[0].p4.Pt() < 130) continue;

    }

    //Definition of the relevant analysis variables
    std::vector<Lepton> tightLeptons;

    if (isMuMu) {
      for (auto const &mu : tightMuons)
        tightLeptons.emplace_back(mu);
    } else if (isEE) {
      for (auto const &e : tightElectrons)
        tightLeptons.emplace_back(e);
    }

    TLorentzVector boson = (isPhotonDatadriven_) ? photons[0].p4 :
      tightLeptons[0].p4 + tightLeptons[1].p4;

    TLorentzVector METVector; METVector.SetPtEtaPhiE(METPtType1XY[0],0.,METPhiType1XY[0],METPtType1XY[0]);
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
        itlow = NVtxWeight_map[tagsR[c]].upper_bound(*EvtVtxCnt); //look at which bin in the map currentEvt.nVtx corresponds
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
      currentEvt.Fill_evt(
        v_jetCat[jetCat], tagsR[c], boson, METVector, selJets, *EvtRunNum,
        *EvtVtxCnt, *EvtFastJetRho, METsig[0], tightLeptons);

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


      //3rd lepton veto (with uncertainty)
      bool passLeptonVeto = true;

      if (isPhotonDatadriven_)
        passLeptonVeto = (looseElectrons.empty() and looseMuons.empty());

      if (isEE)
        passLeptonVeto = (looseElectrons.size() == 2 and looseMuons.empty());

      if (isMuMu)
        passLeptonVeto = (looseElectrons.empty() and looseMuons.size() == 2);

      if (not passLeptonVeto and
          (syst_ == "lepveto_up" or syst_ == "lepveto_down")) {
        int numExtraLeptons;

        if (isPhotonDatadriven_)
          numExtraLeptons = std::min<int>(tightElectrons.size(), 2) +
            std::min<int>(tightMuons.size(), 2);
        else
          numExtraLeptons = looseElectrons.size() + looseMuons.size() - 2;
        
        if (theRandomNumber < std::pow(0.04, numExtraLeptons))
          // Estimate a 4% uncertainty
          passLeptonVeto = true;
        
        if (syst_ == "lepveto_down")
          weight *= -1.;
      }

      if (not passLeptonVeto)
        continue;
      
      if (currentEvt.s_lepCat == "_ll")
        mon.fillHisto("eventflow", "tot", 5, weight);

      // Compute the btagging efficiency
      if (isMC_)
        FillBTagEfficiency(selCentralJets, btags, JetAk04HadFlav, weight, mon);

      //b veto
      bool passBTag = true;
      for(int i =0 ; i < btags.size() ; i++){
        if (btags[i] > 0.5426) passBTag = false;
      }
      if(!passBTag) continue;

      // Apply the btag weights
      if (isMC_) {
        double const w = bTagWeight(selCentralJets, btags, JetAk04HadFlav);
        weight *= w;

        if (currentEvt.s_lepCat == "_ll")
          mon.fillProfile("BTagWeightvsMT", "TEST", currentEvt.MT, w, weight); //FIXME remove after
      }

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
      mon.fillHisto("jetCategory","finalSelection"+currentEvt.s_lepCat,jetCat,weight, divideFinalHistoByBinWidth);
      mon.fillAnalysisHistos(currentEvt, "finalSelection", weight, divideFinalHistoByBinWidth);
      //THE RESULT:
      mon.fillHisto("mT_final"+currentEvt.s_jetCat, tagsR[c].substr(1), currentEvt.MT, weight, divideFinalHistoByBinWidth);
 
      //Uncertainties
      if(syst_ == "QCDscale_up" || syst_ == "QCDscale_down"){
        mon.fillHisto("mT_finalBinning0j"+currentEvt.s_jetCat, tagsR[c].substr(1), currentEvt.MT, weight, divideFinalHistoByBinWidth);
        mon.fillHisto("mT_finalBinning0j"+currentEvt.s_jetCat, tagsR[c].substr(1)+"_nominal", currentEvt.MT, weight/thUncWeight, divideFinalHistoByBinWidth);
      }
      if((syst_ == "pdf_up" || syst_ == "pdf_down") && currentEvt.s_lepCat != "_ll" && EvtWeights.GetSize() >= 110){
        for(int i = 0 ; i < 100 ; i++){
          pdfReplicas.at(jetCat).at(lepCat).at(i)->Fill(currentEvt.MT,weight*EvtWeights[i+10]/EvtWeights[1]);
        }
      }


      //Prepare the correct computation of the stat uncertainty for the mT plots with Instr.MET.
      if(isPhotonDatadriven_){
        int mT = min(h_mT_size[jetCat], h_mT[jetCat]->FindBin(currentEvt.MT));
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
        for(unsigned int mT = 1; mT <= h_mT_size[jetCat]; mT++){
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
          mon.setBinContentAndError("mT_final"+v_jetCat[jetCat], tagsR[lepCat].substr(1), mT, content, sqrt(error2), divideFinalHistoByBinWidth); //erase the mT final plot
        }
      }
    }
  }

  if(EvtWeights.GetSize() >= 110 && (syst_ == "pdf_up" || syst_ == "pdf_down")){ // Loop on the 100 replicas, to compute the pdf uncertainty
    for(unsigned int lepCat = 0; lepCat < tagsR.size()-1; lepCat++){
      for(unsigned int jetCat = 0; jetCat < v_jetCat.size(); jetCat++){
        for(unsigned int bin = 1 ; bin <= h_mT_size[jetCat] ; bin++){
          double pdf_mean = 0.;
          for(unsigned int rep = 0 ; rep < 100 ; rep++){
            pdf_mean+= pdfReplicas.at(jetCat).at(lepCat).at(rep)->GetBinContent(bin) / 100.;
          }
          double pdf_squaredSum = 0.;
          for(unsigned int rep = 0 ; rep < 100 ; rep++){
            pdf_squaredSum+= pow(pdfReplicas.at(jetCat).at(lepCat).at(rep)->GetBinContent(bin) - pdf_mean, 2) /99.;
          }
          double binContent = mon.getHisto("mT_final"+v_jetCat[jetCat], tagsR[lepCat].substr(1), divideFinalHistoByBinWidth)->GetBinContent(bin);
          if(syst_ == "pdf_up") binContent +=  sqrt(pdf_squaredSum);
          else if(syst_ == "pdf_down") binContent -= sqrt(pdf_squaredSum);
          mon.setBinContentAndError("mT_final"+v_jetCat[jetCat], tagsR[lepCat].substr(1), bin, binContent, 0., divideFinalHistoByBinWidth);
        }
      }
    }
  }

  if(syst_ == "QCDscale_up" || syst_ == "QCDscale_down"){ // Stewart-Tackman prescription
    for(unsigned int lepCat = 0; lepCat < tagsR.size()-1; lepCat++){
      for(unsigned int bin = 1 ; bin <= h_mT_size[eq0jets] ; bin++){
        double sigma_0 = 0., sigma_1 = 0., sigma_VBF = 0., sigma_tot = 0., sigma_0_nom = 0., sigma_1_nom = 0., sigma_VBF_nom = 0., sigma_tot_nom = 0., delta_sigma_0 = 0.;
        sigma_0 = mon.getHisto("mT_final_eq0jets", tagsR[lepCat].substr(1), divideFinalHistoByBinWidth)->GetBinContent(bin);
        sigma_1 = mon.getHisto("mT_finalBinning0j_geq1jets", tagsR[lepCat].substr(1), divideFinalHistoByBinWidth)->GetBinContent(bin);
        sigma_VBF = mon.getHisto("mT_finalBinning0j_vbf", tagsR[lepCat].substr(1), divideFinalHistoByBinWidth)->GetBinContent(bin);
        sigma_tot = sigma_0 + sigma_1 + sigma_VBF;
        sigma_0_nom = mon.getHisto("mT_finalBinning0j_eq0jets", tagsR[lepCat].substr(1)+"_nominal", divideFinalHistoByBinWidth)->GetBinContent(bin);
        sigma_1_nom = mon.getHisto("mT_finalBinning0j_geq1jets", tagsR[lepCat].substr(1)+"_nominal", divideFinalHistoByBinWidth)->GetBinContent(bin);
        sigma_VBF_nom = mon.getHisto("mT_finalBinning0j_vbf", tagsR[lepCat].substr(1)+"_nominal", divideFinalHistoByBinWidth)->GetBinContent(bin);
        sigma_tot_nom = sigma_0_nom + sigma_1_nom + sigma_VBF_nom;
        delta_sigma_0 = sqrt(pow(sigma_tot - sigma_tot_nom,2) + pow(sigma_1 - sigma_1_nom,2) + pow(sigma_VBF - sigma_VBF_nom,2));
        double binContent = 0.;
        if(sigma_0_nom == 0){
          if(syst_ == "QCDscale_up") binContent = sigma_0_nom + delta_sigma_0;
          else if(syst_ == "QCDscale_down") binContent = sigma_0_nom - delta_sigma_0;
        }
        else if(sigma_0 > sigma_0_nom) binContent = sigma_0_nom + delta_sigma_0;
        else binContent = sigma_0_nom - delta_sigma_0;
        mon.setBinContentAndError("mT_final_eq0jets", tagsR[lepCat].substr(1), bin, binContent, 0., divideFinalHistoByBinWidth);
      }
    }
  }

  mon.WriteForSysts(syst_,keepAllControlPlots_);
  outFile->Close();

}


void LooperMain::FillBTagEfficiency(
    std::vector<TLorentzVectorWithIndex> selCentralJets,
    std::vector<double> btags, TTreeReaderArray<float> const &JetAk04HadFlav,
    double weight, SmartSelectionMonitor_hzz &mon) const {

  for(unsigned int i = 0 ; i < selCentralJets.size() ; i ++){
    std::string tag = "";
    if(fabs(JetAk04HadFlav[selCentralJets.at(i).GetIndex()])==5) tag = "bjet";
    else if(fabs(JetAk04HadFlav[selCentralJets.at(i).GetIndex()])==4) tag = "cjet";
    else tag = "udsgjet";
    bool tagged_loose = btags.at(i) > 0.5426;
    bool tagged_medium = btags.at(i) > 0.8484;
    bool tagged_tight = btags.at(i) > 0.9535;
    mon.fillHisto("btagEff","den_"+tag,selCentralJets.at(i).Pt(),selCentralJets.at(i).Eta(),weight);
    if(tagged_loose) mon.fillHisto("btagEff","num_"+tag+"_tagged_loose",selCentralJets.at(i).Pt(),selCentralJets.at(i).Eta(),weight);
    if(tagged_medium) mon.fillHisto("btagEff","num_"+tag+"_tagged_medium",selCentralJets.at(i).Pt(),selCentralJets.at(i).Eta(),weight);
    if(tagged_tight) mon.fillHisto("btagEff","num_"+tag+"_tagged_tight",selCentralJets.at(i).Pt(),selCentralJets.at(i).Eta(),weight);
  }
}

