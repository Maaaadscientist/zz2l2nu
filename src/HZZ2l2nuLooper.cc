#define HZZ2l2nuLooper_cxx

#include <cmath>
#include <cstdlib>
#include <limits>

#include <BTagger.h>
#include <BTagWeight.h>
#include <ElectronBuilder.h>
#include <EWCorrectionWeight.h>
#include <GenJetBuilder.h>
#include <GenWeight.h>
#include <JetBuilder.h>
#include <KFactorCorrection.h>
#include <LeptonsEfficiencySF.h>
#include <LooperMain.h>
#include <MelaWeight.h>
#include <MetFilters.h>
#include <MuonBuilder.h>
#include <ObjectSelection.h>
#include <PhotonBuilder.h>
#include <PhotonEfficiencySF.h>
#include <PileUpWeight.h>
#include <PtMissBuilder.h>
#include <SmartSelectionMonitor.h>
#include <Trigger.h>
#include <Utils.h>

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <algorithm>




void LooperMain::Loop()
{
  //Get file info
  int64_t const nentries = dataset_.NumEntries();
  TString const fileName{dataset_.Info().Files().at(0)};
  bool isMC_Wlnu_inclusive = (isMC_ && fileName.Contains("-WJetsToLNu_") && !fileName.Contains("HT"));
  bool isMC_Wlnu_HT100 = (isMC_ && fileName.Contains("-WJetsToLNu_HT-") );
  bool isMC_NLO_ZGTo2NuG_inclusive = (isMC_ && fileName.Contains("-ZGTo2NuG_") && !fileName.Contains("PtG-130"));
  bool isMC_NLO_ZGTo2NuG_Pt130 = (isMC_ && fileName.Contains("-ZGTo2NuG_PtG-130_"));

  //###############################################################
  //################## DECLARATION OF HISTOGRAMS ##################
  //###############################################################
  BTagger bTagger{options_};
  
  ElectronBuilder electronBuilder{dataset_, options_};
  MuonBuilder muonBuilder{dataset_, options_, randomGenerator_};

  PhotonBuilder photonBuilder{dataset_, options_};
  photonBuilder.EnableCleaning({&muonBuilder, &electronBuilder});

  std::unique_ptr<GenJetBuilder> genJetBuilder;
  JetBuilder jetBuilder{dataset_, options_, randomGenerator_};
  jetBuilder.EnableCleaning({&muonBuilder, &electronBuilder, &photonBuilder});
  if(isMC_) {
    genJetBuilder.reset(new GenJetBuilder(dataset_, options_));
    jetBuilder.SetGenJetBuilder(genJetBuilder.get());
  }

  PtMissBuilder ptMissBuilder{dataset_};
  ptMissBuilder.PullCalibration({&muonBuilder, &electronBuilder, &photonBuilder,
                                 &jetBuilder});

  MetFilters metFilters{dataset_};

  std::unique_ptr<GenWeight> genWeight;
  std::unique_ptr<EWCorrectionWeight> ewCorrectionWeight;
  std::unique_ptr<PileUpWeight> pileUpWeight;
  std::unique_ptr<KFactorCorrection> kfactorCorrection;
  if(isMC_) {
    genWeight.reset(new GenWeight(dataset_));
    ewCorrectionWeight.reset(new EWCorrectionWeight(dataset_, options_));
    pileUpWeight.reset(new PileUpWeight(dataset_, options_));
    kfactorCorrection.reset(new KFactorCorrection(dataset_, options_));
  }

  BTagWeight bTagWeight(options_, bTagger);
  MelaWeight melaWeight{dataset_, options_};

  SmartSelectionMonitor_hzz mon;
  mon.declareHistos();

  LOG_DEBUG << "nb of entries in the input file =" << nentries;

  LOG_DEBUG << "fileName is " << fileName;

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
    dataset_.SetEntry(jentry);

    if (jentry % 10000 == 0)
      LOG_INFO << Logger::TimeStamp << " Event " << jentry << " out of " <<
        nentries;

    if (not metFilters())
      continue;

    evt currentEvt;

    double theRandomNumber = randomGenerator_.Rndm(); //Used for the uncertainty on the 3rd lepton veto.

    double weight = 1.;
    //get the MC event weight if exists
    if (isMC_) {
      weight *= (*genWeight)() * intLumi_;

      //get the PU weights
      weight *= (*pileUpWeight)();
    }

    mon.fillHisto("eventflow","tot",0,weight);

    // Remove events with 0 vtx
    if(*PV_npvsGood == 0 ) continue;

    for(int i =0 ; i < Muon_pt.GetSize() ; i++) mon.fillHisto("pT_mu","tot",Muon_pt[i],weight);
    for(int i =0 ; i < Electron_pt.GetSize() ; i++) mon.fillHisto("pT_e","tot",Electron_pt[i],weight);
    mon.fillHisto("nb_mu","tot",Muon_pt.GetSize(),weight);
    mon.fillHisto("nb_e","tot",Electron_pt.GetSize(),weight);
    mon.fillHisto("reco-vtx","tot",*PV_npvsGood,weight);


    //###############################################################
    //##################     OBJECT CORRECTIONS    ##################
    //###############################################################
     
    // MELA weight and kfactor
    weight *= melaWeight();
    if(isMC_) weight *= (*kfactorCorrection)();

    // electroweak corrections
    if(isMC_) weight *= (*ewCorrectionWeight)();

    // Theory uncertainties
    double thUncWeight = 1.;
    
    if (syst_ != "")
      thUncWeight = utils::getTheoryUncertainties(*genWeight, syst_);
    
    if(thUncWeight == 0) continue; // There are some rare cases where a weight is at 0, making an indeterminate form (0/0) in the code. I figured out it was an easy (albeit a bit coward) way to avoid it without changing all the code for an effect of less than 0.01%.
    weight *= thUncWeight;

    LOG_TRACE << "Weight after corrections: " << weight;

    //###############################################################
    //##################     OBJECT SELECTION      ##################
    //###############################################################

    auto const &tightElectrons = electronBuilder.GetTight();
    auto const &looseElectrons = electronBuilder.GetLoose();

    auto const &tightMuons = muonBuilder.GetTight();
    auto const &looseMuons = muonBuilder.GetLoose();

    auto const &photons = photonBuilder.Get();
    auto const &jets = jetBuilder.Get();

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
        //PhotonEfficiencySF phoEff;
        //weight *= phoEff.getPhotonEfficiency(
        //  photons[0].p4.Pt(), photons[0].etaSc, "tight",
        //  utils::CutVersion::Moriond17Cut).first;
        // FIXME Broken since we don't have etaSC for photons. This will need to be fixed.
      }
    }

    //trigger weights for photon data
    if(isPhotonDatadriven_){
      int triggerWeight = 0, triggerType = 0;
      /*
      if(isMC_) triggerType = trigger::MC_Photon;
      else triggerType = trigger::SinglePhoton;
      */

      //triggerWeight = trigger::passTrigger(triggerType, *TrigHltDiMu, *TrigHltMu, *TrigHltDiEl, *TrigHltEl, *TrigHltElMu, *TrigHltPhot, TrigHltDiMu_prescale, TrigHltMu_prescale, TrigHltDiEl_prescale, TrigHltEl_prescale, TrigHltElMu_prescale, TrigHltPhot_prescale, photons[0].p4.Pt());
      triggerWeight = 1.; //FIXME no prescales in NanoAOD
      if(triggerWeight==0) continue; //trigger not found
      weight *= triggerWeight;
    }

    //Avoid double counting for W+jets
    //For some reasons we just have the inclusive sample for the Dilepton region while we have both HT and inclusive samples for the photon region. Hence this cleaning only applies to the photon region.
    if(isPhotonDatadriven_){
      if (isMC_Wlnu_inclusive || isMC_Wlnu_HT100){ //Avoid double counting and make our W#rightarrow l#nu exclusif of the dataset with a cut on HT...
        bool isHT100 = false;

        //Let's create our own gen HT variable
        double vHT = 0;

        for (auto const &genJet : genJetBuilder->Get())
          if (not muonBuilder.GetMomenta().HasOverlap(genJet.p4, 0.4) and
              not electronBuilder.GetMomenta().HasOverlap(genJet.p4, 0.4))
            vHT += genJet.p4.Pt();

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

    int jetCat = geq1jets;

    if (jets.size() == 0)
      jetCat = eq0jets;
    else if (utils::PassVbfCuts(jets, boson))
      jetCat = vbf;
    
    currentEvt.s_jetCat = v_jetCat[jetCat];

    TLorentzVector const ptMissP4 = ptMissBuilder.Get().p4;

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
        itlow = NVtxWeight_map[tagsR[c]].upper_bound(*PV_npvsGood); //look at which bin in the map currentEvt.nVtx corresponds
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
        v_jetCat[jetCat], tagsR[c], boson, ptMissP4, jets, *run,
        *PV_npvsGood, *fixedGridRhoFastjetAll, /**MET_significance, */tightLeptons);

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
        FillBTagEfficiency(jets, weight, mon, bTagger);

      // b veto
      bool passBTag = true;

      for (auto const &jet : jets)
        if (bTagger(jet)) {
          passBTag = false;
          break;
        }

      if (not passBTag)
        continue;

      // Apply the btag weights
      if (isMC_) {
        double const w = bTagWeight(jets);
        weight *= w;

        if (currentEvt.s_lepCat == "_ll")
          mon.fillProfile("BTagWeightvsMT", "TEST", currentEvt.MT, w, weight); //FIXME remove after
      }

      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",6,weight);

      // Phi(jet,MET)
      bool passDeltaPhiJetMET = true;

      for (auto const &jet : jets)
        if (std::abs(utils::deltaPhi(jet.p4, ptMissP4)) < 0.5) {
          passDeltaPhiJetMET = false;
          break;
        }

      if (not passDeltaPhiJetMET)
        continue;

      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",7,weight);

      mon.fillInstrMETControlRegionHisto(currentEvt, "InstrMET_reweighting", weight);
      mon.fillAnalysisHistos(currentEvt, "beforeMETcut", weight);
      mon.fillHisto("jetCategory","beforeMETcut"+currentEvt.s_lepCat,jetCat,weight);

      //MET>80
      if(ptMissP4.Pt()<80) continue;
      if(currentEvt.s_lepCat == "_ll") mon.fillHisto("eventflow","tot",8,weight);

      //MET>125
      if(ptMissP4.Pt()<125) continue;
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

      if ((syst_ == "pdf_up" || syst_ == "pdf_down")
          && currentEvt.s_lepCat != "_ll") {
        for (int i = 0 ; i < 100; i++)
          pdfReplicas.at(jetCat).at(lepCat).at(i)->Fill(
            currentEvt.MT, weight * genWeight->RelWeightPdf(i));
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

  if(syst_ == "pdf_up" || syst_ == "pdf_down"){ // Loop on the 100 replicas, to compute the pdf uncertainty
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
    std::vector<Jet> const &jets, double weight,
    SmartSelectionMonitor_hzz &mon, BTagger const &bTagger) const {

  for (auto const &jet : jets) {
    if ( not bTagger.IsTaggable(jet))
      continue;

    std::string tag;

    if (std::abs(jet.hadronFlavour) == 5)
      tag = "bjet";
    else if (std::abs(jet.hadronFlavour) == 4)
      tag = "cjet";
    else
      tag = "udsgjet";

    auto fillHistogram = [&](std::string const &label) {
      mon.fillHisto("btagEff", label, jet.p4.Pt(), jet.p4.Eta(), weight);
    };

    fillHistogram("den_" + tag);

    if (bTagger(jet))
      fillHistogram("num_" + tag);
  }
}

