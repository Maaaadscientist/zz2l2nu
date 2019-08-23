#define InstrMETLooper_cxx

#include <BTagger.h>
#include <ElectronBuilder.h>
#include <GenJetBuilder.h>
#include <GenWeight.h>
#include <JetBuilder.h>
#include <LooperMain.h>
#include <MeKinFilter.h>
#include <MetFilters.h>
#include <MuonBuilder.h>
#include <ObjectSelection.h>
#include <PhotonBuilder.h>
#include <PhotonEfficiencySF.h>
#include <PileUpWeight.h>
#include <PtMissBuilder.h>
#include <SmartSelectionMonitor.h>
#include <SmartSelectionMonitor_hzz.h>
#include <Trigger.h>
#include <Utils.h>

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

#define MAXIMAL_AMOUNT_OF_HISTOS false //Set this to FALSE to only have the main histos. If you want all debug histo at different level, set it to TRUE

void LooperMain::Loop_InstrMET()
{
  LOG_DEBUG << "Starting the InstrMET Looper...";

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

  MeKinFilter meKinFilter{dataset_};
  MetFilters metFilters{dataset_};

  std::unique_ptr<GenWeight> genWeight;
  std::unique_ptr<PileUpWeight> pileUpWeight;
  if(isMC_) {
    genWeight.reset(new GenWeight(dataset_));
    pileUpWeight.reset(new PileUpWeight(dataset_, options_));
  }

  SmartSelectionMonitor_hzz mon;
  mon.declareHistos_InstrMET();

  std::unique_ptr<TTreeReaderArray<float>> GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass;
  std::unique_ptr<TTreeReaderArray<int>> GenPart_pdgId, GenPart_genPartIdxMother, Photon_genPartIdx;
  std::unique_ptr<TTreeReaderArray<unsigned char>> Photon_genPartFlav;
  if(isMC_) {
   GenPart_pt.reset(new TTreeReaderArray<float>(dataset_.Reader(), "GenPart_pt"));
   GenPart_eta.reset(new TTreeReaderArray<float>(dataset_.Reader(), "GenPart_eta"));
   GenPart_phi.reset(new TTreeReaderArray<float>(dataset_.Reader(), "GenPart_phi"));
   GenPart_mass.reset(new TTreeReaderArray<float>(dataset_.Reader(), "GenPart_mass"));
   GenPart_pdgId.reset(new TTreeReaderArray<int>(dataset_.Reader(), "GenPart_pdgId"));
   GenPart_genPartIdxMother.reset(new TTreeReaderArray<int>(dataset_.Reader(), "GenPart_genPartIdxMother"));
   Photon_genPartIdx.reset(new TTreeReaderArray<int>(dataset_.Reader(), "Photon_genPartIdx"));
   Photon_genPartFlav.reset(new TTreeReaderArray<unsigned char>(dataset_.Reader(), "Photon_genPartFlav"));
  }

  int64_t const nentries = dataset_.NumEntries();
  TString const fileName{dataset_.Info().Files().at(0)};

  bool isMC_QCD = (isMC_ && fileName.Contains("-QCD_"));
  bool isMC_QCD_HT = (isMC_ && fileName.Contains("-QCD_") && fileName.Contains("HT"));
  bool isMC_GJet = (isMC_ && fileName.Contains("-GJets_"));
  bool isMC_GJet_HT = (isMC_ && fileName.Contains("-GJets_HT"));
  bool isMC_WGToLNuG = (isMC_ && fileName.Contains("-WGToLNuG_") );
  bool isMC_LO_ZNuNuGJets = (isMC_ && fileName.Contains("-ZNuNuGJets_"));
  bool isMC_NLO_ZGTo2NuG_inclusive = (isMC_ && fileName.Contains("-ZGTo2NuG_") && !fileName.Contains("PtG-130"));
  bool isMC_NLO_ZGTo2NuG_Pt130 = (isMC_ && fileName.Contains("-ZGTo2NuG_PtG-130_"));
  bool isMC_ZJetsToNuNu = (isMC_ && fileName.Contains("-ZJetsToNuNu_"));

  LOG_DEBUG << "nb of entries in the input file " << fileName << " = " << nentries;

  //Compute once weights for Instr. MET reweighting if needed
  std::vector<string> v_jetCat = {"_eq0jets","_geq1jets","_vbf"};
  string const base_path = string(getenv("HZZ2L2NU_BASE")) + "/";
  bool doClosureTest = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/please_do_closure_test_when_running_InstrMETLooper");
  
  if (doClosureTest)
    LOG_INFO << "/!\\/!\\ CLOSURE TEST ONGOING - not wanted? Then remove 'WeightsAndDatadriven/InstrMET/please_do_closure_test_when_running_InstrMETLooper' /!\\/!\\" << std::endl;

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
    dataset_.SetEntry(jentry);

    if (jentry % 10000 == 0)
      LOG_INFO << Logger::TimeStamp << " Event " << jentry << " out of " <<
        nentries;

    if (not meKinFilter())
      continue;

    photon_evt currentEvt;

    double weight = 1.;
    int eventflowStep=0;

    //get the MC event weight if exists
    if (isMC_)
      weight *= (*genWeight)() * intLumi_;


    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //output of bonzais
    eventflowStep++;

    //Cleaning of low stats MC spikes that are gathering in some specific spot (in MET phi, pt, MET delta phi(MET, spike)...).
    bool isPathologicEvent=false;
    if(isMC_) isPathologicEvent = objectSelection::cleanPathologicEventsInPhotons(fileName, *run, *luminosityBlock, *event);
    if(isPathologicEvent) continue;

    // Remove events with 0 vtx
    if(*PV_npvsGood == 0 ) continue;

    //###############################################################
    //##################     OBJECT SELECTION      ##################
    //###############################################################

    auto const &tightElectrons = electronBuilder.GetTight();
    auto const &looseElectrons = electronBuilder.GetLoose();

    auto const &tightMuons = muonBuilder.GetTight();
    auto const &looseMuons = muonBuilder.GetLoose();

    auto const &photons = photonBuilder.Get();
    auto const &jets = jetBuilder.Get();

    //Ask for a prompt photon
    if(photons.size() != 1) continue;

    //Check trigger and find prescale
    int triggerWeight =0;
    /*
    int triggerType;
    if(isMC_) triggerType = trigger::MC_Photon;
    else triggerType = trigger::SinglePhoton;
    */

    //triggerWeight = trigger::passTrigger(triggerType, *TrigHltDiMu, *TrigHltMu, *TrigHltDiEl, *TrigHltEl, *TrigHltElMu, *TrigHltPhot, TrigHltDiMu_prescale, TrigHltMu_prescale, TrigHltDiEl_prescale, TrigHltEl_prescale, TrigHltElMu_prescale, TrigHltPhot_prescale, photons[0].p4.Pt());
    triggerWeight = 1.; //FIXME no prescales in NanoAOD
    if(triggerWeight==0) continue; //trigger not found

    if(MAXIMAL_AMOUNT_OF_HISTOS) mon.fillHisto("pT_Boson","noPrescale",photons[0].p4.Pt(),weight);
    weight *= triggerWeight;
    if(MAXIMAL_AMOUNT_OF_HISTOS) mon.fillHisto("pT_Boson","withPrescale",photons[0].p4.Pt(),weight);
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after prescale
    eventflowStep++;

    //photon efficiencies
    PhotonEfficiencySF phoEff;
    // if(isMC_) weight *= phoEff.getPhotonEfficiency(photons[0].p4.Pt(), photons[0].etaSc, "tight",utils::CutVersion::Moriond17Cut ).first; // FIXME broken since we don't have etaSC for photons in NanoAOD.
    if(MAXIMAL_AMOUNT_OF_HISTOS) mon.fillHisto("pT_Boson","withPrescale_and_phoEff",photons[0].p4.Pt(),weight);

    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after Photon Efficiency
    eventflowStep++;

    if(isMC_)
      weight *= (*pileUpWeight)();

    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after PU reweighting
    eventflowStep++;

    if (not metFilters())
      continue;

    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after MET filters
    eventflowStep++;


    //Resolve G+jet/QCD mixing (avoid double counting of photons)
    bool passVetoQCDevent = true;
    if(isMC_QCD){
      for(unsigned int i = 0; i < Photon_pt.GetSize(); i++){
        if(Photon_genPartFlav->At(i)==1 && GenPart_pt->At(Photon_genPartIdx->At(i)) >25){ //Gjets generated prompt photon above 25 GeV. QCD above 10 GeV, so the double counting occurs above 25.
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
      for (unsigned int i =0; i < GenPart_pt->GetSize(); i++){
        if(fabs(GenPart_pdgId->At(i))==12 || fabs(GenPart_pdgId->At(i))==14 || fabs(GenPart_pdgId->At(i))==16){
          if(fabs(GenPart_genPartIdxMother->At(i))==23 /*&& genParticle.mother()->status()==62*/){ //after testing, the status is not needed at all.
            tmpVector.SetPtEtaPhiM(GenPart_pt->At(i), GenPart_eta->At(i), GenPart_phi->At(i), GenPart_mass->At(i));
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

    if( isMC_NLO_ZGTo2NuG_inclusive && photons[0].p4.Pt() >= 130) continue;
    if( isMC_NLO_ZGTo2NuG_Pt130 && photons[0].p4.Pt() < 130) continue;

    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after LO-to-NLO k-factor for ZnunuGamma
    eventflowStep++;

    //###############################################################
    //##################       ANALYSIS CUTS       ##################
    //###############################################################


    //Definition of the relevant analysis variables and storage in the currentEvt
    TLorentzVector boson = photons[0].p4;

    //Jet category
    enum {eq0jets,geq1jets,vbf};
    int jetCat = geq1jets;

    if (jets.size() == 0)
      jetCat = eq0jets;
    else if (utils::PassVbfCuts(jets, boson))
      jetCat = vbf;

    auto const &ptMiss = ptMissBuilder.Get();
    TLorentzVector const &ptMissP4 = ptMiss.p4;

    currentEvt.Fill_photonEvt(
      v_jetCat[jetCat], tagsR[0], boson, ptMissP4, jets, *run,
      *PV_npvsGood, *fixedGridRhoFastjetAll, ptMiss.significance);

    mon.fillHisto("jetCategory","afterWeight",jetCat,weight);
    mon.fillAnalysisHistos(currentEvt, "afterWeight", weight);

    //std::cout<<"Event info: " << run<<":"<<luminosityBlock<<":"<<event << "; boson pt = "<<boson.Pt()<<"; weight = "<<weight<<"; triggerPrescale = "<<triggerWeight<<"; met = "<<currentEvt.MET<<"; mt = "<<currentEvt.MT<<"; njets = "<<currentEvt.nJets<<"; vtx = "<<PV_npvsGood<<"; rho = "<<fixedGridRhoFastjetAll<<"; puWeight = "<<weightPU<<std::endl;

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
    if (looseElectrons.size() + looseMuons.size() > 0)
      continue;

    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight); //after no extra leptons
    eventflowStep++;

    // b veto
    bool passBTag = true;

    for (auto const &jet : jets)
      if (bTagger(jet)) {
        passBTag = false;
        break;
      }

    if (not passBTag)
      continue;

    //mon.fillHisto("eventflow","tot"+tagsR[c],eventflowStep,weight); //after b-tag veto
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight);
    eventflowStep++;

    // Phi(jet,MET)
    bool passDeltaPhiJetMET = true;

    for (auto const &jet : jets)
      if (std::abs(utils::deltaPhi(jet.p4, ptMissP4)) < 0.5) {
        passDeltaPhiJetMET = false;
        break;
      }

    if (not passDeltaPhiJetMET)
      continue;

    //mon.fillHisto("eventflow","tot"+tagsR[c],eventflowStep,weight); //after delta phi (jet, met)
    for(unsigned int i = 0; i < tagsR_size; i++) mon.fillHisto("eventflow","tot"+tagsR[i],eventflowStep,weight);
    eventflowStep++;

    mon.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweighting", weight);

    // -- Histograms used to compute weights for the Instr. MET estimation : NVtx part --
    if(ptMissP4.Pt()<125){
      mon.fillHisto("reco-vtx_MET125",    "InstrMET_reweighting"+currentEvt.s_jetCat+currentEvt.s_lepCat, *PV_npvsGood, weight, true);
      mon.fillHisto("reco-vtx_MET125",    "InstrMET_reweighting"+currentEvt.s_lepCat, *PV_npvsGood, weight, true); //for all jet cats
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
        itlow = NVtxWeight_map[tagsR[c]].upper_bound(*PV_npvsGood); //look at which bin in the map currentEvt.rho corresponds
        if(itlow == NVtxWeight_map[tagsR[c]].begin()) throw std::out_of_range("You are trying to access your NVtx reweighting map outside of bin boundaries)");
        itlow--;

        weight *= itlow->second.first; //don't apply for first element of the map which is the normal one without reweighting.
      }

      mon.fillHisto("eventflow","tot"+tagsR[c],eventflowStep,weight); //after ee or mumu reweighting
      eventflowStep++;

      mon.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweightingAfter"+tagsR[c], weight);

      // -- Histograms used to compute weights for the Instr. MET estimation : Pt part --
      if(ptMissP4.Pt()<125){
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
        currentEvt.MT = sqrt(pow(sqrt(pow(boson.Pt(),2)+pow(boson.M(),2))+sqrt(pow(ptMissP4.Pt(),2)+pow(91.1876,2)),2)-pow((boson+ptMissP4).Pt(),2));
        currentEvt.M_Boson = boson.M();
      }

      mon.fillInstrMETControlRegionHisto(currentEvt, "InstrMET_AllWeightsAndLineshapeApplied"+tagsR[c], weight);
      mon.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweightingAfter"+tagsR[c]+"AfterPtR_andMassivePhoton", weight);

      mon.fillAnalysisHistos(currentEvt, "beforeMETcut_After"+tagsR[c], weight);

      //MET>80
      if(ptMissP4.Pt()<80) continue;
      mon.fillHisto("eventflow","tot"+tagsR[c],eventflowStep,weight); //after MET > 80
      eventflowStep++;

      //MET>125
      if(ptMissP4.Pt()<125) continue;
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
