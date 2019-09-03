#include <InstrMetAnalysis.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <TFile.h>

#include <ObjectSelection.h>
#include <PhotonEfficiencySF.h>
#include <Utils.h>

#define MAXIMAL_AMOUNT_OF_HISTOS false //Set this to FALSE to only have the main histos. If you want all debug histo at different level, set it to TRUE


namespace po = boost::program_options;


InstrMetAnalysis::InstrMetAnalysis(Options const &options, Dataset &dataset)
    : dataset_{dataset}, isMC_{dataset_.Info().IsSimulation()},
      intLumi_{options.GetConfig()["luminosity"].as<double>()},
      outputFile_{options.GetAs<std::string>("output")},
      syst_{options.GetAs<std::string>("syst")},
      randomGenerator_{options.GetAs<unsigned>("seed")},
      bTagger_{options},
      electronBuilder_{dataset_, options},
      muonBuilder_{dataset_, options, randomGenerator_},
      photonBuilder_{dataset_, options},
      jetBuilder_{dataset_, options, randomGenerator_},
      ptMissBuilder_{dataset_},
      meKinFilter_{dataset_}, metFilters_{dataset_},
      bTagWeight_{options, bTagger_},
      divideFinalHistoByBinWidth_{false},  //For final plots, we don't divide by the bin width to ease computations of the yields by eye.
      v_jetCat_{"_eq0jets","_geq1jets","_vbf"},
      fileName_{dataset_.Info().Files().at(0)}
{
  if (isMC_) {
    genJetBuilder_.reset(new GenJetBuilder(dataset_, options));
    jetBuilder_.SetGenJetBuilder(genJetBuilder_.get());
  }
  
  // Cross-cleaning
  photonBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});
  jetBuilder_.EnableCleaning(
      {&muonBuilder_, &electronBuilder_, &photonBuilder_});
  
  // Type 1 corrections for ptmiss
  ptMissBuilder_.PullCalibration(
      {&muonBuilder_, &electronBuilder_, &photonBuilder_, &jetBuilder_});

  if (isMC_) {
    genWeight_.reset(new GenWeight{dataset_});
    pileUpWeight_.reset(new PileUpWeight{dataset_, options});
  }

  isMC_QCD_ = (isMC_ && fileName_.Contains("-QCD_"));
  isMC_GJet_HT_ = (isMC_ && fileName_.Contains("-GJets_HT"));
  isMC_LO_ZNuNuGJets_ = (isMC_ && fileName_.Contains("-ZNuNuGJets_"));
  isMC_NLO_ZGTo2NuG_inclusive_ = (isMC_ && fileName_.Contains("-ZGTo2NuG_") && !fileName_.Contains("PtG-130"));
  isMC_NLO_ZGTo2NuG_Pt130_ = (isMC_ && fileName_.Contains("-ZGTo2NuG_PtG-130_"));

  if (isMC_) {
    genPartPt_.reset(new TTreeReaderArray<float>(dataset_.Reader(), "GenPart_pt"));
    genPartEta_.reset(new TTreeReaderArray<float>(dataset_.Reader(), "GenPart_eta"));
    genPartPhi_.reset(new TTreeReaderArray<float>(dataset_.Reader(), "GenPart_phi"));
    genPartMass_.reset(new TTreeReaderArray<float>(dataset_.Reader(), "GenPart_mass"));
    genPartPdgId_.reset(new TTreeReaderArray<int>(dataset_.Reader(), "GenPart_pdgId"));
    genPartMotherIndex_.reset(new TTreeReaderArray<int>(dataset_.Reader(), "GenPart_genPartIdxMother"));
    photonGenPartIndex_.reset(new TTreeReaderArray<int>(dataset_.Reader(), "Photon_genPartIdx"));
    photonGenPartFlavor_.reset(new TTreeReaderArray<unsigned char>(dataset_.Reader(), "Photon_genPartFlav"));
  }

  mon_.declareHistos_InstrMET();

  //Compute once weights for Instr. MET reweighting if needed
  auto const base_path = std::string(std::getenv("HZZ2L2NU_BASE")) + "/";
  bool doClosureTest = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/please_do_closure_test_when_running_InstrMETLooper");
  
  if (doClosureTest)
    LOG_INFO << "/!\\/!\\ CLOSURE TEST ONGOING - not wanted? Then remove 'WeightsAndDatadriven/InstrMET/please_do_closure_test_when_running_InstrMETLooper' /!\\/!\\" << std::endl;

  std::string weightFileType = (doClosureTest) ? "closureTest" : "InstrMET";
  weight_NVtx_exist_ = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_NVtx.root");
  weight_Pt_exist_ = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_weight_pt.root");
  weight_Mass_exist_ = utils::file_exist(base_path+"WeightsAndDatadriven/InstrMET/"+weightFileType+"_lineshape_mass.root");
  utils::loadInstrMETWeights(weight_NVtx_exist_, weight_Pt_exist_, weight_Mass_exist_, nVtxWeight_map_, ptWeight_map_, lineshapeMassWeight_map_, weightFileType, base_path, v_jetCat_);

  tagsR_.push_back("_gamma"); //_gamma, i.e. no reweighting to ee or mumu
  if (weight_NVtx_exist_) {
    tagsR_.push_back("_ee");
    tagsR_.push_back("_mumu");
    tagsR_.push_back("_ll");
  }
  tagsR_size_ = tagsR_.size();

  if (syst_ == "")
    LOG_DEBUG << "Will not apply systematic variations.";
  else
    LOG_DEBUG << "Will apply systematic variation \"" << syst_ << "\".";
}


po::options_description InstrMetAnalysis::OptionsDescription() {
  po::options_description optionsDescription{"Analysis-specific options"};
  optionsDescription.add_options()
    ("syst", po::value<std::string>()->default_value(""),
     "Requested systematic variation")
    ("output,o", po::value<std::string>()->default_value("outputFile.root"),
     "Name for output file with histograms")
    ("seed", po::value<unsigned>()->default_value(0),
     "Seed for random number generator; 0 means a unique seed");
  return optionsDescription;
}


void InstrMetAnalysis::PostProcessing() {
  TFile *outFile = TFile::Open(outputFile_.c_str(), "recreate");
  mon_.Write();
  outFile->Close();
}


bool InstrMetAnalysis::ProcessEvent() {
  if (not meKinFilter_())
    return false;

  photon_evt currentEvt;

  double weight = 1.;
  int eventflowStep=0;

  //get the MC event weight if exists
  if (isMC_)
    weight *= (*genWeight_)() * intLumi_;


  for(unsigned int i = 0; i < tagsR_size_; i++) mon_.fillHisto("eventflow","tot"+tagsR_[i],eventflowStep,weight); //output of bonzais
  eventflowStep++;

  //Cleaning of low stats MC spikes that are gathering in some specific spot (in MET phi, pt, MET delta phi(MET, spike)...).
  bool isPathologicEvent=false;
  if(isMC_) isPathologicEvent = objectSelection::cleanPathologicEventsInPhotons(fileName_, *run_, *luminosityBlock_, *eventNumber_);
  if(isPathologicEvent)
    return false;

  // Remove events with 0 vtx
  if(*numPVGood_ == 0 )
    return false;

  //###############################################################
  //##################     OBJECT SELECTION      ##################
  //###############################################################

  auto const &tightElectrons = electronBuilder_.GetTight();
  auto const &looseElectrons = electronBuilder_.GetLoose();

  auto const &tightMuons = muonBuilder_.GetTight();
  auto const &looseMuons = muonBuilder_.GetLoose();

  auto const &photons = photonBuilder_.Get();
  auto const &jets = jetBuilder_.Get();

  //Ask for a prompt photon
  if(photons.size() != 1)
    return false;

  //Check trigger and find prescale
  int triggerWeight =0;
  /*
  int triggerType;
  if(isMC_) triggerType = trigger::MC_Photon;
  else triggerType = trigger::SinglePhoton;
  */

  //triggerWeight = trigger::passTrigger(triggerType, *TrigHltDiMu, *TrigHltMu, *TrigHltDiEl, *TrigHltEl, *TrigHltElMu, *TrigHltPhot, TrigHltDiMu_prescale, TrigHltMu_prescale, TrigHltDiEl_prescale, TrigHltEl_prescale, TrigHltElMu_prescale, TrigHltPhot_prescale, photons[0].p4.Pt());
  triggerWeight = 1.; //FIXME no prescales in NanoAOD
  if(triggerWeight==0)  //trigger not found
    return false;

  if(MAXIMAL_AMOUNT_OF_HISTOS) mon_.fillHisto("pT_Boson","noPrescale",photons[0].p4.Pt(),weight);
  weight *= triggerWeight;
  if(MAXIMAL_AMOUNT_OF_HISTOS) mon_.fillHisto("pT_Boson","withPrescale",photons[0].p4.Pt(),weight);
  for(unsigned int i = 0; i < tagsR_size_; i++) mon_.fillHisto("eventflow","tot"+tagsR_[i],eventflowStep,weight); //after prescale
  eventflowStep++;

  //photon efficiencies
  PhotonEfficiencySF phoEff;
  // if(isMC_) weight *= phoEff.getPhotonEfficiency(photons[0].p4.Pt(), photons[0].etaSc, "tight",utils::CutVersion::Moriond17Cut ).first; // FIXME broken since we don't have etaSC for photons in NanoAOD.
  if(MAXIMAL_AMOUNT_OF_HISTOS) mon_.fillHisto("pT_Boson","withPrescale_and_phoEff",photons[0].p4.Pt(),weight);

  for(unsigned int i = 0; i < tagsR_size_; i++) mon_.fillHisto("eventflow","tot"+tagsR_[i],eventflowStep,weight); //after Photon Efficiency
  eventflowStep++;

  if(isMC_)
    weight *= (*pileUpWeight_)();

  for(unsigned int i = 0; i < tagsR_size_; i++) mon_.fillHisto("eventflow","tot"+tagsR_[i],eventflowStep,weight); //after PU reweighting
  eventflowStep++;

  if (not metFilters_())
    return false;

  for(unsigned int i = 0; i < tagsR_size_; i++) mon_.fillHisto("eventflow","tot"+tagsR_[i],eventflowStep,weight); //after MET filters
  eventflowStep++;


  //Resolve G+jet/QCD mixing (avoid double counting of photons)
  bool passVetoQCDevent = true;
  if(isMC_QCD_){
    for(unsigned int i = 0; i < photonGenPartFlavor_->GetSize(); i++){
      if(photonGenPartFlavor_->At(i)==1 && genPartPt_->At(photonGenPartIndex_->At(i)) >25){ //Gjets generated prompt photon above 25 GeV. QCD above 10 GeV, so the double counting occurs above 25.
        passVetoQCDevent = false;
      }
    }
  }
  if(!passVetoQCDevent)
    return false;

  for(unsigned int i = 0; i < tagsR_size_; i++) mon_.fillHisto("eventflow","tot"+tagsR_[i],eventflowStep,weight); //after avoiding G+jets and QCD mixing
  eventflowStep++;

  //LO to NLO k-factor for ZNuNuGamma (ref: fig 16 (bottom right) of http://link.springer.com/article/10.1007%2FJHEP02%282016%29057)
  double kFactor_ZNuNuGWeight = 1.;
  if(isMC_LO_ZNuNuGJets_){
    //reconstruct the gen transverse energy
    std::vector<TLorentzVector> genNeutrinosFromZ;
    TLorentzVector tmpVector;
    for (unsigned int i =0; i < genPartPt_->GetSize(); i++){
      if(fabs(genPartPdgId_->At(i))==12 || fabs(genPartPdgId_->At(i))==14 || fabs(genPartPdgId_->At(i))==16){
        if(fabs(genPartMotherIndex_->At(i))==23 /*&& genParticle.mother()->status()==62*/){ //after testing, the status is not needed at all.
          tmpVector.SetPtEtaPhiM(genPartPt_->At(i), genPartEta_->At(i), genPartPhi_->At(i), genPartMass_->At(i));
          genNeutrinosFromZ.push_back(tmpVector);//neutrino originating directly from Z boson
        }
      }
    }
    if(genNeutrinosFromZ.size() < 2)
      return false;
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

  if( isMC_NLO_ZGTo2NuG_inclusive_ && photons[0].p4.Pt() >= 130)
    return false;
  if( isMC_NLO_ZGTo2NuG_Pt130_ && photons[0].p4.Pt() < 130)
    return false;

  for(unsigned int i = 0; i < tagsR_size_; i++) mon_.fillHisto("eventflow","tot"+tagsR_[i],eventflowStep,weight); //after LO-to-NLO k-factor for ZnunuGamma
  eventflowStep++;

  //###############################################################
  //##################       ANALYSIS CUTS       ##################
  //###############################################################


  //Definition of the relevant analysis variables and storage in the currentEvt
  TLorentzVector boson = photons[0].p4;

  //Jet category
  int jetCat = geq1jets;

  if (jets.size() == 0)
    jetCat = eq0jets;
  else if (utils::PassVbfCuts(jets, boson))
    jetCat = vbf;

  auto const &ptMiss = ptMissBuilder_.Get();
  TLorentzVector const &ptMissP4 = ptMiss.p4;

  currentEvt.Fill_photonEvt(
    v_jetCat_[jetCat], tagsR_[0], boson, ptMissP4, jets, *run_,
    *numPVGood_, *rho_, ptMiss.significance);

  mon_.fillHisto("jetCategory","afterWeight",jetCat,weight);
  mon_.fillAnalysisHistos(currentEvt, "afterWeight", weight);

  //std::cout<<"Event info: " << run_<<":"<<luminosityBlock_<<":"<<eventNumber_ << "; boson pt = "<<boson.Pt()<<"; weight = "<<weight<<"; triggerPrescale = "<<triggerWeight<<"; met = "<<currentEvt.MET<<"; mt = "<<currentEvt.MT<<"; njets = "<<currentEvt.nJets<<"; vtx = "<<numPVGood_<<"; rho = "<<fixedGridRhoFastjetAll<<"; puWeight = "<<weightPU<<std::endl;

  if(boson.Pt() < 55.)
    return false;
  for(unsigned int i = 0; i < tagsR_size_; i++) mon_.fillHisto("eventflow","tot"+tagsR_[i],eventflowStep,weight); //after pt cut
  eventflowStep++;

  //LO-to-NLO k-factor for GJets_HT
  if(isMC_GJet_HT_) weight *= std::max(1., 1.716910-0.001221*boson.Pt()); //We use the same weights than the ones used by JME-17-001. However when the weight becomes lower than one (at 587.15 GeV) we keep the weight =1. This looks like the weights we found when comparing our LO samples to our NLO samples.

  //Phi(Z,MET)
  if(currentEvt.deltaPhi_MET_Boson<0.5)
    return false;
  for(unsigned int i = 0; i < tagsR_size_; i++) mon_.fillHisto("eventflow","tot"+tagsR_[i],eventflowStep,weight); //after delta phi (Z, met)
  eventflowStep++;

  //No Extra Lepton
  if (looseElectrons.size() + looseMuons.size() > 0)
    return false;

  for(unsigned int i = 0; i < tagsR_size_; i++) mon_.fillHisto("eventflow","tot"+tagsR_[i],eventflowStep,weight); //after no extra leptons
  eventflowStep++;

  // b veto
  bool passBTag = true;

  for (auto const &jet : jets)
    if (bTagger_(jet)) {
      passBTag = false;
      break;
    }

  if (not passBTag)
    return false;

  //mon_.fillHisto("eventflow","tot"+tagsR_[c],eventflowStep,weight); //after b-tag veto
  for(unsigned int i = 0; i < tagsR_size_; i++) mon_.fillHisto("eventflow","tot"+tagsR_[i],eventflowStep,weight);
  eventflowStep++;

  // Phi(jet,MET)
  bool passDeltaPhiJetMET = true;

  for (auto const &jet : jets)
    if (std::abs(utils::deltaPhi(jet.p4, ptMissP4)) < 0.5) {
      passDeltaPhiJetMET = false;
      break;
    }

  if (not passDeltaPhiJetMET)
    return false;

  //mon_.fillHisto("eventflow","tot"+tagsR_[c],eventflowStep,weight); //after delta phi (jet, met)
  for(unsigned int i = 0; i < tagsR_size_; i++) mon_.fillHisto("eventflow","tot"+tagsR_[i],eventflowStep,weight);
  eventflowStep++;

  mon_.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweighting", weight);

  // -- Histograms used to compute weights for the Instr. MET estimation : NVtx part --
  if(ptMissP4.Pt()<125){
    mon_.fillHisto("reco-vtx_MET125",    "InstrMET_reweighting"+currentEvt.s_jetCat+currentEvt.s_lepCat, *numPVGood_, weight, true);
    mon_.fillHisto("reco-vtx_MET125",    "InstrMET_reweighting"+currentEvt.s_lepCat, *numPVGood_, weight, true); //for all jet cats
  }
    
  //Apply NVtx reweighting if file exist!
  //Starting from here, plots won't be "gamma" anymore but "eeR" or "mumuR". R for Reweighted.   
  double weightBeforeLoop = weight;
  double MTBeforeLoop = currentEvt.MT;
  double MBeforeLoop = currentEvt.M_Boson;
  TLorentzVector bosonBeforeLoop = boson;
  bool eventAccepted = false;

  for(unsigned int c = 0; c < tagsR_size_; c++){
    weight = weightBeforeLoop;
    currentEvt.MT = MTBeforeLoop;
    currentEvt.M_Boson = MBeforeLoop;
    boson = bosonBeforeLoop;


    if(c > 0){ //c=0 corresponds to no reweighting
      std::map<double, std::pair<double,double> >::iterator itlow;
      itlow = nVtxWeight_map_[tagsR_[c]].upper_bound(*numPVGood_); //look at which bin in the map currentEvt.rho corresponds
      if(itlow == nVtxWeight_map_[tagsR_[c]].begin()) throw std::out_of_range("You are trying to access your NVtx reweighting map outside of bin boundaries)");
      itlow--;

      weight *= itlow->second.first; //don't apply for first element of the map which is the normal one without reweighting.
    }

    mon_.fillHisto("eventflow","tot"+tagsR_[c],eventflowStep,weight); //after ee or mumu reweighting
    eventflowStep++;

    mon_.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweightingAfter"+tagsR_[c], weight);

    // -- Histograms used to compute weights for the Instr. MET estimation : Pt part --
    if(ptMissP4.Pt()<125){
      mon_.fillHisto("pT_Boson_MET125", "InstrMET_reweightingAfter"+tagsR_[c]+currentEvt.s_jetCat, boson.Pt(), weight, true); // all jet cats
      mon_.fillHisto("pT_Boson_MET125", "InstrMET_reweightingAfter"+tagsR_[c], boson.Pt(), weight, true); // all jet cats
    }

    //Apply pt weight on top of NVtxWeight... so if i>0:
    if(c > 0 && weight_Pt_exist_){
      std::map<double, std::pair<double,double> >::iterator itlow;
      itlow = ptWeight_map_[tagsR_[c]+currentEvt.s_jetCat].upper_bound(currentEvt.pT_Boson); //look at which bin in the map currentEvt.pT corresponds
      if(itlow == ptWeight_map_[tagsR_[c]+currentEvt.s_jetCat].begin()) throw std::out_of_range("You are trying to access your Pt reweighting map outside of bin boundaries)");
      itlow--;
      weight *= itlow->second.first; //don't apply for first element of the map which is the normal one without reweighting.
    }

    mon_.fillHisto("eventflow","tot"+tagsR_[c],eventflowStep,weight); //after Pt reweighting
    eventflowStep++;

    mon_.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweightingAfter"+tagsR_[c]+"AfterPtR", weight);

    //Apply mass on the photon:
    if(c > 0 && weight_Mass_exist_){
      utils::giveMassToPhoton(boson, lineshapeMassWeight_map_[tagsR_[c]]);
      currentEvt.MT = sqrt(pow(sqrt(pow(boson.Pt(),2)+pow(boson.M(),2))+sqrt(pow(ptMissP4.Pt(),2)+pow(91.1876,2)),2)-pow((boson+ptMissP4).Pt(),2));
      currentEvt.M_Boson = boson.M();
    }

    mon_.fillInstrMETControlRegionHisto(currentEvt, "InstrMET_AllWeightsAndLineshapeApplied"+tagsR_[c], weight);
    mon_.fillPhotonIDHistos_InstrMET(currentEvt, "ReadyForReweightingAfter"+tagsR_[c]+"AfterPtR_andMassivePhoton", weight);

    mon_.fillAnalysisHistos(currentEvt, "beforeMETcut_After"+tagsR_[c], weight);

    //MET>80
    if(ptMissP4.Pt()<80)
      continue;
    mon_.fillHisto("eventflow","tot"+tagsR_[c],eventflowStep,weight); //after MET > 80
    eventflowStep++;

    //MET>125
    if(ptMissP4.Pt()<125)
      continue;
    mon_.fillHisto("eventflow","tot"+tagsR_[c],eventflowStep,weight); //after MET > 125
    eventflowStep++;

    //###############################################################
    //##################     END OF SELECTION      ##################
    //###############################################################
    eventAccepted = true;
    mon_.fillAnalysisHistos(currentEvt, "final"+tagsR_[c], weight, divideFinalHistoByBinWidth_);
    //THE RESULT:
    mon_.fillHisto("mT_final"+currentEvt.s_jetCat, currentEvt.s_lepCat, currentEvt.MT, weight, divideFinalHistoByBinWidth_);

  }

  return eventAccepted;
}

