#include <FileInPath.h>
#include <MainAnalysis.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <TFile.h>

#include <HZZException.h>
#include <PhotonEfficiencySF.h>
#include <Trigger.h>
#include <Utils.h>


namespace po = boost::program_options;


MainAnalysis::MainAnalysis(Options const &options, Dataset &dataset)
    : AnalysisCommon{options, dataset},
      dataset_{dataset},
      outputFile_{options.GetAs<std::string>("output")},
      keepAllControlPlots_{options.Exists("all-control-plots")},
      isPhotonDatadriven_{options.Exists("dd-photon")},
      syst_{options.GetAs<std::string>("syst")},
      tabulatedRng_{tabulatedRngEngine_},
      photonBuilder_{dataset_, options},
      photonPrescales_{dataset, options},
      photonWeight_{dataset, options, &photonBuilder_},
      melaWeight_{dataset_, options},
      mon_{photonPrescales_.GetThresholdsBinning()},
      divideFinalHistoByBinWidth_{false},  //For final plots, we don't divide by the bin width to ease computations of the yields by eye.
      v_jetCat_{"_eq0jets","_geq1jets","_vbf"},
      tagsR_{"_ee", "_mumu", "_ll"}, tagsR_size_{unsigned(tagsR_.size())},
      h_mT_{jetCat_size}, h_mT_size_{jetCat_size}
{
  // Cross-cleaning for photons
  photonBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});

  // Jet and ptmiss builders from AnalysisCommon are not aware of the photon
  // builder. Register it.
  jetBuilder_.EnableCleaning({&photonBuilder_});
  ptMissBuilder_.PullCalibration({&photonBuilder_});

  TString const fileName{dataset_.Info().Files().at(0)};
  isMC_NLO_ZGTo2NuG_inclusive_ = (isSim_ && fileName.Contains("-ZGTo2NuG_") && !fileName.Contains("PtG-130"));
  isMC_NLO_ZGTo2NuG_Pt130_ = (isSim_ && fileName.Contains("-ZGTo2NuG_PtG-130_"));

  InitializeHistograms(options);

  if (syst_ == "")
    LOG_DEBUG << "Will not apply systematic variations.";
  else
    LOG_DEBUG << "Will apply systematic variation \"" << syst_ << "\".";
}


po::options_description MainAnalysis::OptionsDescription() {
  auto optionsDescription = AnalysisCommon::OptionsDescription();
  optionsDescription.add_options()
    ("dd-photon", "Use data-driven photon+jets background")
    ("all-control-plots", "Keep all control plots")
    ("mela-weight", po::value<unsigned>(), "MELA reweighting index");
  return optionsDescription;
}


void MainAnalysis::PostProcessing() {
  TFile *outFile = TFile::Open(outputFile_.c_str(), "recreate");

  //Construct mT plots with correct stat uncertainty for the Instr.MET part.
  if(isPhotonDatadriven_){
    double content =0, error2 = 0, N_events = 0, N_error2 = 0, reweighting = 0, reweighting_error = 0;
    for (int lepCat = 0; lepCat < int(tagsR_.size()); lepCat++) {
      for (int jetCat = 0; jetCat < int(v_jetCat_.size()); jetCat++) {
        for (int mT = 1; mT <= h_mT_size_[jetCat]; mT++) {
          content = 0;
          error2 = 0;
          for (int Vtx = 1; Vtx <= h_Vtx_->GetNbinsX(); Vtx++) {
            for (int pT_threshold = 1;
                pT_threshold <= h_pT_thresholds_->GetNbinsX(); pT_threshold++) {
              for (int pT = 1; pT <= h_pT_->GetNbinsX(); pT++) {
                N_events = (*mT_InstrMET_map_)[lepCat][jetCat][mT][Vtx][pT_threshold][pT].first; //sum of weights
                N_error2 = (*mT_InstrMET_map_)[lepCat][jetCat][mT][Vtx][pT_threshold][pT].second; //sum of weights*weights
                reweighting = (*photon_reweighting_)[lepCat][jetCat][Vtx][pT_threshold][pT].first;
                reweighting_error = (*photon_reweighting_)[lepCat][jetCat][Vtx][pT_threshold][pT].second;
                content += N_events*reweighting;
                error2 += reweighting*reweighting*N_error2 + N_events*N_events*reweighting_error*reweighting_error;
              }
            }
          }
          mon_.setBinContentAndError("mT_final"+v_jetCat_[jetCat], tagsR_[lepCat].substr(1), mT, content, sqrt(error2), divideFinalHistoByBinWidth_); //erase the mT final plot
        }
      }
    }
  }

  if(syst_ == "pdf_up" || syst_ == "pdf_down"){ // Loop on the 100 replicas, to compute the pdf uncertainty
    for (int lepCat = 0; lepCat < int(tagsR_.size()) - 1; lepCat++) {
      for (int jetCat = 0; jetCat < int(v_jetCat_.size()); jetCat++) {
        for (int bin = 1 ; bin <= h_mT_size_[jetCat] ; bin++) {
          double pdf_mean = 0.;
          for(unsigned int rep = 0 ; rep < 100 ; rep++){
            pdf_mean+= pdfReplicas_.at(jetCat).at(lepCat).at(rep)->GetBinContent(bin) / 100.;
          }
          double pdf_squaredSum = 0.;
          for(unsigned int rep = 0 ; rep < 100 ; rep++){
            pdf_squaredSum+= pow(pdfReplicas_.at(jetCat).at(lepCat).at(rep)->GetBinContent(bin) - pdf_mean, 2) /99.;
          }
          double binContent = mon_.getHisto("mT_final"+v_jetCat_[jetCat], tagsR_[lepCat].substr(1), divideFinalHistoByBinWidth_)->GetBinContent(bin);
          if(syst_ == "pdf_up") binContent +=  sqrt(pdf_squaredSum);
          else if(syst_ == "pdf_down") binContent -= sqrt(pdf_squaredSum);
          mon_.setBinContentAndError("mT_final"+v_jetCat_[jetCat], tagsR_[lepCat].substr(1), bin, binContent, 0., divideFinalHistoByBinWidth_);
        }
      }
    }
  }

  if(syst_ == "QCDscale_up" || syst_ == "QCDscale_down"){ // Stewart-Tackman prescription
    for (int lepCat = 0; lepCat < int(tagsR_.size()) - 1; lepCat++) {
      for (int bin = 1 ; bin <= h_mT_size_[eq0jets] ; bin++) {
        double sigma_0 = 0., sigma_1 = 0., sigma_VBF = 0., sigma_tot = 0., sigma_0_nom = 0., sigma_1_nom = 0., sigma_VBF_nom = 0., sigma_tot_nom = 0., delta_sigma_0 = 0.;
        sigma_0 = mon_.getHisto("mT_final_eq0jets", tagsR_[lepCat].substr(1), divideFinalHistoByBinWidth_)->GetBinContent(bin);
        sigma_1 = mon_.getHisto("mT_finalBinning0j_geq1jets", tagsR_[lepCat].substr(1), divideFinalHistoByBinWidth_)->GetBinContent(bin);
        sigma_VBF = mon_.getHisto("mT_finalBinning0j_vbf", tagsR_[lepCat].substr(1), divideFinalHistoByBinWidth_)->GetBinContent(bin);
        sigma_tot = sigma_0 + sigma_1 + sigma_VBF;
        sigma_0_nom = mon_.getHisto("mT_finalBinning0j_eq0jets", tagsR_[lepCat].substr(1)+"_nominal", divideFinalHistoByBinWidth_)->GetBinContent(bin);
        sigma_1_nom = mon_.getHisto("mT_finalBinning0j_geq1jets", tagsR_[lepCat].substr(1)+"_nominal", divideFinalHistoByBinWidth_)->GetBinContent(bin);
        sigma_VBF_nom = mon_.getHisto("mT_finalBinning0j_vbf", tagsR_[lepCat].substr(1)+"_nominal", divideFinalHistoByBinWidth_)->GetBinContent(bin);
        sigma_tot_nom = sigma_0_nom + sigma_1_nom + sigma_VBF_nom;
        delta_sigma_0 = sqrt(pow(sigma_tot - sigma_tot_nom,2) + pow(sigma_1 - sigma_1_nom,2) + pow(sigma_VBF - sigma_VBF_nom,2));
        double binContent = 0.;
        if(sigma_0_nom == 0){
          if(syst_ == "QCDscale_up") binContent = sigma_0_nom + delta_sigma_0;
          else if(syst_ == "QCDscale_down") binContent = sigma_0_nom - delta_sigma_0;
        }
        else if(sigma_0 > sigma_0_nom) binContent = sigma_0_nom + delta_sigma_0;
        else binContent = sigma_0_nom - delta_sigma_0;
        mon_.setBinContentAndError("mT_final_eq0jets", tagsR_[lepCat].substr(1), bin, binContent, 0., divideFinalHistoByBinWidth_);
      }
    }
  }

  mon_.WriteForSysts(syst_,keepAllControlPlots_);
  outFile->Close();

}


bool MainAnalysis::ProcessEvent() {
  if (not ApplyCommonFilters())
    return false;

  evt currentEvt;

  double weight = 1.;
  //get the MC event weight if exists
  if (isSim_) {
    weight *= (*genWeight_)() * intLumi_;

    //get the PU weights
    weight *= (*pileUpWeight_)();
  }

  mon_.fillHisto("eventflow","tot",0,weight);

  // Remove events with 0 vtx
  if(*numPVGood_ == 0 )
    return false;

  for (int i = 0; i < int(muonPt_.GetSize()); i++)
    mon_.fillHisto("pT_mu", "tot", muonPt_[i], weight);
  for (int i = 0; i < int(electronPt_.GetSize()); i++)
    mon_.fillHisto("pT_e", "tot", electronPt_[i], weight);
  mon_.fillHisto("nb_mu","tot",muonPt_.GetSize(),weight);
  mon_.fillHisto("nb_e","tot",electronPt_.GetSize(),weight);
  mon_.fillHisto("reco-vtx","tot",*numPVGood_,weight);


  //###############################################################
  //##################     OBJECT CORRECTIONS    ##################
  //###############################################################
   
  // MELA weight and kfactor
  weight *= melaWeight_();
  if(isSim_) weight *= (*kFactorCorrection_)();

  // electroweak corrections
  if(isSim_) weight *= (*ewCorrectionWeight_)();

  // Theory uncertainties
  double thUncWeight = 1.;

  if (syst_ != "")
    thUncWeight = utils::getTheoryUncertainties(*genWeight_, syst_);

  if(thUncWeight == 0) return false; // There are some rare cases where a weight is at 0, making an indeterminate form (0/0) in the code. I figured out it was an easy (albeit a bit coward) way to avoid it without changing all the code for an effect of less than 0.01%.
  weight *= thUncWeight;

  LOG_TRACE << "Weight after corrections: " << weight;

  //###############################################################
  //##################     OBJECT SELECTION      ##################
  //###############################################################

  auto const &tightElectrons = electronBuilder_.GetTight();
  auto const &looseElectrons = electronBuilder_.GetLoose();

  auto const &tightMuons = muonBuilder_.GetTight();
  auto const &looseMuons = muonBuilder_.GetLoose();

  auto const &photons = photonBuilder_.Get();
  auto const &jets = jetBuilder_.Get();

  //Discriminate ee and mumu
  bool isEE = (tightElectrons.size() >= 2 && !isPhotonDatadriven_); //2 good electrons
  bool isMuMu = (tightMuons.size() >= 2 && !isPhotonDatadriven_); //2 good muons
  bool isGamma = (photons.size() == 1 && isPhotonDatadriven_); //1 good photon

  mon_.fillHisto("nb_mu", "sel",
                std::min<int>(tightMuons.size(), 2), weight);
  mon_.fillHisto("nb_e", "sel",
                std::min<int>(tightElectrons.size(), 2), weight);
  mon_.fillHisto("nb_mu", "extra",
                looseMuons.size() - std::min<int>(tightMuons.size(), 2),
                weight);
  mon_.fillHisto("nb_e", "extra",
                looseElectrons.size() -
                  std::min<int>(tightElectrons.size(), 2),
                weight);

  //###############################################################
  //##################       ANALYSIS CUTS       ##################
  //###############################################################

  if(!isEE && !isMuMu && !isGamma) //not a good lepton pair or photon (if datadriven)
    return false;
  mon_.fillHisto("eventflow","tot",1,weight);

  if(isEE) currentEvt.s_lepCat = "_ee";
  else if(isMuMu) currentEvt.s_lepCat = "_mumu";

  //compute and apply the efficiency SFs
  if (isSim_) {
    if (not isPhotonDatadriven_) {  // Leptons
      weight *= leptonWeight_();
    }
    else {  // Photons
      weight *= photonWeight_();
    }
  }

  //trigger weights for photon data
  if(isPhotonDatadriven_){
    int triggerPrescale = photonPrescales_.GetPhotonPrescale(photons[0].p4.Pt());
    if (triggerPrescale == 0)  //trigger not found
      return false;

    weight *= triggerPrescale;
  }

  if(isPhotonDatadriven_){
    //Avoid double counting for NLO ZvvG:
    if (isMC_NLO_ZGTo2NuG_inclusive_ && photons[0].p4.Pt() >= 130)
      return false;
    if (isMC_NLO_ZGTo2NuG_Pt130_ && photons[0].p4.Pt() < 130)
      return false;

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

  currentEvt.s_jetCat = v_jetCat_[jetCat];

  auto const &ptMiss = ptMissBuilder_.Get();
  TLorentzVector const ptMissP4 = ptMissBuilder_.Get().p4;

  //Loop on lepton type. This is important also to apply Instr.MET if needed:
  double weightBeforeLoop = weight;
  TLorentzVector bosonBeforeLoop = boson;
  double photon_reweighting_tot = 1.;
  bool eventAccepted = false;

  for(unsigned int c = 0; c < tagsR_size_; c++){
    weight = weightBeforeLoop;
    boson = bosonBeforeLoop;
    photon_reweighting_tot = 1.;
    int lepCat = c;

    if(!isPhotonDatadriven_){
      if(tagsR_[c] == "_ee" && !isEE)
        continue;
      else if(tagsR_[c] == "_mumu" && !isMuMu)
        continue;
    }
    else{
      //Apply photon reweighting
      //1. #Vtx
      std::map<std::pair<double,double>, std::pair<double, double> >::iterator itlow;
      itlow = nVtxWeight_map_[tagsR_[c]].upper_bound(std::make_pair(*numPVGood_,boson.Pt())); //look at which bin in the map currentEvt.rho corresponds
      if (itlow == nVtxWeight_map_[tagsR_[c]].begin())
        throw HZZException{
          "You are trying to access your NVtx reweighting map outside of bin "
          "boundaries."
        };
      itlow--;
      weight *= itlow->second.first; //(itlow->second.first = reweighting value; itlow->second.second = reweighting error)
      //2. Pt
      std::map<double, std::pair<double, double> >::iterator itlow2;
      itlow2 = ptWeight_map_[tagsR_[c]+v_jetCat_[jetCat]].upper_bound(boson.Pt()); //look at which bin in the map currentEvt.pT corresponds
      if (itlow2 == ptWeight_map_[tagsR_[c] + v_jetCat_[jetCat]].begin())
        throw HZZException{
          "You are trying to access your Pt reweighting map outside of bin "
          "boundaries."
        };
      itlow2--;
      weight *= itlow2->second.first; 
      //3. Mass lineshape
      utils::giveMassToPhoton(boson, lineshapeMassWeight_map_[tagsR_[c]]);

      photon_reweighting_tot = 1.*itlow->second.first*itlow2->second.first;
      if(photon_reweighting_tot == 0)
        continue;
    }
    //Jet category

    //Warning, starting from here ALL plots have to have the currentEvt.s_lepCat in their name, otherwise the reweighting will go crazy
    currentEvt.Fill_evt(
      v_jetCat_[jetCat], tagsR_[c], boson, ptMissP4, jets, *run_,
      *numPVGood_, *rho_, ptMiss.significance, tightLeptons);

    mon_.fillHisto("jetCategory","tot"+currentEvt.s_lepCat,jetCat,weight);
    mon_.fillHisto("nJets","tot"+currentEvt.s_lepCat,currentEvt.nJets,weight);

    if(currentEvt.s_lepCat == "_ll") mon_.fillAnalysisHistos(currentEvt, "tot", weight);

    if(fabs(boson.M()-91.1876) > zMassWindow_) //We keep this cut for the photon... because we've already given a mass to the photon :)
      continue;
    if(currentEvt.s_lepCat == "_ll") mon_.fillHisto("eventflow","tot",2,weight);
    mon_.fillAnalysisHistos(currentEvt, "inZpeak", weight);

    if(boson.Pt() < minPtLL_)
      continue;
    if(currentEvt.s_lepCat == "_ll") mon_.fillHisto("eventflow","tot",3,weight);

    //Phi(Z,MET)
    if(currentEvt.deltaPhi_MET_Boson < minDphiLLPtMiss_)
      continue;
    if(currentEvt.s_lepCat == "_ll") mon_.fillHisto("eventflow","tot",4,weight);


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

      if (tabulatedRng_.Rndm(0) < std::pow(0.04, numExtraLeptons))
        // Estimate a 4% uncertainty
        passLeptonVeto = true;

      if (syst_ == "lepveto_down")
        weight *= -1.;
    }

    if (not passLeptonVeto)
      continue;

    // Reject event with any extra iso-track
    if (isotrkBuilder_.Get().size() > 0)
      return false;

    if (currentEvt.s_lepCat == "_ll")
      mon_.fillHisto("eventflow", "tot", 5, weight);

    // Compute the btagging efficiency
    if (isSim_)
      FillBTagEfficiency(jets, weight);

    // b veto
    bool passBTag = true;

    for (auto const &jet : jets)
      if (bTagger_(jet)) {
        passBTag = false;
        break;
      }

    if (not passBTag)
      continue;

    // Apply the btag weights
    if (isSim_) {
      double const w = bTagWeight_();
      weight *= w;

      if (currentEvt.s_lepCat == "_ll")
        mon_.fillProfile("BTagWeightvsMT", "TEST", currentEvt.MT, w, weight);
    }

    if(currentEvt.s_lepCat == "_ll") mon_.fillHisto("eventflow","tot",6,weight);

    // Phi(jet,MET)
    bool passDeltaPhiJetMET = true;

    for (auto const &jet : jets)
      if (std::abs(utils::deltaPhi(jet.p4, ptMissP4)) < minDphiJetsPtMiss_) {
        passDeltaPhiJetMET = false;
        break;
      }

    if (not passDeltaPhiJetMET)
      continue;

    if (DPhiPtMiss({&jetBuilder_, &muonBuilder_, &electronBuilder_})
        < minDphiLeptonsJetsPtMiss_)
      continue;

    if(currentEvt.s_lepCat == "_ll") mon_.fillHisto("eventflow","tot",7,weight);

    mon_.fillInstrMETControlRegionHisto(currentEvt, "InstrMET_reweighting", weight);
    mon_.fillAnalysisHistos(currentEvt, "beforeMETcut", weight);
    mon_.fillHisto("jetCategory","beforeMETcut"+currentEvt.s_lepCat,jetCat,weight);

    //MET>80
    if(ptMissP4.Pt()<80)
      continue;
    if(currentEvt.s_lepCat == "_ll") mon_.fillHisto("eventflow","tot",8,weight);

    mon_.fillAnalysisHistos(currentEvt, "baselineSelection", weight, divideFinalHistoByBinWidth_);

    //MET>125
    if(ptMissP4.Pt()<125)
      continue;
    if(currentEvt.s_lepCat == "_ll") mon_.fillHisto("eventflow","tot",9,weight);

    eventAccepted = true;

    //###############################################################
    //##################     END OF SELECTION      ##################
    //###############################################################
    mon_.fillHisto("jetCategory","finalSelection"+currentEvt.s_lepCat,jetCat,weight, divideFinalHistoByBinWidth_);
    mon_.fillAnalysisHistos(currentEvt, "finalSelection", weight, divideFinalHistoByBinWidth_);
    //THE RESULT:
    mon_.fillHisto("mT_final"+currentEvt.s_jetCat, tagsR_[c].substr(1), currentEvt.MT, weight, divideFinalHistoByBinWidth_);

    //Uncertainties
    if(syst_ == "QCDscale_up" || syst_ == "QCDscale_down"){
      mon_.fillHisto("mT_finalBinning0j"+currentEvt.s_jetCat, tagsR_[c].substr(1), currentEvt.MT, weight, divideFinalHistoByBinWidth_);
      mon_.fillHisto("mT_finalBinning0j"+currentEvt.s_jetCat, tagsR_[c].substr(1)+"_nominal", currentEvt.MT, weight/thUncWeight, divideFinalHistoByBinWidth_);
    }

    if ((syst_ == "pdf_up" || syst_ == "pdf_down")
        && currentEvt.s_lepCat != "_ll") {
      for (int i = 0 ; i < 100; i++)
        pdfReplicas_.at(jetCat).at(lepCat).at(i)->Fill(
          currentEvt.MT, weight * genWeight_->RelWeightPdf(i));
    }


    //Prepare the correct computation of the stat uncertainty for the mT plots with Instr.MET.
    if(isPhotonDatadriven_){
      int mT = std::min(h_mT_size_[jetCat], h_mT_[jetCat]->FindBin(currentEvt.MT));
      int Vtx = std::min(h_Vtx_->GetNbinsX(), h_Vtx_->FindBin(currentEvt.nVtx));
      int pT_threshold = std::min(h_pT_thresholds_->GetNbinsX() , h_pT_thresholds_->FindBin(currentEvt.pT_Boson));
      int pT = std::min(h_pT_->GetNbinsX() , h_pT_->FindBin(currentEvt.pT_Boson));
      (*mT_InstrMET_map_)[lepCat][jetCat][mT][Vtx][pT_threshold][pT].first += 1.*weight/photon_reweighting_tot; //Fill with the weight before photon reweighting
      (*mT_InstrMET_map_)[lepCat][jetCat][mT][Vtx][pT_threshold][pT].second += 1.*weight*weight/(photon_reweighting_tot*photon_reweighting_tot); //Errors due to weights without photon reweighting
    }
  }

  return eventAccepted;
}


void MainAnalysis::FillBTagEfficiency(
    std::vector<Jet> const &jets, double weight) const {

  for (auto const &jet : jets) {
    if (not bTagger_.IsTaggable(jet))
      continue;

    std::string tag;

    if (std::abs(jet.hadronFlavour) == 5)
      tag = "bjet";
    else if (std::abs(jet.hadronFlavour) == 4)
      tag = "cjet";
    else
      tag = "udsgjet";

    auto fillHistogram = [&](std::string const &label) {
      mon_.fillHisto("btagEff", label, jet.p4.Pt(), jet.p4.Eta(), weight);
    };

    fillHistogram("den_" + tag);

    if (bTagger_(jet))
      fillHistogram("num_" + tag);
  }
}


void MainAnalysis::InitializeHistograms(Options const &options)
{
  mon_.declareHistos();

  //Definition of the final histos (and in particular of the mT binning
  h_mT_[eq0jets] = (TH1*) mon_.getHisto("mT_final_eq0jets", tagsR_[ee].substr(1), divideFinalHistoByBinWidth_); h_mT_size_[eq0jets] = h_mT_[eq0jets]->GetNbinsX();
  h_mT_[geq1jets] = (TH1*) mon_.getHisto("mT_final_geq1jets", tagsR_[ee].substr(1), divideFinalHistoByBinWidth_); h_mT_size_[geq1jets] = h_mT_[geq1jets]->GetNbinsX();
  h_mT_[vbf] = (TH1*) mon_.getHisto("mT_final_vbf", tagsR_[ee].substr(1), divideFinalHistoByBinWidth_); h_mT_size_[vbf] = h_mT_[vbf]->GetNbinsX();
  mon_.getHisto("mT_final_eq0jets", tagsR_[mumu].substr(1), divideFinalHistoByBinWidth_); //The .substr(1) removes the annoying _ in the tagsR_ definition.
  mon_.getHisto("mT_final_geq1jets", tagsR_[mumu].substr(1), divideFinalHistoByBinWidth_);
  mon_.getHisto("mT_final_vbf", tagsR_[mumu].substr(1), divideFinalHistoByBinWidth_);
  int h_mT_maxSize = std::max({h_mT_size_[eq0jets], h_mT_size_[geq1jets], h_mT_size_[vbf]}); 

  // Histograms for PDF replicas
  if(syst_=="pdf_up" || syst_=="pdf_down"){
    for(unsigned int i = 0 ; i < v_jetCat_.size() ; i++){
      std::vector<std::vector<TH1*>> currentVectorOfVectors;
      for(unsigned int j = 0 ; j < tagsR_.size() -1 ; j++){
        std::vector<TH1*> currentVectorOfReplicas;
        for(unsigned int k = 0 ; k < 100 ; k++){
          currentVectorOfReplicas.push_back((TH1*) mon_.getHisto("mT_final"+v_jetCat_.at(i), tagsR_.at(j).substr(1), divideFinalHistoByBinWidth_)->Clone((TString) "mT"+v_jetCat_.at(i)+tagsR_.at(j)+"_"+std::to_string(k)));
        }
        currentVectorOfVectors.push_back(currentVectorOfReplicas);
      }
      pdfReplicas_.push_back(currentVectorOfVectors);
    }
  }

  // ***--- Instr. MET building ---***
  //Compute once weights for Instr. MET reweighting if needed
  applyNvtxWeights_ = Options::NodeAs<bool>(
    options.GetConfig(), {"photon_reweighting", "apply_nvtx_reweighting"});
  applyPtWeights_ = Options::NodeAs<bool>(
    options.GetConfig(), {"photon_reweighting", "apply_pt_reweighting"});
  applyMassLineshape_ = Options::NodeAs<bool>(
    options.GetConfig(), {"photon_reweighting", "apply_mass_lineshape"});
  
  h_Vtx_ = (TH1*) mon_.getHisto("reco-vtx", "toGetBins", false);
  int const h_Vtx_size = h_Vtx_->GetNbinsX();
  h_pT_thresholds_ = (TH1*) mon_.getHisto("pT_Boson_thresholds", "toGetBins", false);
  int const h_pT_thresholds_size = h_pT_thresholds_->GetNbinsX();
  h_pT_ = (TH1*) mon_.getHisto("pT_Boson", "toGetBins", false);
  int const h_pT_size = h_pT_->GetNbinsX();
  mT_InstrMET_map_.reset(new decltype(mT_InstrMET_map_)::element_type (lepCat_size, std::vector<std::vector<std::vector<std::vector<std::vector<std::pair<double, double> > > > > >(jetCat_size, std::vector<std::vector<std::vector<std::vector<std::pair<double, double> > > > >(h_mT_maxSize+1, std::vector<std::vector<std::vector<std::pair<double, double> > > >(h_Vtx_size+1, std::vector< std::vector<std::pair<double, double> > >(h_pT_thresholds_size+1, std::vector<std::pair<double, double> >(h_pT_size+1)))))));
  photon_reweighting_.reset(new decltype(photon_reweighting_)::element_type (
      lepCat_size, std::vector<std::vector<std::vector<std::vector<std::pair<double, double> > > > >(jetCat_size, std::vector<std::vector<std::vector<std::pair<double, double> > > >(h_Vtx_size+1, std::vector<std::vector<std::pair<double, double> > >(h_pT_thresholds_size+1, std::vector<std::pair<double, double> >(h_pT_size+1))))));

  if (isPhotonDatadriven_
      && (!applyNvtxWeights_ || !applyPtWeights_ || !applyMassLineshape_))
    throw HZZException{
      "You tried to run datadriven method without having weights for "
      "Instr.MET. This is bad :-) Please compute weights first!"
    };
  if(isPhotonDatadriven_){
    utils::loadInstrMETWeights(applyNvtxWeights_, applyPtWeights_, applyMassLineshape_, nVtxWeight_map_, ptWeight_map_, lineshapeMassWeight_map_, v_jetCat_, options);
    for (int lepCat = 0; lepCat < int(tagsR_.size()); lepCat++) {
      for (int jetCat = 0; jetCat < int(v_jetCat_.size()); jetCat++) {
        for (int Vtx = 1; Vtx <= h_Vtx_size; Vtx++){
          for (int pT_threshold = 1; pT_threshold <= h_pT_thresholds_size; pT_threshold++) {
            //1. #Vtx
            std::map<std::pair<double,double>, std::pair<double, double> >::iterator Vtx_low;
            Vtx_low = nVtxWeight_map_[tagsR_[lepCat]].upper_bound(std::make_pair(h_Vtx_->GetBinLowEdge(Vtx),h_pT_thresholds_->GetBinCenter(pT_threshold))); //look at which bin in the map this nVtx corresponds
            if(Vtx_low == nVtxWeight_map_[tagsR_[lepCat]].begin()) continue;
            Vtx_low--;
            for (int pT = 1; pT <= h_pT_size; pT++) {
              //2. Pt
              std::map<double, std::pair<double, double> >::iterator pT_low;
              pT_low = ptWeight_map_[tagsR_[lepCat]+v_jetCat_[jetCat]].upper_bound(h_pT_->GetBinCenter(pT)); //look at which bin in the map this pT corresponds
              if(pT_low == ptWeight_map_[tagsR_[lepCat]+v_jetCat_[jetCat]].begin()) continue;
              pT_low--;
              (*photon_reweighting_)[lepCat][jetCat][Vtx][pT_threshold][pT].first = Vtx_low->second.first * pT_low->second.first;
              (*photon_reweighting_)[lepCat][jetCat][Vtx][pT_threshold][pT].second = sqrt(Vtx_low->second.second*Vtx_low->second.second*pT_low->second.first*pT_low->second.first + pT_low->second.second*pT_low->second.second*Vtx_low->second.first*Vtx_low->second.first);
            }
          }
        }
      }
    }
  }
}
