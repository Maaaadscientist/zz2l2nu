#include <NrbAnalysis.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <TFile.h>

#include <Utils.h>


namespace po = boost::program_options;


NrbAnalysis::NrbAnalysis(Options const &options, Dataset &dataset)
    : AnalysisCommon{options, dataset},
      dataset_{dataset},
      outputFile_{options.GetAs<std::string>("output")},
      keepAllControlPlots_{options.Exists("all-control-plots")},
      syst_{options.GetAs<std::string>("syst")},
      divideFinalHistoByBinWidth_{false},  //For final plots, we don't divide by the bin width to ease computations of the yields by eye.
      v_jetCat_{"_eq0jets","_geq1jets","_vbf"},
      tagsR_{"_ee", "_mumu", "_ll"}, tagsR_size_{unsigned(tagsR_.size())},
      fileName_{dataset_.Info().Files().at(0)}
{
  if (isSim_) {
    genPartPdgId_.reset(new TTreeReaderArray<int>(
        dataset_.Reader(), "GenPart_pdgId"));
    genPartMotherIndex_.reset(new TTreeReaderArray<int>(
        dataset_.Reader(), "GenPart_genPartIdxMother"));
  }

  InitializeHistograms();

  if (syst_ == "")
    LOG_DEBUG << "Will not apply systematic variations.";
  else
    LOG_DEBUG << "Will apply systematic variation \"" << syst_ << "\".";
}


po::options_description NrbAnalysis::OptionsDescription() {
  po::options_description optionsDescription{"Analysis-specific options"};
  optionsDescription.add_options()
    ("syst", po::value<std::string>()->default_value(""),
     "Requested systematic variation")
    ("all-control-plots", "Keep all control plots")
    ("output,o", po::value<std::string>()->default_value("outputFile.root"),
     "Name for output file with histograms");
  return optionsDescription;
}


void NrbAnalysis::PostProcessing() {
  TFile *outFile = TFile::Open(outputFile_.c_str(), "recreate");
  mon_.WriteForSysts(syst_, keepAllControlPlots_);
  outFile->Close();
}


void NrbAnalysis::InitializeHistograms() {
  mon_.declareHistos_NRB();

  optim_Cuts1_met_.push_back(0); //add a bin in the shapes with a MET cut of 0
  for(double met=50;met<140;met+=5) {  optim_Cuts1_met_.push_back(met);  }
  TH2F *h_2D=(TH2F *) mon_.addHistogram( new TH2F ("mt_shapes_NRBctrl",";cut Index;Selection region;Events",optim_Cuts1_met_.size(),0,optim_Cuts1_met_.size(),6,0,6) );

  h_2D->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
  h_2D->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
  h_2D->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
  h_2D->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
  h_2D->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
  h_2D->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");
  
  //Definition of the final histos (and in particular of the mT binning
  std::vector<TH1*> h_mT(jetCat_size); std::vector<int> h_mT_size(jetCat_size);
  h_mT[eq0jets] = (TH1*) mon_.getHisto("mT_final_eq0jets", "ee", divideFinalHistoByBinWidth_); h_mT_size[eq0jets] = h_mT[eq0jets]->GetNbinsX();
  h_mT[geq1jets] = (TH1*) mon_.getHisto("mT_final_geq1jets", "ee", divideFinalHistoByBinWidth_); h_mT_size[geq1jets] = h_mT[geq1jets]->GetNbinsX();
  h_mT[vbf] = (TH1*) mon_.getHisto("mT_final_vbf", "ee", divideFinalHistoByBinWidth_); h_mT_size[vbf] = h_mT[vbf]->GetNbinsX();
  mon_.getHisto("mT_final_eq0jets", "mumu", divideFinalHistoByBinWidth_); //The .substr(1) removes the annoying _ in the tagsR_ definition.
  mon_.getHisto("mT_final_geq1jets", "mumu", divideFinalHistoByBinWidth_);
  mon_.getHisto("mT_final_vbf", "mumu", divideFinalHistoByBinWidth_);
}


bool NrbAnalysis::ProcessEvent() {
  if (not meKinFilter_())
    return false;

  if (not metFilters_())
    return false;

  evt currentEvt;
  
  double weight = 1.;
  //get the MC event weight if exists
  if (isSim_) {
    //get the MC event weight if exists
    weight *= (*genWeight_)() * intLumi_;

    //get the PU weights
    weight *= (*pileUpWeight_)();
  }

  // Remove events with 0 vtx
  if (*numPVGood_ == 0 )
    return false;

  mon_.fillHisto("eventflow","tot",0,weight);


  //###############################################################
  //##################     OBJECT CORRECTIONS    ##################
  //###############################################################
  // electroweak corrections
  if(isSim_) weight *= (*ewCorrectionWeight_)();


  //###############################################################
  //##################     OBJECT SELECTION      ##################
  //###############################################################

  auto const &tightElectrons = electronBuilder_.GetTight();
  auto const &looseElectrons = electronBuilder_.GetLoose();

  auto const &tightMuons = muonBuilder_.GetTight();
  auto const &looseMuons = muonBuilder_.GetLoose();

  auto const &jets = jetBuilder_.Get();

  //Discriminate ee and mumu
  bool isEE = (tightElectrons.size() >= 2); //2 good electrons
  bool isMuMu = (tightMuons.size() >= 2); //2 good muons
  bool isEMu = (tightMuons.size() == 1 and tightElectrons.size() == 1);


  //###############################################################
  //##################       ANALYSIS CUTS       ##################
  //###############################################################

  
  for (int i = 0; i < int(muonPt_.GetSize()); i++)
    mon_.fillHisto("pT_mu", "tot", muonPt_[i], weight);
  for (int i = 0; i < int(electronPt_.GetSize()); i++)
    mon_.fillHisto("pT_e", "tot", electronPt_[i], weight);
  mon_.fillHisto("nb_mu","tot",muonPt_.GetSize(),weight);
  mon_.fillHisto("nb_e","tot",electronPt_.GetSize(),weight);  

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

  if(!isEE && !isMuMu && !isEMu)
    return false;
  
  
  if(isEE) currentEvt.s_lepCat = "_ee";
  else if(isMuMu) currentEvt.s_lepCat = "_mumu";
  else if (isEMu) currentEvt.s_lepCat = "_emu";
  if(isSim_&& (fileName_.Contains("DY")||fileName_.Contains("ZZTo2L")||fileName_.Contains("ZZToTauTau"))){
    int GLepId = 1;
    for (int i = 0; i < int(genPartPdgId_->GetSize()); i++) {
      if(genPartMotherIndex_->At(i) == 23) GLepId *= fabs(genPartPdgId_->At(i));
    }
    if (fileName_.Contains("DYJetsToTauTau")  &&   GLepId %5 != 0  )
      return true ;
    if (fileName_.Contains("ZZToTauTau2Nu") && GLepId % 5 != 0 )
      return true;
    if (fileName_.Contains("ZZToTauTau2Q") &&  GLepId % 5 != 0 )
      return true;
    if (fileName_.Contains("DYJetsToLL")  &&       GLepId %5 == 0  )
      return true ;
    if (fileName_.Contains("ZZTo2L2Nu")&&      GLepId % 5 == 0 )
      return true;
    if (fileName_.Contains("ZZTo2L2Q")&&       GLepId % 5 == 0 )
      return true;    
    //cout << " genLep Id product: "<< GLepId << endl;
  }
  
  //compute and apply the efficiency SFs
  if (isSim_)
    weight *= leptonWeight_();

  //Definition of the relevant analysis variables
  std::vector<Lepton> tightLeptons;

  if (isMuMu or isEMu) {
    for (auto const &mu : tightMuons)
      tightLeptons.emplace_back(mu);
  } 
  if (isEE or isEMu) {
    for (auto const &e : tightElectrons)
      tightLeptons.emplace_back(e);
  }

  if (isEMu)
    std::sort(tightLeptons.begin(), tightLeptons.end(), PtOrdered);

  TLorentzVector boson = tightLeptons[0].p4 + tightLeptons[1].p4;

  auto const &ptMiss = ptMissBuilder_.Get();
  TLorentzVector const ptMissP4 = ptMissBuilder_.Get().p4;

  //Loop on lepton type
  double weightBeforeLoop = weight;
  TLorentzVector bosonBeforeLoop = boson;
  bool eventAccepted = false;

  for(unsigned int c = 0; c < tagsR_size_; c++){
    weight = weightBeforeLoop;
    boson = bosonBeforeLoop;

    if(tagsR_[c] == "_ee" && !isEE) continue;
    else if(tagsR_[c] == "_mumu" && !isMuMu) continue;
    else if(tagsR_[c] == "_emu" && !isEMu) continue;
    else if(tagsR_[c] == "_ll" && !(isMuMu || isEE)) continue;

    //Jet category
    int jetCat = geq1jets;

    if (jets.size() == 0)
      jetCat = eq0jets;
    else if (utils::PassVbfCuts(jets, boson))
      jetCat = vbf;

    currentEvt.s_jetCat = v_jetCat_[jetCat];

    //Warning, starting from here ALL plots have to have the currentEvt.s_lepCat in their name, otherwise the reweighting will go crazy
    currentEvt.Fill_evt(
      v_jetCat_[jetCat], tagsR_[c], boson, ptMissP4, jets, *run_,
      *numPVGood_, *rho_, ptMiss.significance, tightLeptons);

    // Apply the btag weights
    if (isSim_)
      weight *= bTagWeight_();
    
    mon_.fillAnalysisHistos(currentEvt, "tot", weight);

    // b veto
    bool passBTag = true;

    for (auto const &jet : jets)
      if (bTagger_(jet)) {
        passBTag = false;
        break;
      }

    // Phi(jet,MET)
    bool passDeltaPhiJetMET = true;

    for (auto const &jet : jets)
      if (std::abs(utils::deltaPhi(jet.p4, ptMissP4)) < 0.5) {
        passDeltaPhiJetMET = false;
        break;
      }

    //DPhi
    bool passDphi(currentEvt.deltaPhi_MET_Boson>0.5);
    //boson
    bool passMass(fabs(currentEvt.M_Boson-91)<15);
    bool isZ_SB ( (currentEvt.M_Boson>40  && currentEvt.M_Boson<70) || (currentEvt.M_Boson>110 && currentEvt.M_Boson<200) );
    bool isZ_upSB ( (currentEvt.M_Boson>110 && currentEvt.M_Boson<200) );
    bool passQt (currentEvt.pT_Boson > 55.);

    unsigned const numExtraLeptons =
      looseMuons.size() - std::min<unsigned>(tightMuons.size(), 2) +
      looseElectrons.size() - std::min<unsigned>(tightElectrons.size(), 2);
    bool passThirdLeptonveto = ((isEE and tightMuons.empty()) or
      (isMuMu and tightElectrons.empty()) or isEMu) and
      numExtraLeptons == 0;
    
    TString tags = "tot"+currentEvt.s_lepCat; 

    if(currentEvt.M_Boson>40 && currentEvt.M_Boson<200 && passQt && passThirdLeptonveto  && passDeltaPhiJetMET && passDphi){
      if(passBTag)
      {
         if(ptMissP4.Pt()>50 )mon_.fillHisto("zmass_bveto50" , tags,currentEvt.M_Boson,weight);
         if(ptMissP4.Pt()>80 )mon_.fillHisto("zmass_bveto80" , tags,currentEvt.M_Boson,weight);
         if(ptMissP4.Pt()>125)mon_.fillHisto("zmass_bveto125", tags,currentEvt.M_Boson,weight);
         if(passMass)
         {
            mon_.fillHisto( "met_Inbveto",tags,ptMissP4.Pt(),weight);
            if(ptMissP4.Pt()>50 )mon_.fillHisto("mt_Inbveto50" , tags,currentEvt.MT,weight);
            if(ptMissP4.Pt()>80 )mon_.fillHisto("mt_Inbveto80" , tags,currentEvt.MT,weight);
            if(ptMissP4.Pt()>125){
              mon_.fillHisto("mt_Inbveto125", tags,currentEvt.MT,weight);
              mon_.fillHisto("mT_final"+currentEvt.s_jetCat, tagsR_[c].substr(1), currentEvt.MT, weight, divideFinalHistoByBinWidth_);
            }
         }
         else if(isZ_SB)
         {
            mon_.fillHisto( "met_Outbveto",tags,ptMissP4.Pt(),weight);
            if(ptMissP4.Pt()>50 )mon_.fillHisto("mt_Outbveto50" , tags,currentEvt.MT,weight);
            if(ptMissP4.Pt()>80 )mon_.fillHisto("mt_Outbveto80" , tags,currentEvt.MT,weight);
            if(ptMissP4.Pt()>125)mon_.fillHisto("mt_Outbveto125", tags,currentEvt.MT,weight);
         }
      }
      else
      {
        if(ptMissP4.Pt()>50 )mon_.fillHisto("zmass_btag50" , tags,currentEvt.M_Boson,weight);
        if(ptMissP4.Pt()>80 )mon_.fillHisto("zmass_btag80" , tags,currentEvt.M_Boson,weight);
        if(ptMissP4.Pt()>125)mon_.fillHisto("zmass_btag125", tags,currentEvt.M_Boson,weight);
        if(passMass)
        {
          mon_.fillHisto( "met_Inbtag",tags,ptMissP4.Pt(),weight);
          if(ptMissP4.Pt()>50 )mon_.fillHisto("mt_Inbtag50" , tags,currentEvt.MT,weight);
          if(ptMissP4.Pt()>80 )mon_.fillHisto("mt_Inbtag80" , tags,currentEvt.MT,weight);
          if(ptMissP4.Pt()>125)mon_.fillHisto("mt_Inbtag125", tags,currentEvt.MT,weight);
        }
        else if(isZ_SB)
        {
          mon_.fillHisto( "met_Outbtag",tags,ptMissP4.Pt(),weight);
          if(ptMissP4.Pt()>50 )mon_.fillHisto("mt_Outbtag50" , tags,currentEvt.MT,weight);
          if(ptMissP4.Pt()>80 )mon_.fillHisto("mt_Outbtag80" , tags,currentEvt.MT,weight);
          if(ptMissP4.Pt()>125)mon_.fillHisto("mt_Outbtag125", tags,currentEvt.MT,weight);
        }
      }

    }
    
    if(currentEvt.M_Boson>40 && currentEvt.M_Boson<200 && passQt && passThirdLeptonveto  && passDeltaPhiJetMET && passDphi)
    {
      for(unsigned int Index=0;Index<optim_Cuts1_met_.size();Index++)
      {
        if(ptMissP4.Pt()>optim_Cuts1_met_[Index])
        {
          if(passBTag && passMass)mon_.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 0.5,weight);
          if(passBTag && isZ_SB)mon_.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 1.5,weight);
          if(passBTag && isZ_upSB)mon_.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 2.5,weight);
          if(!passBTag && passMass)mon_.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 3.5,weight);
          if(!passBTag && isZ_SB)mon_.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 4.5,weight);
          if(!passBTag && isZ_upSB)mon_.fillHisto(TString("mt_shapes_NRBctrl"),tags,Index, 5.5,weight);
        }
      }
    }
    mon_.fillHisto("eventflow","tot",1,weight);
    mon_.fillHisto("eventflow",tags,1,weight);
    if(!passMass) continue;
    mon_.fillHisto("eventflow","tot",2,weight);
    mon_.fillHisto("eventflow",tags,2,weight);
    if(!passQt) continue;
    mon_.fillHisto("eventflow","tot",3,weight);
    mon_.fillHisto("eventflow",tags,3,weight);
     //Phi(Z,MET)
    if(!passThirdLeptonveto) continue;
    mon_.fillHisto("eventflow","tot",4,weight);
    mon_.fillHisto("eventflow",tags,4,weight);
    if(!passBTag) continue;
    mon_.fillHisto("eventflow","tot",5,weight);
    mon_.fillHisto("eventflow",tags,5,weight);
    if(!passDeltaPhiJetMET) continue;
    mon_.fillHisto("eventflow","tot",6,weight);
    mon_.fillHisto("eventflow",tags,6,weight);
    if(!passDphi) continue;
    mon_.fillHisto("eventflow","tot",7,weight);
    mon_.fillHisto("eventflow",tags,7,weight);

    mon_.fillAnalysisHistos(currentEvt, "beforeMETcut", weight);

    //MET>80
    if(ptMissP4.Pt()<80) continue;
    mon_.fillHisto("eventflow","tot",8,weight);
    mon_.fillHisto("eventflow",tags,8,weight);
    //MET>125
    if(ptMissP4.Pt()<125) continue;
    mon_.fillHisto("eventflow","tot",9,weight);
    mon_.fillHisto("eventflow",tags,9,weight);
    eventAccepted = true;
  }

  return eventAccepted;
}

