#include <ObjectSelection.h>

#include <algorithm>
#include <limits>

#include <Utils.h>

namespace objectSelection
{

  bool selectPhotons(std::vector<TLorentzVectorWithIndex> & selPhotons, TTreeReaderArray<float> const &PhotPt, TTreeReaderArray<float> const &PhotEta, TTreeReaderArray<float> const &PhotPhi, TTreeReaderArray<unsigned> const &PhotId, TTreeReaderArray<float> const &PhotScEta, TTreeReaderValue<std::vector<bool>> &PhotHasPixelSeed, TTreeReaderArray<float> const &PhotSigmaIetaIeta, TTreeReaderArray<float> const &PhotSigmaIphiIphi, std::vector<Muon> const &selMuons, std::vector<Electron> const &selElectrons)
  {
    for(unsigned int i = 0 ; i<PhotPt.GetSize() ; i++){
      bool passId = false, passPt = false, passEta = false, passLeptonCleaning = false, passSpikes = false;
      TLorentzVectorWithIndex currentPhoton = TLorentzVectorWithIndex::PtEtaPhiEIndex(PhotPt[i],PhotEta[i],PhotPhi[i],utils::getPhotonEnergy(PhotPt[i],PhotEta[i]), i); //photon energy is completely given by Pt and Eta.
      passId = PhotId[i] & (1<<2); //tight, according to llvv_fwk the code. FIXME: check that it's not better to redefine everything ourselves.
      passPt = (currentPhoton.Pt() >= 55);
      passEta = (fabs(PhotScEta[i])<=1.4442);
      passSpikes = PhotSigmaIetaIeta[i] > 0.001 && PhotSigmaIphiIphi[i] > 0.001; //added to the ID atfer PhotonCR study.

      double minDRlg = std::numeric_limits<double>::infinity();

      for (auto const &l : selMuons)
        minDRlg = std::min(minDRlg, utils::deltaR(currentPhoton, l.p4));

      for (auto const &l : selElectrons)
        minDRlg = std::min(minDRlg, utils::deltaR(currentPhoton, l.p4));

      passLeptonCleaning = (minDRlg>=0.1); //according to the llvv_fwk code.
      if(passId && passPt && passEta && passLeptonCleaning && !(*PhotHasPixelSeed)[i] && passSpikes) selPhotons.push_back(currentPhoton); //We ask for no pixel seed for the photons.
    }
    return true;
  }

  bool selectPartiallyPhotons(std::vector<TLorentzVectorWithIndex> & selPhotons, TTreeReaderArray<float> const &PhotPt, TTreeReaderArray<float> const &PhotEta, TTreeReaderArray<float> const &PhotPhi, TTreeReaderArray<float> const &PhotScEta, TTreeReaderValue<std::vector<bool>> &PhotHasPixelSeed, std::vector<Muon> const &selMuons, std::vector<Electron> const &selElectrons, TTreeReaderArray<float> const &PhotHoE, TTreeReaderArray<float> const &PhotSigmaIetaIeta, TTreeReaderArray<float> const &PhotPfIsoChHad, TTreeReaderArray<float> const &PhotPfIsoNeutralHad, TTreeReaderArray<float> const &PhotPfIsoPhot, Float_t EvtFastJetRho, TString selectionLevel)
  {
    for(unsigned int i = 0 ; i<PhotPt.GetSize() ; i++){
      bool passId = false, passPt = false, passEta = false, passLeptonCleaning = false;
      TLorentzVectorWithIndex currentPhoton = TLorentzVectorWithIndex::PtEtaPhiEIndex(PhotPt[i],PhotEta[i],PhotPhi[i],utils::getPhotonEnergy(PhotPt[i],PhotEta[i]), i); //photon energy is completely given by Pt and Eta.
      //Tight ID by hand with available variables for barrel photon
      double pt = currentPhoton.Pt();
      double sceta = PhotScEta[i];
      passId = true;
      if(selectionLevel.Contains("1") && PhotHoE[i] > 0.0269) passId = false; //HoE
      if(selectionLevel.Contains("2") && PhotSigmaIetaIeta[i] > 0.00994) passId = false;
      if(selectionLevel.Contains("3") && utils::photon_rhoCorrectedIso(PhotPfIsoChHad[i], EvtFastJetRho, sceta, "chIso") > 0.202) passId = false;
      if(selectionLevel.Contains("4") && utils::photon_rhoCorrectedIso(PhotPfIsoNeutralHad[i], EvtFastJetRho, sceta, "nhIso") > 0.264+0.0148*pt+0.000017*pt*pt) passId = false;
      if(selectionLevel.Contains("5") && utils::photon_rhoCorrectedIso(PhotPfIsoPhot[i], EvtFastJetRho, sceta, "gIso") > 2.362+0.0047*pt) passId = false;
      passPt = (pt >= 55);
      passEta = (fabs(sceta)<=1.4442);
      
      double minDRlg = std::numeric_limits<double>::infinity();

      for (auto const &l : selMuons)
        minDRlg = std::min(minDRlg, utils::deltaR(currentPhoton, l.p4));

      for (auto const &l : selElectrons)
        minDRlg = std::min(minDRlg, utils::deltaR(currentPhoton, l.p4));

      passLeptonCleaning = (minDRlg>=0.1); //according to the llvv_fwk code.
      if(passId && passPt && passEta && passLeptonCleaning && !(*PhotHasPixelSeed)[i]) selPhotons.push_back(currentPhoton); //We ask for no pixel seed for the photons.
    }
    return true;
  }

  bool selectJets(std::vector<TLorentzVectorWithIndex> & selJets, std::vector<TLorentzVectorWithIndex> & selCentralJets, std::vector<double> & btags, TTreeReaderArray<float> const &JetAk04Pt, TTreeReaderArray<float> const &JetAk04Eta, TTreeReaderArray<float> const &JetAk04Phi, TTreeReaderArray<float> const &JetAk04E, TTreeReaderArray<float> const &JetAk04Id, TTreeReaderArray<float> const &JetAk04NeutralEmFrac, TTreeReaderArray<float> const &JetAk04NeutralHadAndHfFrac, TTreeReaderArray<float> const &JetAk04NeutMult, TTreeReaderArray<float> const &JetAk04BDiscCisvV2, std::vector<Muon> const &selMuons, std::vector<Electron> const &selElectrons, const std::vector<TLorentzVectorWithIndex> & selPhotons)
  {
    for(unsigned int i =0 ; i<JetAk04Pt.GetSize() ; i++){
      bool passSelPt = false, passEta = false, passTightEta = false, passId = false, passLeptonCleaning = false, passPhotonCleaning = false;
      TLorentzVectorWithIndex currentJet = TLorentzVectorWithIndex::PtEtaPhiEIndex(JetAk04Pt[i],JetAk04Eta[i],JetAk04Phi[i],JetAk04E[i], i);
      passSelPt = (currentJet.Pt() >=30);
      double eta = fabs(currentJet.Eta());
      passEta = (eta <= 4.7);
      passTightEta = (eta <= 2.5);
      if(eta<2.7){
        passId = (JetAk04Id[i] >= 1); //This case is simple, it corresponds exactly to what we apply for now
        }
      float nef = JetAk04NeutralEmFrac[i];
      float nhf = JetAk04NeutralHadAndHfFrac[i];
      float nnp = JetAk04NeutMult[i];
      if(eta<3.0 && eta >=2.7){
        passId = (nef > 0.01 && nhf < 0.98 && nnp > 2);
        //passId = (JetAk04Id[i] >= 2); //Simpler criterium, but not equivalent to what is mentionned in the AN. Was applied before having access to nnp.
        }
      if(eta>=3.0){
        passId = (nef<0.90 && nnp > 10);
        //passId = (JetAk04Id[i] >= 3); //Simpler criterium, but not equivalent to what is mentionned in the AN. Was applied before having access to nnp.
        }

      double minDRlj = std::numeric_limits<double>::infinity();

      for (auto const &l : selMuons)
        minDRlj = std::min(minDRlj, utils::deltaR(currentJet, l.p4));

      for (auto const &l : selElectrons)
        minDRlj = std::min(minDRlj, utils::deltaR(currentJet, l.p4));

      passLeptonCleaning = (minDRlj>=0.4);
      double minDRgj(9999.); for(unsigned int ilep=0; ilep<selPhotons.size(); ilep++) minDRgj = TMath::Min( minDRgj, utils::deltaR(currentJet,selPhotons[ilep]) );
      passPhotonCleaning = (minDRgj>=0.4);
      if(passSelPt && passEta && passId && passLeptonCleaning && passPhotonCleaning) selJets.push_back(currentJet);
      if(passSelPt && passTightEta && passId && passLeptonCleaning && passPhotonCleaning){
        selCentralJets.push_back(currentJet);
        btags.push_back(JetAk04BDiscCisvV2[i]);
        }
      }
    return true;
  }

  bool cleanPathologicEventsInPhotons(TString datasetName, float EvtRunNum, float EvtLumiNum, float EvtNum){
    bool eventShouldBeCleaned = false;
    //We remove events in the spike for those samples 
    if(datasetName.Contains("QCD_HT100to200")){
      //Spike at 190~200GeV (cut if pt > 190GeV)
      if( EvtRunNum == 1 && EvtLumiNum == 21997  && EvtNum == 32986438 ) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 123682 && EvtNum == 185472705) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 133696 && EvtNum == 200489234) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 301998 && EvtNum == 452875030) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 237717 && EvtNum == 356480214) eventShouldBeCleaned = true;
      
      //Spike at ~330-340GeV (cut if pt > 330GeV)
      if( EvtRunNum == 1 && EvtLumiNum == 405615 && EvtNum == 608258936) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 238040 && EvtNum == 356963627) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 192575 && EvtNum == 288784917) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 405110 && EvtNum == 607502440) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 398170 && EvtNum == 597094584) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 242217 && EvtNum == 363227739) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 175934 && EvtNum == 263829468) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 336765 && EvtNum == 505011533) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 252239 && EvtNum == 378257731) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 239877 && EvtNum == 359718643) eventShouldBeCleaned = true;
    }
    if(datasetName.Contains("QCD_Pt-30to50_EMEnriched")){
      //Spike at ~330-340GeV (cut if pt > 330GeV)
      if( EvtRunNum == 1 && EvtLumiNum == 5960 && EvtNum == 16066353) eventShouldBeCleaned = true;
    }
    if(datasetName.Contains("QCD_Pt-20toInf_MuEnrichedPt15")){
      //Spike at ~330-340GeV (cut if pt > 330GeV)
      if( EvtRunNum == 1 && EvtLumiNum == 181694 && EvtNum == 2097414398) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 40564  && EvtNum == 3384544677) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 49551  && EvtNum == 1742764771) eventShouldBeCleaned = true;
      if( EvtRunNum == 1 && EvtLumiNum == 94563  && EvtNum == 2145642832) eventShouldBeCleaned = true;
    }
    return eventShouldBeCleaned;
  }

}
