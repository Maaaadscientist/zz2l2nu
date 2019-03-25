#include <ObjectSelection.h>

#include <algorithm>
#include <limits>

#include <Utils.h>

namespace objectSelection
{

  bool selectJets(std::vector<TLorentzVectorWithIndex> & selJets, std::vector<TLorentzVectorWithIndex> & selCentralJets, std::vector<double> & btags, TTreeReaderArray<float> const &JetAk04Pt, TTreeReaderArray<float> const &JetAk04Eta, TTreeReaderArray<float> const &JetAk04Phi, TTreeReaderArray<float> const &JetAk04E, TTreeReaderArray<float> const &JetAk04Id, TTreeReaderArray<float> const &JetAk04NeutralEmFrac, TTreeReaderArray<float> const &JetAk04NeutralHadAndHfFrac, TTreeReaderArray<float> const &JetAk04NeutMult, TTreeReaderArray<float> const &JetAk04BDiscCisvV2, std::vector<Muon> const &selMuons, std::vector<Electron> const &selElectrons, std::vector<Photon> const &selPhotons)
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

      double minDRgj = std::numeric_limits<double>::infinity();

      for (auto const &ph : selPhotons)
        minDRgj = std::min(minDRgj, utils::deltaR(currentJet, ph.p4));

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
