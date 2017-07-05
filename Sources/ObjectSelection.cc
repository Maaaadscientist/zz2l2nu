#include "../Includes/ObjectSelection.h"

namespace objectSelection
{

  bool selectElectrons(std::vector<TLorentzVector> & selElectrons, std::vector<TLorentzVector> & extraElectrons, std::vector<float> *ElPt, std::vector<float> *ElEta, std::vector<float> *ElPhi, std::vector<float> *ElE, std::vector<unsigned int> *ElId, std::vector<float> *ElEtaSc, std::vector<float> *ElPfIsoRho)
  {
    for(int i = 0 ; i<ElPt->size() ; i++){
      bool passEta = false, passIso = false, passId = false, passPt = false, passLoosePt = false, passLooseId = false;
      TLorentzVector currentLepton; currentLepton.SetPtEtaPhiE(ElPt->at(i),ElEta->at(i),ElPhi->at(i),ElE->at(i));
       //Id //Very temporary!!! Used without much cross-checking.
      passId = ElId->at(i) & (1<<17);
      passLooseId = ElId->at(i) & (1<<16);
      int eta = fabs(ElEtaSc->at(i));//I took the supercluster eta since it's really the geometry which is taken here.
      passEta = (eta<=2.5 && (eta>=1.5660 || eta<=1.4442));
      //Iso //We use ElPfIsoRho for now, we'll see after if it's mandatory to refine it. Iso is applied only for the "tight" selection, not for the extra lepton veto.
      if(eta>=1.5660 && ElPfIsoRho->at(i)<0.0646) passIso = true;
      if(eta<=1.4442 && ElPfIsoRho->at(i)<0.0354) passIso = true; //Numbers are taken from llvv_fwk and have not been checked.
      passPt = (currentLepton.Pt() >=25);
      passLoosePt = (currentLepton.Pt() >=10);
      if(passEta && passLooseId && passLoosePt && selElectrons.size()==2) extraElectrons.push_back(currentLepton); //No iso criteria for extra leptons.
      if(passEta && passIso && passId && passPt && selElectrons.size()<2) selElectrons.push_back(currentLepton);
    }
    return true;
  }

  bool selectMuons(std::vector<TLorentzVector> & selMuons, std::vector<TLorentzVector> & extraMuons, std::vector<float> *MuPt, std::vector<float> *MuEta, std::vector<float> *MuPhi, std::vector<float> *MuE, std::vector<unsigned int> *MuId, std::vector<unsigned int> *MuIdTight, std::vector<float> *MuPfIso)
  {
    for(int i = 0 ; i<MuPt->size() ; i++){
      bool passEta = false, passIso = false, passId = false, passPt = false, passLoosePt = false, passLooseId = false, passSoftId = false, passSoftPt = false;
      TLorentzVector currentLepton; currentLepton.SetPtEtaPhiE(MuPt->at(i),MuEta->at(i),MuPhi->at(i),MuE->at(i));
      //Id //Very temporary!!! Used without much cross-checking.
      passId = MuIdTight->at(i) & (1<<0); //Look at the first vertex, hence the bit 0.
      passLooseId = MuId->at(i) & (1<<0);
      passSoftId = false; //Don't know yet how to implement this thing. FIXME
      int eta = fabs(MuEta->at(i));
      passEta = (eta<=2.4);
      //Iso //We use MuPfIso for now, we'll see after if it's mandatory to refine it. Iso is applied only for the "tight" selection, not for the extra lepton veto.
      passIso = (MuPfIso->at(i)<0.15); //Numbers are taken from llvv_fwk and have not been checked.
      passPt = (currentLepton.Pt() >=25);
      passLoosePt = (currentLepton.Pt() >=10);
      passSoftPt = (currentLepton.Pt() >=3);
      if(passEta && passLooseId && passLoosePt && selMuons.size()==2) extraMuons.push_back(currentLepton); //No iso criteria for extra leptons.
      if(passEta && !(passLooseId && passLoosePt) && passSoftId && passSoftPt && selMuons.size()==2) extraMuons.push_back(currentLepton); //Soft leptons. Need a particular cut?
      if(passEta && passIso && passId && passPt && selMuons.size()<2) selMuons.push_back(currentLepton);
    }
    return true;
  }

  bool selectPhotons(std::vector<TLorentzVector> & selPhotons, std::vector<float> *PhotPt, std::vector<float> *PhotEta, std::vector<float> *PhotPhi, std::vector<unsigned int> *PhotId, std::vector<float> *PhotScEta, std::vector<TLorentzVector> & selMuons, std::vector<TLorentzVector> & selElectrons)
  {
    for(int i = 0 ; i<PhotPt->size() ; i++){
      bool passId = false, passPt = false, passEta = false, passLeptonCleaning = false;
      TLorentzVector currentPhoton; currentPhoton.SetPtEtaPhiE(PhotPt->at(i),PhotEta->at(i),PhotPhi->at(i),utils::getPhotonEnergy(PhotPt->at(i),PhotEta->at(i))); //photon energy is completely given by Pt and Eta.
      passId = PhotId->at(i) & (1<<2); //tight, according to llvv_fwk the code. FIXME: check that it's not better to redefine everything ourselves.
      passPt = (currentPhoton.Pt() >= 55);
      passEta = (fabs(PhotScEta->at(i))<=1.4442);
      double minDRlg(9999.); for(int ilep=0; ilep<selMuons.size(); ilep++) minDRlg = TMath::Min( minDRlg, utils::deltaR(currentPhoton,selMuons[ilep]) );
      for(int ilep=0; ilep<selElectrons.size(); ilep++) minDRlg = TMath::Min( minDRlg, utils::deltaR(currentPhoton,selElectrons[ilep]) );
      passLeptonCleaning = (minDRlg>=0.1); //according to the llvv_fwk code.
      if(passId && passPt && passEta && passLeptonCleaning) selPhotons.push_back(currentPhoton);
    }
    return true;
  }

  bool selectJets(std::vector<TLorentzVector> & tagJets, std::vector<TLorentzVector> & selJets, std::vector<double> & btags, std::vector<float> *JetAk04Pt, std::vector<float> *JetAk04Eta, std::vector<float> *JetAk04Phi, std::vector<float> *JetAk04E, std::vector<float> *JetAk04Id, std::vector<float> *JetAk04NeutralEmFrac, std::vector<float> *JetAk04NeutralHadAndHfFrac, std::vector<float> *JetAk04BDiscCisvV2, const std::vector<TLorentzVector> & selMuons, const std::vector<TLorentzVector> & selElectrons, const std::vector<TLorentzVector> & selPhotons)
  {
    for(int i =0 ; i<JetAk04Pt->size() ; i++){
      bool passPt = false, passSelPt = false, passEta = false, passTightEta = false, passId = false, passLeptonCleaning = false, passPhotonCleaning = false;
      TLorentzVector currentJet; currentJet.SetPtEtaPhiE(JetAk04Pt->at(i),JetAk04Eta->at(i),JetAk04Phi->at(i),JetAk04E->at(i));
      passPt = (currentJet.Pt() >=15);
      passSelPt = (currentJet.Pt() >=30);
      double eta = fabs(currentJet.Eta());
      passEta = (eta <= 4.7);
      passTightEta = (eta <= 2.5);
      if(eta<2.7){
        passId = (JetAk04Id->at(i) >= 1); //This case is simple, it corresponds exactly to what we apply for now
        }
      float nef = JetAk04NeutralEmFrac->at(i);
      float nhf = JetAk04NeutralHadAndHfFrac->at(i);
      float nnp = 20; //FIXME this variable could not be found in the tree, it has obviously to be fixed if we are to reproduce what stands in the AN.
      if(eta<3.0 && eta >=2.7){
        //passId = (nef > 0.01 && nhf < 0.98 && nnp > 2);
        passId = (JetAk04Id->at(i) >= 2); //Simpler criterium, but not equivalent to what is mentionned in the AN
        }
      if(eta>=3.0){
        //passId = (nef<0.90 && nnp > 10);
        passId = (JetAk04Id->at(i) >= 3); //Simpler criterium, but not equivalent to what is mentionned in the AN
        }
      double minDRlj(9999.); for(int ilep=0; ilep<selMuons.size(); ilep++) minDRlj = TMath::Min( minDRlj, utils::deltaR(currentJet,selMuons[ilep]) );
      for(int ilep=0; ilep<selElectrons.size(); ilep++) minDRlj = TMath::Min( minDRlj, utils::deltaR(currentJet,selElectrons[ilep]) );
      passLeptonCleaning = (minDRlj>=0.4);
      double minDRgj(9999.); for(int ilep=0; ilep<selPhotons.size(); ilep++) minDRgj = TMath::Min( minDRgj, utils::deltaR(currentJet,selPhotons[ilep]) );
      passPhotonCleaning = (minDRgj>=0.4);
      if(passPt && passEta && passId && passLeptonCleaning && passPhotonCleaning) tagJets.push_back(currentJet);
      if(passSelPt && passEta && passId && passLeptonCleaning && passPhotonCleaning) selJets.push_back(currentJet);
      if(passSelPt && passTightEta && passId && passLeptonCleaning && passPhotonCleaning) btags.push_back(JetAk04BDiscCisvV2->at(i));
      }
    return true;
  }

}
