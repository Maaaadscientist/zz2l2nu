#ifndef objectselection_h
#define objectselection_h

#include <iostream>
#include <string>
#include <TMath.h>
#include <vector>
#include "TLorentzVector.h"
#include "Utils.h"

namespace objectSelection
{

  bool selectElectrons(std::vector<TLorentzVector> & selElectrons, std::vector<TLorentzVector> & extraElectrons, std::vector<float> *ElPt, std::vector<float> *ElEta, std::vector<float> *ElPhi, std::vector<float> *ElE, std::vector<unsigned int> *ElId, std::vector<float> *ElEtaSc, std::vector<float> *ElPfIsoRho);

  bool selectMuons(std::vector<TLorentzVector> & selMuons, std::vector<TLorentzVector> & extraMuons, std::vector<float> *MuPt, std::vector<float> *MuEta, std::vector<float> *MuPhi, std::vector<float> *MuE, std::vector<unsigned int> *MuId, std::vector<unsigned int> *MuIdTight, std::vector<float> *MuPfIso);

  bool selectPhotons(std::vector<TLorentzVector> & selPhotons, std::vector<float> *PhotPt, std::vector<float> *PhotEta, std::vector<float> *PhotPhi, std::vector<unsigned int> *PhotId, std::vector<float> *PhotScEta, std::vector<TLorentzVector> & selMuons, std::vector<TLorentzVector> & selElectrons);

  bool selectJets(std::vector<TLorentzVector> & tagJets, std::vector<TLorentzVector> & selJets, std::vector<double> & btags, std::vector<float> *JetAk04Pt, std::vector<float> *JetAk04Eta, std::vector<float> *JetAk04Phi, std::vector<float> *JetAk04E, std::vector<float> *JetAk04Id, std::vector<float> *JetAk04NeutralEmFrac, std::vector<float> *JetAk04NeutralHadAndHfFrac, std::vector<float> *JetAk04BDiscCisvV2, const std::vector<TLorentzVector> & selMuons, const std::vector<TLorentzVector> & selElectrons, const std::vector<TLorentzVector> & selPhotons);

}

#endif
