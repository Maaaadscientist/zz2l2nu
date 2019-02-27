#ifndef objectselection_h
#define objectselection_h

#include <string>
#include <vector>

#include <TTreeReaderArray.h>

#include <TLorentzVectorWithIndex.h>
#include <RoccoR.h>
#include <Utils.h>


namespace objectSelection
{

  bool selectElectrons(std::vector<TLorentzVectorWithIndex> & selElectrons, std::vector<TLorentzVectorWithIndex> & extraElectrons, TTreeReaderArray<float> const &ElPt, TTreeReaderArray<float> const &ElEta, TTreeReaderArray<float> const &ElPhi, TTreeReaderArray<float> const &ElE, TTreeReaderArray<unsigned> const &ElId, TTreeReaderArray<float> const &ElEtaSc);

  bool selectMuons(std::vector<TLorentzVectorWithIndex> & selMuons, std::vector<TLorentzVectorWithIndex> & extraMuons, std::vector<float> const &MuPt, TTreeReaderArray<float> const &MuEta, TTreeReaderArray<float> const &MuPhi, TTreeReaderArray<float> const &MuE, TTreeReaderArray<unsigned> const &MuId, TTreeReaderArray<unsigned> const &MuIdTight, TTreeReaderArray<unsigned> const &MuIdSoft, TTreeReaderArray<float> const &MuPfIso);

  bool selectPhotons(std::vector<TLorentzVectorWithIndex> & selPhotons, TTreeReaderArray<float> const &PhotPt, TTreeReaderArray<float> const &PhotEta, TTreeReaderArray<float> const &PhotPhi, TTreeReaderArray<unsigned> const &PhotId, TTreeReaderArray<float> const &PhotScEta, TTreeReaderValue<std::vector<bool>> &PhotHasPixelSeed, TTreeReaderArray<float> const &PhotSigmaIetaIeta, TTreeReaderArray<float> const &PhotSigmaIphiIphi, std::vector<TLorentzVectorWithIndex> & selMuons, std::vector<TLorentzVectorWithIndex> & selElectrons);

  bool selectPartiallyPhotons(std::vector<TLorentzVectorWithIndex> & selPhotons, TTreeReaderArray<float> const &PhotPt, TTreeReaderArray<float> const &PhotEta, TTreeReaderArray<float> const &PhotPhi, TTreeReaderArray<float> const &PhotScEta, TTreeReaderValue<std::vector<bool>> &PhotHasPixelSeed, std::vector<TLorentzVectorWithIndex> & selMuons, std::vector<TLorentzVectorWithIndex> & selElectrons, TTreeReaderArray<float> const &PhotHoE, TTreeReaderArray<float> const &PhotSigmaIetaIeta, TTreeReaderArray<float> const &PhotPfIsoChHad, TTreeReaderArray<float> const &PhotPfIsoNeutralHad, TTreeReaderArray<float> const &PhotPfIsoPhot, Float_t EvtFastJetRho, TString selectionLevel);

  bool selectJets(std::vector<TLorentzVectorWithIndex> & selJets, std::vector<TLorentzVectorWithIndex> & selCentralJets, std::vector<double> & btags, TTreeReaderArray<float> const &JetAk04Pt, TTreeReaderArray<float> const &JetAk04Eta, TTreeReaderArray<float> const &JetAk04Phi, TTreeReaderArray<float> const &JetAk04E, TTreeReaderArray<float> const &JetAk04Id, TTreeReaderArray<float> const &JetAk04NeutralEmFrac, TTreeReaderArray<float> const &JetAk04NeutralHadAndHfFrac, TTreeReaderArray<float> const &JetAk04NeutMult, TTreeReaderArray<float> const &JetAk04BDiscCisvV2, const std::vector<TLorentzVectorWithIndex> & selMuons, const std::vector<TLorentzVectorWithIndex> & selElectrons, const std::vector<TLorentzVectorWithIndex> & selPhotons);

  bool cleanPathologicEventsInPhotons(TString datasetName, float EvtRunNum, float EvtLumiNum, float EvtNum);
}

#endif
