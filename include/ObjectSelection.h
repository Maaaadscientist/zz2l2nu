#ifndef objectselection_h
#define objectselection_h

#include <string>
#include <vector>

#include <TTreeReaderArray.h>

#include <PhysicsObjects.h>
#include <TLorentzVectorWithIndex.h>


namespace objectSelection
{

  bool selectPhotons(std::vector<TLorentzVectorWithIndex> & selPhotons, TTreeReaderArray<float> const &PhotPt, TTreeReaderArray<float> const &PhotEta, TTreeReaderArray<float> const &PhotPhi, TTreeReaderArray<unsigned> const &PhotId, TTreeReaderArray<float> const &PhotScEta, TTreeReaderValue<std::vector<bool>> &PhotHasPixelSeed, TTreeReaderArray<float> const &PhotSigmaIetaIeta, TTreeReaderArray<float> const &PhotSigmaIphiIphi, std::vector<Muon> const &selMuons, std::vector<Electron> const &selElectrons);

  bool selectJets(std::vector<TLorentzVectorWithIndex> & selJets, std::vector<TLorentzVectorWithIndex> & selCentralJets, std::vector<double> & btags, TTreeReaderArray<float> const &JetAk04Pt, TTreeReaderArray<float> const &JetAk04Eta, TTreeReaderArray<float> const &JetAk04Phi, TTreeReaderArray<float> const &JetAk04E, TTreeReaderArray<float> const &JetAk04Id, TTreeReaderArray<float> const &JetAk04NeutralEmFrac, TTreeReaderArray<float> const &JetAk04NeutralHadAndHfFrac, TTreeReaderArray<float> const &JetAk04NeutMult, TTreeReaderArray<float> const &JetAk04BDiscCisvV2, std::vector<Muon> const &selMuons, std::vector<Electron> const &selElectrons, const std::vector<TLorentzVectorWithIndex> & selPhotons);

  bool cleanPathologicEventsInPhotons(TString datasetName, float EvtRunNum, float EvtLumiNum, float EvtNum);
}

#endif
