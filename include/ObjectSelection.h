#ifndef objectselection_h
#define objectselection_h

#include <string>
#include <vector>

#include <TTreeReaderArray.h>

#include <PhysicsObjects.h>
#include <TLorentzVectorWithIndex.h>


namespace objectSelection
{

  bool selectJets(std::vector<TLorentzVectorWithIndex> & selJets, std::vector<TLorentzVectorWithIndex> & selCentralJets, std::vector<double> & btags, TTreeReaderArray<float> const &JetAk04Pt, TTreeReaderArray<float> const &JetAk04Eta, TTreeReaderArray<float> const &JetAk04Phi, TTreeReaderArray<float> const &JetAk04E, TTreeReaderArray<float> const &JetAk04Id, TTreeReaderArray<float> const &JetAk04NeutralEmFrac, TTreeReaderArray<float> const &JetAk04NeutralHadAndHfFrac, TTreeReaderArray<float> const &JetAk04NeutMult, TTreeReaderArray<float> const &JetAk04BDiscCisvV2, std::vector<Muon> const &selMuons, std::vector<Electron> const &selElectrons, std::vector<Photon> const &selPhotons);

  bool cleanPathologicEventsInPhotons(TString datasetName, float EvtRunNum, float EvtLumiNum, float EvtNum);
}

#endif
