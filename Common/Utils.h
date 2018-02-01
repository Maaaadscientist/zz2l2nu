#ifndef utils_h
#define utils_h

#include <iostream>
#include <string>
#include <TMath.h>
#include "TLorentzVectorWithIndex.h"

namespace utils
{
  double deltaPhi (TLorentzVector v1, TLorentzVector v2);

  double deltaR (TLorentzVector v1, TLorentzVector v2);

  double getPhotonEnergy (double pT, double eta);

  bool passVBFcuts (std::vector<TLorentzVectorWithIndex> selJets, TLorentzVector boson);

  bool passMetFilter(ULong64_t TrigMET, std::vector<std::pair<int, int> > & listMETFilter, bool isMC);

  namespace CutVersion { enum CutSet {Spring15Cut25ns, ICHEP16Cut, Moriond17Cut}; }
}

#endif
