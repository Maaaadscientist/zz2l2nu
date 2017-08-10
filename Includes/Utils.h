#ifndef utils_h
#define utils_h

#include <iostream>
#include <string>
#include <TMath.h>
#include "TLorentzVector.h"

namespace utils
{
  double deltaPhi (TLorentzVector v1, TLorentzVector v2);

  double deltaR (TLorentzVector v1, TLorentzVector v2);

  double getPhotonEnergy (double pT, double eta);

  bool passVBFcuts (std::vector<TLorentzVector> selJets, TLorentzVector boson);
}

#endif
