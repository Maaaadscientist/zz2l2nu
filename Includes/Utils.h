#ifndef utils_h
#define utils_h

#include <iostream>
#include <string>
#include <TMath.h>
#include "TLorentzVector.h"

namespace utils
{
  enum {DoubleMu, SingleMu, DoubleE, HighPtE, SingleE, EMu, SinglePhoton, Ntrig}; //List of triggers used for our analysis

  double deltaPhi (TLorentzVector v1, TLorentzVector v2);

  double deltaR (TLorentzVector v1, TLorentzVector v2);

  double getPhotonEnergy (double pT, double eta);

  bool passVBFcuts (std::vector<TLorentzVector> selJets, TLorentzVector boson);

  int passTrigger(int trig, ULong64_t TrigHltDiMu, ULong64_t TrigHltMu, ULong64_t TrigHltDiEl, ULong64_t TrigHltEl, ULong64_t TrigHltElMu, ULong64_t TrigHltPhot, std::vector<unsigned int> *TrigHltDiMu_prescale, std::vector<unsigned int> *TrigHltMu_prescale, std::vector<unsigned int> *TrigHltDiEl_prescale, std::vector<unsigned int> *TrigHltEl_prescale, std::vector<unsigned int> *TrigHltElMu_prescale, std::vector<unsigned int> *TrigHltPhot_prescale, double selectedPhotonPt = 0);

  bool passMetFilter(ULong64_t TrigMET, std::vector<std::pair<int, int> > & listMETFilter, bool isMC);

  namespace CutVersion { enum CutSet {Spring15Cut25ns, ICHEP16Cut, Moriond17Cut}; }
}

#endif
