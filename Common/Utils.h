#ifndef utils_h
#define utils_h

#include <iostream>
#include <string>
#include <TMath.h>
#include <map>
#include <string>
#include <TH1.h>
#include "TLorentzVectorWithIndex.h"

namespace utils
{
  double deltaPhi (TLorentzVector v1, TLorentzVector v2);

  double deltaPhi (float phi1, float phi2);

  double deltaR (TLorentzVector v1, TLorentzVector v2);
  
  double deltaR (float eta1, float phi1, float eta2, float phi2);

  double getPhotonEnergy (double pT, double eta);

  bool passVBFcuts (std::vector<TLorentzVectorWithIndex> selJets, TLorentzVector boson);

  double photon_rhoCorrectedIso(double pfIso, double rho, double sceta, TString isoType);

  double photonID_effArea(double sceta, TString isoType);

  bool passMetFilter(ULong64_t TrigMET, std::vector<std::pair<int, int> > & listMETFilter, bool isMC);

  bool file_exist(std::string name);

  std::map<double, double> TH1toMap(TH1D *h_weight);
  
  std::map<double, double> TH1toMap(std::string fileName, std::string histoName);

  namespace CutVersion { enum CutSet {Spring15Cut25ns, ICHEP16Cut, Moriond17Cut, Moriond17CutRunGH}; }
}

#endif
