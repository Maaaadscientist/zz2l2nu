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
  double deltaPhi(TLorentzVector const &v1, TLorentzVector const &v2);

  double deltaPhi (float phi1, float phi2);

  double deltaR(TLorentzVector const &v1, TLorentzVector const &v2);
  
  double deltaR (float eta1, float phi1, float eta2, float phi2);

  double getPhotonEnergy (double pT, double eta);

  bool passVBFcuts(std::vector<TLorentzVectorWithIndex> const &selJets,
                   TLorentzVector const &boson);

  double photon_rhoCorrectedIso(double pfIso, double rho, double sceta,
                                TString const &isoType);

  double photonID_effArea(double sceta, TString const &isoType);

  bool passMetFilter(ULong64_t TrigMET, std::vector<std::pair<int, int> > & listMETFilter, bool isMC);

  bool file_exist(std::string const &name);

  std::map<double, double> TH1toMap(TH1D *h_weight);
  
  std::map<double, std::pair<double, double>> TH1toMap(
      std::string const &fileName, std::string const &histoName);

  void giveMassToPhoton(TLorentzVector & boson, TH1D *h_weight);

  void loadInstrMETWeights(
      bool weight_NVtx_exist, bool weight_Pt_exist, bool weight_Mass_exist,
      std::map<TString, std::map<double, std::pair<double, double>>> &NVtxWeight_map,
      std::map<TString, std::map<double, std::pair<double, double>>> &PtWeight_map,
      std::map<TString, TH1D*> &LineshapeMassWeight_map,
      std::string const &weightFileType, std::string const &base_path,
      std::vector<std::string> const &v_jetCat);

  double getTheoryUncertainties(std::vector<double> *EvtWeights,
                                TString const &syst);

  // DEPRECATED. This function is not used anymore, only kept for comparison
  // with this old method.
  double getPdfUncertainty(std::vector<double> const *EvtWeights, bool isUp);

  double getQCDScaleUncertainty(std::vector<double> const *EvtWeights,
                                bool isUp);

  double getAlphaUncertainty(std::vector<double> const *EvtWeights, bool isUp);

  namespace CutVersion { enum CutSet {Spring15Cut25ns, ICHEP16Cut, Moriond17Cut, Moriond17CutRunGH}; }
}

#endif
