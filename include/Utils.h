#ifndef utils_h
#define utils_h

#include <cmath>
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TTreeReaderArray.h>
#include <TVector2.h>

#include <Options.h>
#include <PhysicsObjects.h>


class GenWeight;


namespace utils {

/// Computes squared distance in (eta, phi) metric
inline double DeltaR2(double eta1, double phi1, double eta2, double phi2) {
  double const dPhi = TVector2::Phi_mpi_pi(phi1 - phi2);
  return std::pow(eta1 - eta2, 2) + std::pow(dPhi, 2);
}

/// Computes squared distance in (eta, phi) metric
inline double DeltaR2(TLorentzVector const &p1, TLorentzVector const &p2) {
  return DeltaR2(p1.Eta(), p1.Phi(), p2.Eta(), p2.Phi());
}

double deltaPhi(TLorentzVector const &v1, TLorentzVector const &v2);

double deltaPhi (float phi1, float phi2);

bool PassVbfCuts(std::vector<Jet> const &jets, TLorentzVector const &boson);

std::map<double, double> TH1toMap(TH1D *h_weight);

std::map<double, std::pair<double, double>> TH1toMap(
    std::string const &fileName, std::string const &histoName);

std::map<std::pair<double, double>, std::pair<double, double>> TH2toMap(
    std::string const &fileName, std::string const &histoName);

void loadInstrMETWeights(
    bool applyNvtxWeights, bool applyPtWeights, bool applyMassLineshape,
    std::map<TString, std::map<double, std::pair<double, double>>> &NvtxWeight_map,
    std::map<TString, std::map<double, std::pair<double, double>>> &PtWeight_map,
    std::map<TString, TH1*> &LineshapeMassWeight_map,
    std::vector<std::string> const &v_jetCat,
    Options const &options);

void loadMeanWeights(
    bool applyMeanWeights,
    std::map<TString, std::map<double, double>> &meanWeight_map,
    std::vector<std::string> const &v_jetCat,
    Options const &options);

}  // namespace utils

#endif
