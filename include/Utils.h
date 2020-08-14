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

double photon_rhoCorrectedIso(double pfIso, double rho, double sceta,
                              TString const &isoType);

double photonID_effArea(double sceta, TString const &isoType);

bool file_exist(std::string const &name);

std::map<double, double> TH1toMap(TH1D *h_weight);

std::map<double, std::pair<double, double>> TH1toMap(
    std::string const &fileName, std::string const &histoName);

std::map<std::pair<double, double>, std::pair<double, double>> TH2toMap(
    std::string const &fileName, std::string const &histoName);

void giveMassToPhoton(TLorentzVector & boson, TH1 *h_weight);

void loadInstrMETWeights(
    bool weight_NVtx_exist, bool weight_Pt_exist, bool weight_Mass_exist,
    std::string const &nvtxFile, std::string const &ptFile, std::string const &massFile,
    std::map<TString, std::map<std::pair<double,double>, std::pair<double, double>>> &NVtxWeight_map,
    std::map<TString, std::map<double, std::pair<double, double>>> &PtWeight_map,
    std::map<TString, TH1*> &LineshapeMassWeight_map,
    std::vector<std::string> const &v_jetCat);

void loadMeanWeights(
    bool meanWeight_exist, std::string const &meanWeightsFile,
    std::map<TString, std::map<double, double>> &meanWeight_map,
    std::vector<std::string> const &v_jetCat);

/**
 * \brief Forwards computation of systematic variations in generator weights to
 * class GenWeight
 */
double getTheoryUncertainties(GenWeight const &genWeight,
                              std::string_view syst);

namespace CutVersion { enum CutSet {Spring15Cut25ns, ICHEP16Cut, Moriond17Cut, Moriond17CutRunGH}; }

}  // namespace utils

#endif
