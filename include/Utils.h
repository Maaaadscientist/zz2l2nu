#ifndef HZZ2L2NU_INCLUDE_UTILS_H_
#define HZZ2L2NU_INCLUDE_UTILS_H_

#include <cmath>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TTreeReaderArray.h>
#include <TVector2.h>

#include <HZZException.h>
#include <Options.h>
#include <PhysicsObjects.h>


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
    bool applyNvtxWeights, bool applyEtaWeights_, bool applyPtWeights,
    bool applyMassLineshape,
    std::map<TString, std::map<double, std::pair<double, double>>> &NvtxWeight_map,
    std::map<TString, std::map<double, std::pair<double, double>>> &EtaWeight_map,
    std::map<TString, std::map<double, std::pair<double, double>>> &PtWeight_map,
    std::map<TString, TH1*> &LineshapeMassWeight_map,
    std::vector<std::string> const &v_jetCat,
    Options const &options);

void loadMeanWeights(
    bool applyMeanWeights,
    std::map<TString, std::map<double, double>> &meanWeight_map,
    std::vector<std::string> const &v_jetCat,
    Options const &options);

/**
 * \brief Reads a histogram with given name from a ROOT file
 *
 * Checks for and reports errors. The check for missing historgram can be
 * disabled using the last argument. The returned histogram is owned by the
 * caller.
 */
template<typename T = TH1>
std::unique_ptr<T> ReadHistogram(
    std::filesystem::path const &path, std::string const &name,
    bool checkMissing = true) {
  TFile inputFile{path.c_str()};
  if (inputFile.IsZombie()) {
    HZZException exception;
    exception << "Could not open file " << path << ".";
    throw exception;
  }

  auto hist = dynamic_cast<T *>(inputFile.Get(name.c_str()));
  if (hist)
    hist->SetDirectory(nullptr);
  inputFile.Close();

  if (checkMissing and not hist) {
    HZZException exception;
    exception << "File " << path << " does not contain required histogram \"" <<
      name << "\".";
    throw exception;
  }

  return std::unique_ptr<T>{hist};
}

}  // namespace utils

#endif  // HZZ2L2NU_INCLUDE_UTILS_H_
