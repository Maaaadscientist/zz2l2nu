#ifndef EWCORRECTIONWEIGHT_H_
#define EWCORRECTIONWEIGHT_H_

#include <map>
#include <string>
#include <vector>
#include <utility>

#include <TLorentzVector.h>
#include <TString.h>
#include <TTreeReaderArray.h>

#include <Dataset.h>
#include <Options.h>


/**
 * \brief Applies precomputed higher-order electroweak corrections in the form
 * of event weights
 *
 * The type of the applied correction is determined by parameter "ew_correction"
 * in per-dataset configuration. If this parameter is not found, no correction
 * is applied and a weight of 1 is returned for every event.
 */
class EWCorrectionWeight {
 public:
  /**
   * \brief Constructor
   *
   * \param[in] dataset  Dataset that will be processed.
   * \param[in] options  Configuration options.
   */
  EWCorrectionWeight(Dataset &dataset, Options const &options);

  /**
   * \brief Computes and returns the weight for the current event
   *
   * The weight is 1. if the correction is disabled.
   */
  double operator()() const;

 private:
  /// Supported types of EW corrections
  enum class Type {
    None,
    ZZ,
    WZ
  };

  /// Reads correction table
  void readFile_and_loadEwkTable();
  
  /// Finds the right correction in the file
  std::vector<float> findCorrection(float sqrt_s_hat, float t_hat) const;
  
  std::map<std::string,std::pair<TLorentzVector,TLorentzVector>> reconstructGenLevelBosons() const;
  
  /// The main function, returns the kfactor
  double getEwkCorrections(std::map<std::string,std::pair<TLorentzVector,TLorentzVector>> genLevelLeptons, double & error) const;

  /// Selected type of EW correction
  Type correctionType_;

  /// Label of requested systematic variation
  std::string syst_;

  std::vector<std::vector<float>> ewTable_;

  mutable TTreeReaderArray<float> GLepBarePt, GLepBareEta, GLepBarePhi,
    GLepBareE;
  mutable TTreeReaderArray<int> GLepBareId, GLepBareSt, GLepBareMomId;
  mutable TTreeReaderArray<float> GPdfx1, GPdfx2;
  mutable TTreeReaderArray<int> GPdfId1, GPdfId2;
};

#endif  // EWCORRECTIONWEIGHT_H_

