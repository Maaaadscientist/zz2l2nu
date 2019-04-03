#ifndef EWCORRECTIONWEIGHT_H_
#define EWCORRECTIONWEIGHT_H_

#include <map>
#include <string>
#include <vector>
#include <utility>

#include <TLorentzVector.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <Options.h>


/**
 * \brief Applies precomputed higher-order electroweak corrections in the form
 *   of event weights
 * 
 * Tries to guess physics content of the input data set based on the name of its
 * catalog file. The correction is enabled for targeted processes. Otherwise a
 * weight of 1. is returned.
 */
class EWCorrectionWeight {
 public:
  /**
   * \brief Constructor
   *
   * \param[in] reader   TTreeReader object that provides read access to the
   *   current event in the input data set.
   * \param[in] options  Configuration options.
   */
  EWCorrectionWeight(TTreeReader &reader, Options const &options);

  /**
   * \brief Computes and returns the weight for the current event
   *
   * The weight is 1. if the correction is disabled.
   */
  double operator()() const;

 private:
  /// Reads correction table
  void readFile_and_loadEwkTable();
  
  /// Finds the right correction in the file
  std::vector<float> findCorrection(float sqrt_s_hat, float t_hat) const;
  
  std::map<std::string,std::pair<TLorentzVector,TLorentzVector>> reconstructGenLevelBosons() const;
  
  /// The main function, returns the kfactor
  double getEwkCorrections(std::map<std::string,std::pair<TLorentzVector,TLorentzVector>> genLevelLeptons, double & error) const;

  /**
   * \brief Path to the catalog file
   *
   * Used to guess the physics content of the input data set.
   */
  TString catalogPath_;

  /// Label of requested systematic variation
  std::string syst_;

  /// Indicates whether the reweighting is enabled
  bool enabled_;

  std::vector<std::vector<float>> ewTable_;

  mutable TTreeReaderArray<float> GLepBarePt, GLepBareEta, GLepBarePhi,
    GLepBareE;
  mutable TTreeReaderArray<int> GLepBareId, GLepBareSt, GLepBareMomId;
  mutable TTreeReaderArray<float> GPdfx1, GPdfx2;
  mutable TTreeReaderArray<int> GPdfId1, GPdfId2;
};

#endif  // EWCORRECTIONWEIGHT_H_

