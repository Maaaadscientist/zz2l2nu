#ifndef HZZ2L2NU_INCLUDE_EWCORRECTIONWEIGHT_H_
#define HZZ2L2NU_INCLUDE_EWCORRECTIONWEIGHT_H_

#include <WeightBase.h>

#include <map>
#include <string>
#include <vector>
#include <utility>

#include <TLorentzVector.h>
#include <TString.h>
#include <TTreeReaderArray.h>

#include <Dataset.h>
#include <EventCache.h>
#include <Options.h>


/**
 * \brief Applies precomputed higher-order electroweak corrections in the form
 * of event weights
 *
 * The type of the applied correction is determined by parameter "ew_correction"
 * in per-dataset configuration. If this parameter is not found, no correction
 * is applied and a weight of 1 is returned for every event.
 */
class EWCorrectionWeight : public WeightBase {
 public:
  /**
   * \brief Constructor
   *
   * \param[in] dataset  Dataset that will be processed.
   * \param[in] options  Configuration options.
   */
  EWCorrectionWeight(Dataset &dataset, Options const &options);

  /**
   * \brief Returns the nominal weight for the current event
   *
   * The weight is 1. if the correction is disabled.
   */
  virtual double NominalWeight() const override;

  /**
   * \brief Returns the number of systematic variations
   *
   * If the correction is disabled for the current dataset, the number of
   * variations is 0.
   */
  virtual int NumVariations() const override;

  /**
   * \brief Returns weight for the systematic variation requested by command
   * line options
   *
   * If the correction is disabled, the weight is 1.
   */
  virtual double operator()() const override;

  virtual double RelWeight(int variation) const override;
  virtual std::string_view VariationName(int variation) const override;

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

  /// Updates cached \ref weightNominal_ and \ref weightError_;
  void Update() const;

  /// Selected type of EW correction
  Type correctionType_;

  EventCache cache_;

  /// Cached nominal event weight
  mutable double weightNominal_;

  /// Cached absolute error of the weight
  mutable double weightError_;

  /**
   * \brief Requested direction of systematic variation
   *
   * Allowed values are 0, +1, and -1.
   */
  int systDirection_;

  std::vector<std::vector<float>> ewTable_;

  mutable TTreeReaderArray<float> genPartPt_, genPartEta_, genPartPhi_,
    genPartMass_;
  mutable TTreeReaderArray<int> genPartPdgId_, genPartIdxMother_;
  mutable TTreeReaderValue<Float_t> generatorX1_, generatorX2_;
  mutable TTreeReaderValue<Int_t> generatorId1_, generatorId2_;
};

#endif  // HZZ2L2NU_INCLUDE_EWCORRECTIONWEIGHT_H_

