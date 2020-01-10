#ifndef HZZ2L2NU_INCLUDE_KFACTORCORRECTION_H_
#define HZZ2L2NU_INCLUDE_KFACTORCORRECTION_H_

#include <WeightBase.h>

#include <filesystem>

#include <TGraph.h>
#include <TTreeReaderArray.h>

#include <Dataset.h>
#include <Options.h>

/**
 * \brief Computes k factor for gluon fusion production
 *
 * The k factor is computed if the per-dataset configuration contains parameter
 * "k_factor" and its value is "ggF". Otherwise a weight of 1 is returned for
 * every event.
 */
class KFactorCorrection : public WeightBase {
 public:
  KFactorCorrection(Dataset &dataset, Options const &options);

  /**
   * \brief Computes mass of generator-level Higgs boson
   *
   * The Higgs boson is reconstructed from generarol-level leptons.
   */
  double HiggsMass() const;

  /**
   * \brief Computes and returns the k factor for the current event
   *
   * The k factor is 1. if the correction is disabled.
   */
  virtual double NominalWeight() const override;

 private:
  bool enabled_;

  /// The k factor as a function of the mass of the Higgs boson
  std::unique_ptr<TGraph> kfactorGraph_;

  mutable TTreeReaderArray<float> genPartPt_, genPartEta_, genPartPhi_,
    genPartMass_;
  mutable TTreeReaderArray<int> genPartStatus_, genPartStatusFlags_;
};

#endif  // HZZ2L2NU_INCLUDE_KFACTORCORRECTION_H_

