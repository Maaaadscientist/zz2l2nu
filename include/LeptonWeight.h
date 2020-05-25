#ifndef HZZ2L2NU_INCLUDE_LEPTONWEIGHT_H_
#define HZZ2L2NU_INCLUDE_LEPTONWEIGHT_H_

#include <WeightBase.h>

#include <Dataset.h>
#include <ElectronBuilder.h>
#include <MuonBuilder.h>
#include <Options.h>
#include <PhysicsObjects.h>


class PtEtaHistogram;


/**
 * \brief Applies lepton identification efficiency scale factors
 *
 * The computation is done using tight leptons provided by the builders given to
 * the constructor. The scale factors are read from files specified in the
 * master configuration. Multiple scale factors for the same lepton flavour are
 * multiplied together. See documentation for class PtEtaHistogram for a
 * description of the configuration for individual scale factors.
 *
 * Trigger scale factors are missing.
 *
 * Systematic variations are not provided.
 */
class LeptonWeight : public WeightBase {
 public:
  /// Constructor
  LeptonWeight(Dataset &dataset, Options const &options,
               ElectronBuilder const *electronBuilder,
               MuonBuilder const *muonBuilder);

  ~LeptonWeight();

  /// Computes the total lepton efficiency weight for the current event
  double NominalWeight() const override;

  /**
   * \brief Gives scale factor for an electron
   *
   * Uses pseudorapidity of the ECAL supercluster to look up the scale factor.
   */
  double ElectronSF(Electron const &electron) const;

  /**
   * \brief Gives scale factor for a muon
   *
   * Uses uncorrected momentum to look up the scale factor.
   */
  double MuonSF(Muon const &muon) const;

 private:
  /// Individual scale factors for muons and electrons
  std::vector<PtEtaHistogram> muonScaleFactors_, electronScaleFactors_;

  /// Non-owning pointer to object that provides collection of electrons
  ElectronBuilder const *electronBuilder_;

  /// Non-owning pointer to object that provides collection of muons
  MuonBuilder const *muonBuilder_;
};

#endif  // HZZ2L2NU_INCLUDE_LEPTONWEIGHT_H_

