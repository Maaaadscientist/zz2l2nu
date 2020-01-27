#ifndef HZZ2L2NU_INCLUDE_JETCORRECTOR_H_
#define HZZ2L2NU_INCLUDE_JETCORRECTOR_H_

#include <array>
#include <filesystem>
#include <memory>
#include <vector>

#include <TLorentzVector.h>
#include <TTreeReaderValue.h>

#include <Dataset.h>
#include <Options.h>
#include <PhysicsObjects.h>
#include <TabulatedRandomGenerator.h>


// Classes from the JME POG that implement jet corrections are hidden from user
// as in the Pimpl idiom
class FactorizedJetCorrector;
class JetCorrectionUncertainty;

namespace JME {
class JetResolution;
class JetResolutionScaleFactor;
}


/**
 * \brief Implements jet pt scale and resolution corrections
 *
 * This class computes JEC (full nominal and L1-only), JEC uncertainty, and JER
 * smearing factors. Systematic variations are provided if requested via the
 * \c syst option.
 *
 * Paths to files that define JEC and JER are read from sections
 * \c jets/corrections and \c jets/resolution of the master configuration.
 */
class JetCorrector {
 public:
  JetCorrector(Dataset &dataset, Options const &options,
               TabulatedRngEngine &rngEngine);
  ~JetCorrector() noexcept;

  /**
   * \brief Computes JEC with all levels specified in the configuration applied
   *
   * In every event, \ref UpdateIov must be called before the first call to this
   * method. The effect of JEC uncertainties is not included.
   */
  double GetJecFull(TLorentzVector const &rawP4, double area) const;

  /**
   * \brief Computes L1 JEC
   *
   * In every event, \ref UpdateIov must be called before the first call to this
   * method.
   */
  double GetJecL1(TLorentzVector const &rawP4, double area) const;

  /**
   * \brief Computes correction factor to account for JEC uncertainty
   *
   * \param[in] corrP4  Corrected four-momentum of a jet.
   * \return Correction factor to rescale jet four-momentum to reproduce
   *   requested systematic variation in JEC, or 1 if no such variation has been
   *   requested.
   */
  double GetJecUncFactor(TLorentzVector const &corrP4) const;
  
  /**
   * \brief Computes correction factor to account for JER smearing
   *
   * \param[in] corrP4        Corrected four-momentum of a jet.
   * \param[in] genJet        Non-owning pointer to the generator-level jet
   *   matched to the given reconstructed jet. Can be nullptr if there is no
   *   match or when the stochastic smearing is desired.
   * \param[in] ptResolution  Relative pt resolution as computed by
   *   \ref GetPtResolution.
   * \param[in] rngChannel    Channel to be used for the tabulated random number
   *   generator.
   * \return Correction factor to rescale jet four-momentum.
   *
   * The input four-momentum must have JEC applied. Normally, only the nominal
   * JEC should be applied, even when a JEC variation has been requested. This
   * is consistent with how JER smearing is applied in CMSSW.
   *
   * The returned correction factor accounts for a systematic shift in JER if
   * requested in the options. This method should only be called for simulation.
   */
  double GetJerFactor(TLorentzVector const &corrP4, GenJet const *genJet,
                      double ptResultion, int rngChannel) const;

  /**
   * \brief Returns relative jet pt resolution in simulation
   *
   * \param[in] corrP4  Corrected four-momentum of a jet.
   * \return Relative pt resolution.
   */
  double GetPtResolution(TLorentzVector const &corrP4) const;

  /**
   * \brief If needed, update IOV based on the current run
   *
   * Must be called in every event before the first call to \ref GetJecFull or
   * \ref GetJecL1.
   */
  void UpdateIov() const;

 private:
  /// Supported types of systematic variations
  enum class Syst {
    None,
    JEC,
    JER
  };

  /// Possible directions for systematic variations
  enum class SystDirection {
    Up,
    Down
  };

  /// Type for run number
  using run_t = uint64_t;

  /// Auxiliary structure aggregating parameters that depent on IOV
  struct IovParams {
    /// Checks if the given run is included in this IOV
    bool Contains(run_t run) const {
      return run >= runRange[0] and run <= runRange[1];
    }

    /**
     * \brief Run range that defines the IOV
     *
     * Endpoints are included.
     */
    std::array<run_t, 2> runRange;

    /// Fully qualified paths to files defining JEC
    std::vector<std::filesystem::path> jecLevels;
  };

  /// Constructs an object to apply JEC using parameters of current IOV
  void LoadJec() const;

  /**
   * \brief Reads IOV parameters from the given node of the master configuration
   * file
   *
   * The given node must be a sequence of one or more nodes, each of which
   * contains keys \c run_range and \c levels. The node \c run_range is a
   * sequence of two run numbers that define the IOV; it is optional, and if
   * missing, a catch-all IOV is constructed. The node \c levels is a sequence
   * of paths to files defining different JEC levels.
   */
  void ReadIovParams(YAML::Node const config);

  /// Type of requested systematic variation
  Syst syst_;

  /// Direction of requested systematic variation
  SystDirection systDirection_;

  /// Registered IOV-dependent parameters
  std::vector<IovParams> iovs_;

  /// Currently active element of \ref iovs_
  mutable IovParams const *currentIov_;

  /// Run seen most recently
  mutable run_t cachedRun_;

  /**
   * \brief Object to compute JEC
   *
   * It is constructed using run-dependent parameters, which are provided via
   * the IOVs.
   */
  mutable std::unique_ptr<FactorizedJetCorrector> jetEnergyCorrector_;

  /**
   * \brief Object to provide JEC uncertainty
   *
   * Only created when a systematic variation in JEC has been requested and only
   * when processing simulation.
   */
  std::unique_ptr<JetCorrectionUncertainty> jecUncProvider_;

  /**
   * \brief Object that provides pt resolution in simulation
   *
   * Only created when running on simulation.
   */
  std::unique_ptr<JME::JetResolution> jerProvider_;

  /**
   * \brief Object that provides data-to-simulation scale factors for jet pt
   * resolution
   *
   * Only created when running on simulation.
   */
  std::unique_ptr<JME::JetResolutionScaleFactor> jerSFProvider_;

  /// Random number generator
  TabulatedRandomGenerator tabulatedRng_;

  /// Reader to access the current run
  mutable TTreeReaderValue<UInt_t> run_;

  /// Median angular pt density
  mutable TTreeReaderValue<float> rho_;
};

#endif  // HZZ2L2NU_INCLUDE_JETCORRECTOR_H_

