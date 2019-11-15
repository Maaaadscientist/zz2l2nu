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

class FactorizedJetCorrector;


/**
 * \brief Computes JEC
 *
 * Paths to files defining JEC are read from section jets/corrections of the
 * master configuration file.
 */
class JetCorrector {
 public:
  JetCorrector(Dataset &dataset, Options const &options);
  ~JetCorrector() noexcept;

  /**
   * \brief Computes JEC with all levels specified in the configuration applied
   *
   * In every event, \ref UpdateIov must be called before the first call to this
   * method.
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
   * \brief If needed, update IOV based on the current run
   *
   * Must be called in every event before the first call to \ref GetJecFull or
   * \ref GetJecL1.
   */
  void UpdateIov() const;

 private:
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

  /// Reader to access the current run
  mutable TTreeReaderValue<UInt_t> run_;

  /// Median angular pt density
  mutable TTreeReaderValue<float> rho_;

  /// Registered IOV-dependent parameters
  std::vector<IovParams> iovs_;

  /// Currently active element of \ref iovs_
  mutable IovParams const *currentIov_;

  /// Run seen most recently
  mutable run_t cachedRun_;

  /// Object to compute JEC
  mutable std::unique_ptr<FactorizedJetCorrector> jetEnergyCorrector_;
};

#endif  // HZZ2L2NU_INCLUDE_JETCORRECTOR_H_

