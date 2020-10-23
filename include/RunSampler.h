#ifndef HZZ2L2NU_INCLUDE_RUNSAMPLER_H_
#define HZZ2L2NU_INCLUDE_RUNSAMPLER_H_

#include <optional>
#include <vector>
#include <utility>

#include <TTreeReaderValue.h>

#include <Dataset.h>
#include <EventCache.h>
#include <Options.h>
#include <TabulatedRandomGenerator.h>


/**
 * \brief Samples data-like run numbers in simulation
 *
 * When processing simulation, this class samples randomly run numbers from the
 * corresponding period in real data. The probability to sample each run number
 * is proportional to the luminosity recorded in it. When running on data, the
 * actual run number is returned. This allows accounting for time dependence in
 * simulation. For example, different corrections can be applied for different
 * run ranges.
 *
 * This class reads section \c run_sampler from the master configuration. The
 * following fields in this section are checked:
 * - \c luminosity  Path to YAML file with integrated luminosities for different
 *   runs. Resolved with FileInPath.
 * - \c range  Range of runs to filter data from the YAML file. Only runs within
 *     this range (including boundaries) will be selected. This field is
 *     optional. If missing, all runs from the YAML file will be used.
 */
class RunSampler {
 public:
  /// Type representing run number
  using run_t = int32_t;

  RunSampler(Dataset &dataset, Options const &options,
             TabulatedRngEngine &rngEngine);

  /**
   * \brief Returns the run number for the current event
   *
   * In simulation the returned run number is sampled randomly; in data this is
   * the actual run number of the current event. The result is cached per event,
   * and subsequent calls of this method for the same event will return the same
   * value.
   */
  run_t operator()() const;

 private:
  /// Updates per-event cache by sampling or reading new run number
  void Build() const;

  /**
   * \brief Loads integrated luminosities for different runs
   *
   * \param[in] config  YAML node "run_sampler" read from the master
   *   configuration.
   */
  void LoadData(YAML::Node const &config);

  /**
   * \brief Whether the run sampling is enabled
   *
   * If not, the run read from the input file is returned.
   */
  bool samplingEnabled_;

  /// Implements per-event caching
  EventCache cache_;

  /// Sampled or read run number cached for the current event
  mutable run_t currentRun_;

  /**
   * \brief Reader to access the run number from the input file
   *
   * Only used when \ref samplingEnabled_ is false.
   */
  mutable std::optional<TTreeReaderValue<UInt_t>> srcRun_;

  /// Random number generator to do the sampling
  TabulatedRandomGenerator tabulatedRng_;

  /**
   * \brief Cumulative probabilities for the sampling
   *
   * The first element in each pair is the cumulative probability for the run
   * given by the second element. The vector is ordered in the cumulative
   * probability. The cumulative probability is computed post-run, i.e. it's a
   * non-zero value for the first entry in the vector and 1 for the last entry.
   */
  std::vector<std::pair<double, run_t>> cumulProb_;
};

#endif  // HZZ2L2NU_INCLUDE_RUNSAMPLER_H_
