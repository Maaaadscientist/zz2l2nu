#ifndef HZZ2L2NU_INCLUDE_PILEUPWEIGHT_H_
#define HZZ2L2NU_INCLUDE_PILEUPWEIGHT_H_

#include <array>
#include <filesystem>
#include <memory>
#include <string>

#include <TH1.h>
#include <TTreeReaderValue.h>

#include <Dataset.h>
#include <EventCache.h>
#include <Options.h>
#include <RunSampler.h>
#include <WeightBase.h>


/**
 * \brief Performs reweighting for pileup profile
 *
 * Paths to ROOT files containing the pileup profile used to produce simulation
 * and the target profiles in data are read from the analysis configuration
 * file. Systematic uncertainty in pileup is implemented using alternative data
 * profiles, which have been constructed assuming different values for the
 * pileup cross section.
 *
 * The file with target pileup profiles must contain directories for one or more
 * data taking eras. In each directory the following entries are expected:
 * - \c run_range  Range of runs included in this era. Represented by
 *   \c TVectorD of size 2. The boundaries of the range are included.
 * - \c nominal, \c up, \c down  \c TH1 representing nominal pileup profile and
 *   profiles for systamatic varations. Expect an integer binning.
 *
 * Computed weights are cached on per-event basis.
 */
class PileUpWeight : public WeightBase {
 public:
  /// Constructor
  PileUpWeight(Dataset &dataset, Options const &options,
               RunSampler const *runSampler);

  double NominalWeight() const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[0];
  }

  int NumVariations() const override {
    return 2;
  }

  /// Computes the pileup weight for the current event
  double operator()() const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[defaultWeightIndex_];
  }

  double RelWeight(int variation) const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[variation + 1] / weights_[0];
  }

  std::string_view VariationName(int variation) const override;

 private:
  /// Data pileup profiles for a single era
  struct Era {
    /// Range of runs defining the era (boundaries included)
    int minRun, maxRun;

    /**
     * \brief Data pileup profiles
     *
     * The nominal profile and the profiles for up and down systematic
     * varations, in that order. All profiles are normalized to represent
     * probability density. Under- and overflow bins in the histograms are
     * empty.
     */
    std::array<std::unique_ptr<TH1>, 3> dataProfiles;
  };

  /// Loads data pileup profiles for all eras
  void LoadDataProfiles(std::filesystem::path const &path);

  /// Loads pileup profile in simulation for a given dataset or default one
  void LoadSimProfile(YAML::Node const &config, std::string const &datasetName);

  /// Computes all weights for the current event
  void Update() const;

  /**
   * \brief Histogram representing pileup profile in simulation
   *
   * The histogram is normalized to represent probability density. Under- and
   * overflow bins are empty.
   */
  std::unique_ptr<TH1> simProfile_;

  /**
   * \brief Eras with pileup profiles in data
   *
   * The collection is sorted by Era::minRun. There assumed to be no overlap
   * between different eras.
   */
  std::vector<Era> eras_;

  EventCache cache_;

  /// Non-owning pointer to an object that samples representative run numbers
  RunSampler const *runSampler_;

  /**
   * \brief Cached weights
   *
   * They are given in the order nominal, up, down.
   */
  mutable std::array<double, 3> weights_;

  /**
   * \brief Index of the default weight to be returned by operator()
   *
   * Corresponds to array \ref weights_.
   */
  int defaultWeightIndex_;

  /// Interface to read the expected number of pileup interactions
  mutable TTreeReaderValue<float> mu_;
};

#endif  // HZZ2L2NU_INCLUDE_PILEUPWEIGHT_H_
