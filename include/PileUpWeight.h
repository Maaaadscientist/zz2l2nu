#ifndef HZZ2L2NU_INCLUDE_PILEUPWEIGHT_H_
#define HZZ2L2NU_INCLUDE_PILEUPWEIGHT_H_

#include <WeightBase.h>

#include <array>
#include <filesystem>
#include <memory>
#include <string>

#include <TH1.h>
#include <TTreeReaderValue.h>

#include <Dataset.h>
#include <EventCache.h>
#include <Options.h>


/**
 * \brief Performs reweighting for pileup profile
 *
 * Paths to ROOT files containing the pileup profile used to produce simulation
 * and the target profile in data are read from the analysis configuration file.
 * Systematic uncertainty in pileup is implemented using alternative data
 * profiles, which have been constructed assuming different values for the
 * pileup cross section. Computed weights are cached on per-event basis.
 */
class PileUpWeight : public WeightBase {
 public:
  /// Constructor
  PileUpWeight(Dataset &dataset, Options const &options);

  virtual double NominalWeight() const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[0];
  }

  virtual int NumVariations() const override {
    return 2;
  }

  /// Computes the pileup weight for the current event
  virtual double operator()() const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[defaultWeightIndex_];
  }

  virtual double RelWeight(int variation) const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[variation + 1] / weights_[0];
  }

  virtual std::string_view VariationName(int variation) const override;

 private:
  /**
   * \brief Reads a histogram with given name from a ROOT file
   *
   * Checks for and reports errors. The returned histogram is owned by the
   * caller.
   */
  static TH1 *ReadHistogram(std::filesystem::path const &path,
                            std::string const &name);

  /// Computes all weights for the current event
  void Update() const;

  ///@{
  /**
   * \brief Histograms representing pileup profiles in data and simulation
   *
   * They are normalized to represent probability density. Under- and overflow
   * bins are empty. The profiles in data are given in the order nominal, up,
   * down.
   */
  std::array<std::unique_ptr<TH1>, 3> dataProfiles_;
  std::unique_ptr<TH1> simProfile_;
  ///@}

  EventCache cache_;

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
  
  /**
   * \brief Interface to read the expected number of pileup interactions
   */
  mutable TTreeReaderValue<float> mu_;
};

#endif  // HZZ2L2NU_INCLUDE_PILEUPWEIGHT_H_

