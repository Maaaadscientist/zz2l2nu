#ifndef HZZ2L2NU_INCLUDE_BTAGWEIGHT_H_
#define HZZ2L2NU_INCLUDE_BTAGWEIGHT_H_

#include <WeightBase.h>

#include <array>
#include <memory>
#include <string>

#include <TTreeReaderArray.h>

#include <BTagger.h>
#include <Dataset.h>
#include <EventCache.h>
#include <JetBuilder.h>
#include <Options.h>
#include <PhysicsObjects.h>
#include <Tables.h>


class BTagCalibrationReader;


/**
 * \brief Computes weights for b tag scale factors
 *
 * Jets considered in the computation are read from a JetBuilder. Computed
 * weights are cached on per-event basis.
 */
class BTagWeight : public WeightBase {
 public:
  /**
   * \brief Construct
   *
   * \param[in] dataset     Current dataset
   * \param[in] options     Configuration options
   * \param[in] bTagger     Object to tag jets
   * \param[in] jetBuilder  Object providing collection of jets
   *
   * Objects \c bTagger and \c jetBuilder must exist for the lifetime of this
   * BTagWeight.
   */
  BTagWeight(Dataset &dataset, Options const &options, BTagger const *bTagger,
             JetBuilder const *jetBuilder);

  ~BTagWeight() noexcept;

  virtual double NominalWeight() const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[0];
  }

  virtual int NumVariations() const override {
    return 4;
  }

  virtual double operator()() const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[int(defaultVariation_)];
  }

  virtual double RelWeight(int variation) const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[variation + 1] / weights_[0];
  }

  virtual std::string_view VariationName(int variation) const override;

 private:
  /**
   * \brief Supported systematic variations
   *
   * Cast to int, a variable of this type can index array \ref weights_.
   */
  enum class Variation : int {
    kNominal = 0,
    kTagUp = 1,
    kTagDown = 2,
    kMistagUp = 3,
    kMistagDown = 4
  };

  /// Computes event weight for the given systematic variation
  double ComputeWeight(Variation variation) const;

  /// Loads b tag efficiencies
  void LoadEffTables();

  /// Returns b tag efficiency computed for the jet with given properties
  double GetEfficiency(double pt, double eta, int flavour) const;

  /// Return b tag scale factor computed for the jet with given properties
  double GetScaleFactor(double pt, double eta, int flavour,
                        Variation variation) const;

  /// Computes event weights for all systematic variations
  void Update() const;

  /**
   * \brief Non-owning pointer to object that provides numeric value of the
   * b tag discriminator and thresholds from configuration file
   */
  BTagger const *bTagger_;

  /// Non-owning pointer to JetBuilder that provides jets
  JetBuilder const *jetBuilder_;

  /// Path of b tag efficiencies tables
  std::string const effTablePath_;

  /// Tables with b tag efficiencies
  utils::tables efficiencyTables_;

  /// Object that provies values of b tag scale factors
  std::unique_ptr<BTagCalibrationReader> scaleFactorReader_;

  /// Requested systematic variation
  Variation defaultVariation_;

  EventCache cache_;

  /**
   * \brief Cached weights for all systematic variations
   *
   * The order is the same as in enum \ref Variation.
   */
  mutable std::array<double, 5> weights_;
};

#endif  // HZZ2L2NU_INCLUDE_BTAGWEIGHT_H_

