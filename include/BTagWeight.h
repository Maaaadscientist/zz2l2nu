#ifndef HZZ2L2NU_INCLUDE_BTAGWEIGHT_H_
#define HZZ2L2NU_INCLUDE_BTAGWEIGHT_H_

#include <memory>
#include <string>

#include <TTreeReaderArray.h>

#include <BTagger.h>
#include <JetBuilder.h>
#include <Options.h>
#include <PhysicsObjects.h>
#include <Tables.h>


class BTagCalibrationReader;


/**
 * \brief Computes weights for b tag scale factors
 *
 * Jets considered in the computation are read from a JetBuilder.
 */
class BTagWeight {
 public:
  /**
   * \brief Construct
   *
   * \param[in] options     Configuration options
   * \param[in] bTagger     Object to tag jets
   * \param[in] jetBuilder  Object providing collection of jets
   *
   * Objects \c bTagger and \c jetBuilder must exist for the lifetime of this
   * BTagWeight.
   */
  BTagWeight(Options const &options, BTagger const *bTagger,
             JetBuilder const *jetBuilder);

  ~BTagWeight() noexcept;

  /// Computes weight for current simulated event
  double operator()() const;

 private:
  /// Loads b tag efficiencies
  void LoadEffTables();

  /// Returns b tag efficiency computed for the jet with given properties
  double GetEfficiency(double pt, double eta, int flavour) const;

  /// Return b tag scale factor computed for the jet with given properties
  double GetScaleFactor(double pt, double eta, int flavour) const;

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
  std::string syst_;
};

#endif  // HZZ2L2NU_INCLUDE_BTAGWEIGHT_H_

