#ifndef BTAGWEIGHT_H_
#define BTAGWEIGHT_H_

#include <memory>
#include <string>
#include <vector>

#include <TTreeReaderArray.h>

#include <BTagger.h>
#include <Options.h>
#include <PhysicsObjects.h>
#include <Tables.h>


class BTagCalibrationReader;


/**
 * \brief Computes weights for b tag scale factors
 */
class BTagWeight {
 public:
  /// Constructor from configuration options
  BTagWeight(Options const &options, BTagger const &bTagger);

  ~BTagWeight() noexcept;

  /**
   * \brief Computes weight for a simulated event with given jets
   *
   * \param[in] jets  Reconstructed jets in the event.
   */
  double operator()(std::vector<Jet> const &jets) const;

 private:
  /// Loads b tag efficiencies
  void LoadEffTables();

  /// Returns b tag efficiency computed for the jet with given properties
  double GetEfficiency(double pt, double eta, int flavour) const;

  /// Return b tag scale factor computed for the jet with given properties
  double GetScaleFactor(double pt, double eta, int flavour) const;

  /**
   * \brief Object that provides numeric value of the b tag discriminator and
   *   thresholds from configuration file
   */
  BTagger const &bTagger_;

  /// Path of b tag efficiencies tables
  std::string const effTablePath_;

  /// Tables with b tag efficiencies
  utils::tables efficiencyTables_;

  /// Object that provies values of b tag scale factors
  std::unique_ptr<BTagCalibrationReader> scaleFactorReader_;

  /// Requested systematic variation
  std::string syst_;
};

#endif  // BTAGWEIGHT_H_

