#ifndef BTAGWEIGHT_H_
#define BTAGWEIGHT_H_

#include <string>
#include <vector>

#include <TTreeReaderArray.h>

#include <BTagCalibrationStandalone.h>
#include <Options.h>
#include <PhysicsObjects.h>
#include <Tables.h>


/**
 * \brief Computes weights for b tag scale factors
 */
class BTagWeight {
 public:
  /// Constructor from configuration options
  BTagWeight(Options const &options);

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

  /// Returns string representing given working point
  static std::string WPToText(BTagEntry::OperatingPoint wp);

  /// Chosen b-tagging algorithm
  std::string bTagAlgorithm_{"CSVv2"};

  /// Chosen working point for the b-tagging algorithm
  BTagEntry::OperatingPoint bTagWorkingPoint_{BTagEntry::OP_LOOSE};
  
  /**
   * \brief Numeric value of the b tag discriminator that corresponds to the
   *   chosen working point
   */
  double bTagCut_{0.5426};

  /// Tables with b tag efficiencies
  utils::tables efficiencyTables_;

  /// Object that provies values of b tag scale factors
  BTagCalibrationReader scaleFactorReader_;

  /// Requested systematic variation
  std::string syst_;
};

#endif  // BTAGWEIGHT_H_

