#ifndef HZZ2L2NU_INCLUDE_JETGEOMETRICVETO_H_
#define HZZ2L2NU_INCLUDE_JETGEOMETRICVETO_H_

#include <TTreeReaderValue.h>

#include <Dataset.h>
#include <JetBuilder.h>
#include <Options.h>
#include <TabulatedRandomGenerator.h>


/**
 * \brief Checks if event contains a jet in a specified (eta, phi) window
 *
 * The configuration is read from optional section \c jet_geometric_veto of the
 * master configuration. If this section is missing, the filtering is disabled,
 * and \ref operator() will always return true.
 *
 * In data the filtering is restricted to a specified run range. In simulation
 * the decision to apply the filtering or not is taken randomly, with the
 * probability given by the parameter \c lumi_fraction.
 */
class JetGeometricVeto {
 public:
  JetGeometricVeto(
      Dataset &dataset, Options const &options, JetBuilder const *jetBuilder,
      TabulatedRngEngine &rngEngine);

  /**
   * \brief Returns true if there are no jets in the specified window
   *
   * The event is considered good if there are no jets in the window.
   */
  bool operator()() const;

 private:
  /**
   * \brief Range of runs for which the veto should be evaluated
   *
   * Both boundaries are included in the range. Only used for real data.
   */
  int minRun_, maxRun_;

  /**
   * \brief Probability with which to perform the check in simulation
   *
   * This should correspond to the fraction of the integrated luminosity
   * contained within the run range in data given by \ref minRun_ and
   * \ref maxRun_.
   */
  double lumiFraction_;

  /// Window in (eta, phi) to be vetoed
  double minEta_, maxEta_, minPhi_, maxPhi_;

  bool enabled_;
  bool isSim_;
  JetBuilder const *jetBuilder_;
  TabulatedRandomGenerator tabulatedRng_;
  mutable TTreeReaderValue<UInt_t> srcRun_;
};

#endif  // HZZ2L2NU_INCLUDE_JETGEOMETRICVETO_H_
