#ifndef PILEUPWEIGHT_H_
#define PILEUPWEIGHT_H_

#include <vector>


/**
 * \brief Computes the weight for difference in pileup profile
 *
 * This is a crude implementation with hard-coded weights and no systematic
 * uncertainties. To be redone.
 */
class PileUpWeight {
 public:
  /// Constructor
  PileUpWeight();

  /**
   * \brief Computes weight for the given expected number of pileup interactions
   */
  double operator()(double mu) const;

 private:
  /**
   * \brief Precomputed weights
   *
   * The index corresponds to the expected number of pileup interactions in a
   * simulated event, rounded to integers.
   */
  std::vector<double> weights;
};

#endif  // PILEUPWEIGHT_H_

