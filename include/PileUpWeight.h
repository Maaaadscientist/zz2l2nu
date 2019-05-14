#ifndef PILEUPWEIGHT_H_
#define PILEUPWEIGHT_H_

#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderValue.h>


/**
 * \brief Computes the weight for difference in pileup profile
 *
 * This is a crude implementation with hard-coded weights and no systematic
 * uncertainties. To be redone.
 */
class PileUpWeight {
 public:
  /// Constructor
  PileUpWeight(TTreeReader &reader);

  /**
   * \brief Computes the weight for the current event
   */
  double operator()() const;

 private:
  /**
   * \brief Precomputed weights
   *
   * The index corresponds to the expected number of pileup interactions in a
   * simulated event, rounded to integers.
   */
  std::vector<double> weights;
  
  /**
   * \brief Interface to read the expected number of pileup interactions
   *
   * Although in reality this is a floating-point number, in baobabs it is
   * silently converted to an integer [1]. When doing so, the fractional part is
   * discarded.
   *
   * [1] https://gitlab.cern.ch/HZZ-IIHE/shears/blob/a914b9448a769bce9ece367b10074ea6b721583b/Baobabs/src/Tupel.cc#L1212
   */
  mutable TTreeReaderValue<Int_t> mu_;
};

#endif  // PILEUPWEIGHT_H_

