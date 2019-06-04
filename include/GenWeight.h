#ifndef GENWEIGHT_H_
#define GENWEIGHT_H_

#include <map>
#include <utility>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>


/**
 * \brief Interface to access generator weights
 *
 * The nominal weight is provided with \ref operator(). Weights corresponding to
 * several systematic variations are provided with dedicated methods. These are
 * always relative weights, which should be applied in addition to the nominal
 * one.
 */
class GenWeight {
public:
  /// Supported directions for systematic variations
  enum class Var {
    Nominal,
    Up,
    Down
  };

  GenWeight(TTreeReader &reader);

  /**
   * \brief Finds the requested boundary of the envelope of variations in ME
   * scales
   *
   * The returned weight is relative with respect to the nominal one. All
   * variations of the two ME scales are considered, except for the cases when
   * they go in the opposite directions. The "up" boundary corresponds to the
   * largest relative weight (which becomes the smallest absolute weight if the
   * nominal weight is negative), and vice versa.
   *
   * Note that using such an envelope in a shape-based analysis leads to
   * incorrect correlations between different bins of the distribution.
   */
  double EnvelopeMEScale(Var direction) const;

  /// Returns nominal generator weight
  double operator()() const;

  /// Returns relative weight for requested variation in alpha_s in PDF
  double RelWeightAlphaS(Var direction) const;

  /// Returns relative weight for requested variation in ME scales
  double RelWeightMEScale(Var renorm, Var factor) const;

  /**
   * \brief Returns relative weight for given PDF replica
   *
   * It is expected but not checked that the provided replica index satisfies
   * 0 <= replica < 100.
   */
  double RelWeightPdf(int replica) const;

private:
  /**
   * \brief Mapping from variations in the renormalization and factorization
   * scales in the ME to the index of the corresponding weight in the source
   * array.
   */
  std::map<std::pair<Var, Var>, int> meScaleIndices_;

  /**
   * \brief Source array of event weights
   *
   * If the array is not empty, the first weight corresponds to the result of
   * GenEventInfoProduct::weight(), and it is followed by weights from
   * LHEEventProduct::weights().
   */
  mutable TTreeReaderArray<double> srcWeights_;
};

#endif  // GENWEIGHT_H_

