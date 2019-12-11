#ifndef HZZ2L2NU_INCLUDE_GENWEIGHT_H_
#define HZZ2L2NU_INCLUDE_GENWEIGHT_H_

#include <WeightBase.h>

#include <map>
#include <utility>

#include <TTreeReaderArray.h>

#include <Dataset.h>


/**
 * \brief Interface to access generator weights
 *
 * The nominal weight is provided with \ref operator(). It includes the
 * normalization to the cross section and the sum of the nominal weights in the
 * full dataset.
 *
 * Weights corresponding to several systematic variations are provided with
 * dedicated methods. These are always relative weights, which should be applied
 * in addition to the nominal one.
 */
class GenWeight : public WeightBase {
public:
  /// Supported directions for systematic variations
  enum class Var {
    Nominal,
    Up,
    Down
  };

  GenWeight(Dataset &dataset);

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

  /**
   * \brief Returns nominal generator weight
   *
   * It includes the normalization to the cross section and the sum of the
   * nominal weights in the full dataset.
   */
  virtual double NominalWeight() const override;

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
   * \brief Common dataset weight, in pb
   *
   * Computed as the cross section divided by the sum of the nominal event
   * weights in the full dataset.
   */
  double datasetWeight_;

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
  mutable TTreeReaderValue<float> srcLheNominalWeight_;
  mutable TTreeReaderValue<float> srcGenNominalWeight_;
  mutable TTreeReaderArray<float> srcPdfWeights_;
  mutable TTreeReaderArray<float> srcScaleWeights_;
  //mutable TTreeReaderArray<double> srcAlphaSWeights_; // Not there.
};

#endif  // HZZ2L2NU_INCLUDE_GENWEIGHT_H_

