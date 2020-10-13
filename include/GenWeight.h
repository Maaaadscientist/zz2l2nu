#ifndef HZZ2L2NU_INCLUDE_GENWEIGHT_H_
#define HZZ2L2NU_INCLUDE_GENWEIGHT_H_

#include <map>
#include <utility>

#include <TTreeReaderArray.h>

#include <Dataset.h>
#include <Options.h>
#include <WeightBase.h>


/**
 * \brief Interface to access generator weights
 *
 * The nominal weight includes the normalization to the cross section and the
 * sum of the nominal weights in the full dataset.
 *
 * Weights corresponding to several systematic variations are provided with
 * dedicated methods. These are always relative weights, which should be applied
 * in addition to the nominal one. The standard weight interface is also
 * implemented.
 */
class GenWeight : public WeightBase {
public:
  /// Supported directions for systematic variations
  enum class Var {
    Nominal,
    Up,
    Down
  };

  GenWeight(Dataset &dataset, Options const &options);

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
  double NominalWeight() const override;

  int NumVariations() const override {
    return (lheScaleWeightsPresent_) ? 4 : 0;
  }

  double operator()() const override {
    if (defaultVariationIndex_ == -1)
      return NominalWeight();
    else
      return NominalWeight() * RelWeight(defaultVariationIndex_);
  }

  /**
   * \brief Returns relative weight for requested systematic variation
   *
   * The order of the variations, provided that all are available, is as
   * follows:
   *   -# ME renormalization scale up and down,
   *   -# Factorization scale up and down.
   */
  double RelWeight(int variation) const override;

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

  std::string_view VariationName(int variation) const override;

private:
  void InitializeLheScale(Dataset &dataset);

  /**
   * \brief Common dataset weight, in pb
   *
   * Computed as the cross section divided by the sum of the nominal event
   * weights in the full dataset.
   */
  double datasetWeight_;

  /// Indicates whether LHE scale variations are available
  bool lheScaleWeightsPresent_;

 /**
  * \brief Indicates if the sample has bugged LHE weights (currently Pt-binned
  * WG)
  *
  * Problem seems to happen during translation from MiniAOD to NanoAOD: see
  * https://github.com/cms-nanoAOD/cmssw/issues/520
  */
  bool hasBuggedLheWeights_;

  /**
   * \brief Default variation index to be used with \ref operator()
   *
   * Corresponds to the argument of RelWeight, with a special value of -1
   * denoting nominal weight.
   */
  int defaultVariationIndex_;

  /**
   * \brief Mapping from variations in the renormalization and factorization
   * scales in the ME to the index of the corresponding weight in the source
   * array.
   */
  std::map<std::pair<Var, Var>, int> meScaleIndices_;

  mutable TTreeReaderValue<float> srcLheNominalWeight_;
  mutable TTreeReaderValue<float> srcGenNominalWeight_;
  mutable TTreeReaderArray<float> srcPdfWeights_;
  mutable TTreeReaderArray<float> srcScaleWeights_;
};

#endif  // HZZ2L2NU_INCLUDE_GENWEIGHT_H_
