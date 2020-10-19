#ifndef HZZ2L2NU_INCLUDE_GENWEIGHT_H_
#define HZZ2L2NU_INCLUDE_GENWEIGHT_H_

#include <array>
#include <map>
#include <optional>
#include <tuple>
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
    return availableVariations_.size();
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
   *   -# Factorization scale up and down,
   *   -# Combined PDF up and down,
   *   -# alpha_s in PDF up and down.
   */
  double RelWeight(int variation) const override;

  /**
   * \brief Returns relative weight for requested variation in alpha_s in PDF
   *
   * The variation is scaled linearly to the recommended size.
   */
  double RelWeightAlphaS(Var direction) const;

  /// Returns relative weight for requested variation in ME scales
  double RelWeightMEScale(Var renorm, Var factor) const;

  /**
   * \brief Returns relative weight for combined PDF variation.
   *
   * This is a pure PDF variation, i.e. it does not include alpha_s. Different
   * weights are combined on the fly on the per-even basis, as described
   * <a href="https://indico.cern.ch/event/961210/#199-pdf-variations">here</a>.
   * The resulting variation is conservative.
   */
  double RelWeightPdf(Var direction) const;

  std::string_view VariationName(int variation) const override {
    return availableVariations_.at(variation);
  }

private:
  /// Type of PDF variations
  enum class PdfVarType {
    Hessian,  ///< Hessian variations
    MC  ///< MC replicas
  };

  void InitializeLheScale(Dataset &dataset);
  void InitializePdf(Dataset &dataset);

  /// Helper function to look up metadata about a PDF set from YAML file
  static std::tuple<PdfVarType, std::optional<std::array<int, 2>>, double>
  LookUpPdfSet(int lhapdf);

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
   * \brief Mapping from variations in the renormalization and factorization
   * scales in the ME to the index of the corresponding weight in the source
   * array.
   */
  std::map<std::pair<Var, Var>, int> meScaleIndices_;

  /// Indicates whether PDF weights are available
  bool pdfWeightsPresent_;

  /**
   * \brief Indices for PDF weights with variations
   *
   * These indices correspond to the weights in branch LHEPdfWeight. They don't
   * include the nominal weight nor the alpha_s variations. The range of indices
   * to be used is [first, second), i.e. semiinclusive.
   */
  std::pair<int, int> pdfWeightsIndices_;

  /// Type of PDF variations
  PdfVarType pdfVarType_;

  /// Indicates whether PDF weights include also alpha_s variations
  bool alphaSWeightsPresent_;

  /**
   * \brief Indices of weights for alpha_s variation
   *
   * Correspond to the weights in branch LHEPdfWeight. The first index is for
   * an increase in alpha_s, the second one is for a decrease. The array is not
   * initialized if alphaSWeightsPresent_ is false.
   */
  std::array<int, 2> alphaSWeightsIndices_;

  /// How much to rescale alpha_s variation given by weights in input file
  double alphaSVarScaleFactor_;

  /// Names of available variations
  std::vector<std::string> availableVariations_;

  /**
   * \brief Default variation index to be used with \ref operator()
   *
   * Corresponds to the argument of RelWeight, with a special value of -1
   * denoting nominal weight.
   */
  int defaultVariationIndex_;

  mutable TTreeReaderValue<float> srcLheNominalWeight_;
  mutable TTreeReaderValue<float> srcGenNominalWeight_;
  mutable TTreeReaderArray<float> srcScaleWeights_;
  mutable TTreeReaderArray<float> srcPdfWeights_;
};

#endif  // HZZ2L2NU_INCLUDE_GENWEIGHT_H_
