#ifndef HZZ2L2NU_INCLUDE_PILEUPIDWEIGHT_H_
#define HZZ2L2NU_INCLUDE_PILEUPIDWEIGHT_H_

#include <WeightBase.h>

#include <array>
#include <memory>
#include <optional>

#include <TH2.h>
#include <TTreeReaderValue.h>

#include <Dataset.h>
#include <JetBuilder.h>
#include <Options.h>
#include <PhysicsObjects.h>
#include <PileUpIdFilter.h>
#include <XGBoostPredictor.h>


/**
 * \brief Computes event weights that account for pileup ID scale factors
 *
 * A multivariate parameterization is used for the efficiency of the pileup ID
 * in simulation.
 *
 * The behaviour is controlled by the section \c pileup_id of the master
 * configuration. Some information from it is extracted via PileUpIdFilter --
 * see its documentation on how the jet selection should be specified. In class
 * PileUpIdWeight it is assumed that section \c pileup_id exists; otherwise an
 * object of this class should not be constructed.
 *
 * No systematic uncertainties are provided.
 */
class PileUpIdWeight : public WeightBase {
 public:
  /**
   * \brief Constructor
   *
   * \param[in] dataset  Dataset that will be processed.
   * \param[in] options  Configuration options.
   * \param[in] pileUpIdFilter  Pointer to a valid PileUpIdFilter, which
   *   describes jet selection based on pileup ID.
   * \param[in] jetBuilder  Pointer to a valid JetBuilder.
   */
  PileUpIdWeight(
      Dataset &dataset, Options const &options,
      PileUpIdFilter const *pileUpIdFilter,
      JetBuilder const *jetBuilder);

  double NominalWeight() const override;

 private:
  /// Auxiliary structure that provides details specific to each |eta| range
  struct Context {
    Context(Jet::PileUpId wp)
        : workingPoint{wp} {}

    /// Working point for pileup ID used in this range
    Jet::PileUpId workingPoint;

    /// Histograms with pileup ID scale factors for matched and pileup jets
    std::shared_ptr<TH2> sfMatched, sfPileUp;
  };

  /**
   * \brief Finds which |eta| range the jet falls into and returns the
   * corresponding context object
   */
  Context const &FindContext(Jet const &jet) const;

  /// Computes pileup ID efficiency in simulation for given jet
  double GetEfficiency(Context const &context, Jet const &jet) const;

  /// Finds pileup ID scale factor for given jet
  double GetScaleFactor(Context const &context, Jet const &jet) const;

  /**
   * \brief Reads histograms with scale factors for all used working points
   *
   * \param[in] config  Master configuration.
   * \param[in] year  Year of data taking the scale factors for which should be
   *   read.
   */
  void LoadScaleFactors(YAML::Node const config, int year);

  PileUpIdFilter const *pileUpIdFilter_;
  JetBuilder const *jetBuilder_;

  /**
   * \brief Edges between different |eta| regions
   *
   * Copied from pileUpIdFilter_.
   */
  std::vector<double> absEtaEdges_;

  /// \ref Context for each |eta| region
  std::vector<Context> contexts_;

  /**
   * \brief Array with input features for the computation of the pileup ID
   * efficiency in simulation
   *
   * It will be reused for each computation. Features that don't change are set
   * in the constructor.
   */
  mutable std::array<float, 13> effFeatures_;

  /// Model that parameterizes the pileup ID efficiency in simulation
  std::optional<XGBoostPredictor> effCalc_;

  /// Interface to read the expected number of pileup interactions
  mutable TTreeReaderValue<float> expPileUp_;
};

#endif  // HZZ2L2NU_INCLUDE_PILEUPIDWEIGHT_H_

