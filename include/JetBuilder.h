#ifndef JETBUILDER_H_
#define JETBUILDER_H_

#include <initializer_list>
#include <memory>
#include <vector>

#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <CollectionBuilder.h>
#include <Dataset.h>
#include <GenJetBuilder.h>
#include <JetCorrector.h>
#include <Options.h>
#include <PhysicsObjects.h>
#include <TabulatedRandomGenerator.h>


// Classes from the JME POG that implement jet corrections are hidden from user
// as in the Pimpl idiom
class JetCorrectionUncertainty;

namespace JME {
class JetResolution;
class JetResolutionScaleFactor;
};


/**
 * \brief Lazily builds collection of reconstructed jets
 *
 * Reads parameters from section "jets" of the master configuration. When
 * running over simulation, JER smearing is applied. To follow the standard
 * smearing algorithm, this builder needs to be made aware of generator-level
 * jets via method \ref SetGenJetBuilder. Systematic variations in JEC and JER
 * are supported. The changes in momenta of jets with pt > 15 GeV are aggregated
 * for \ref GetSumMomentumShift.
 */
class JetBuilder : public CollectionBuilder<Jet> {
 public:
  JetBuilder(Dataset &dataset, Options const &options,
             TabulatedRngEngine &rngEngine);
  ~JetBuilder() noexcept;

  /// Returns collection of jets
  std::vector<Jet> const &Get() const override;

  /**
   * \brief Specifies an object that provides generator-level jets
   *
   * Generator-level jets are used for JER smearing. However, if no such object
   * is given, stochastic version of the smearing will be performed. When given,
   * the object must have an approriate life time.
   */
  void SetGenJetBuilder(GenJetBuilder const *genJetBuilder);

 private:
  /// Supported types of systematic variations
  enum class Syst {
    None,
    JEC,
    JER
  };

  /// Possible directions for systematic variations
  enum class SystDirection {
    Up,
    Down
  };

  /// Constructs jets in the current event
  void Build() const override;

  /**
   * \brief Computes correction factor to account for JEC uncertainty
   *
   * \param[in] corrP4  Corrected four-momentum of a jet.
   * \return Correction factor to rescale jet four-momentum to reproduce
   *   requested systematic variation in JEC, or 1 if no such variation has been
   *   requested.
   */
  double ComputeJecUncFactor(TLorentzVector const &corrP4) const;

  /**
   * \brief Computes correction factor to account for JER smearing
   *
   * \param[in] corrP4      Corrected four-momentum of a jet.
   * \param[in] rho         Pileup rho in the current event.
   * \param[in] rngChannel  Channel to be used for the tabulated random number
   *   generator.
   * \return Correction factor to rescale jet four-momentum.
   *
   * The input four-momentum must have JEC applied. Normally, only the nominal
   * JEC should be applied, even when a JEC variation has been requested. This
   * is consistent with how JER smearing is applied in CMSSW.
   *
   * The returned correction factor accounts for a systematic shift in JER if
   * requested.
   */
  double ComputeJerFactor(TLorentzVector const &corrP4, double rho,
                          int rngChannel) const;

  /**
   * \brief Finds matching generator-level jet
   *
   * Returns a nullptr if no match is found within the allowed cone.
   */
  GenJet const *FindGenMatch(TLorentzVector const &p4,
                             double ptResolution) const;

  /**
   * \brief Non-owning pointer to an object that produces generator-level jets
   *
   * May be nullptr.
   */
  GenJetBuilder const *genJetBuilder_;

  /// Minimal pt for jets, GeV
  double minPt_;

  /// Maximal |eta| for jets
  double maxAbsEta_;

  /// Collection of jets
  mutable std::vector<Jet> jets_;

  /// Indicates whether running on simulation or data
  bool isSim_;

  /// Type of requested systematic variation
  Syst syst_;

  /// Direction of requested systematic variation
  SystDirection systDirection_;

  /// Object that computes JEC
  JetCorrector jetCorrector_;

  /**
   * \brief An object to provide JEC uncertainty
   *
   * Only created when a systematic variation in JEC has been requested.
   */
  std::unique_ptr<JetCorrectionUncertainty> jecUncProvider_;

  /**
   * \brief An object that provides pt resolution in simulation
   *
   * Only created when running on simulation.
   */
  std::unique_ptr<JME::JetResolution> jerProvider_;

  /**
   * \brief An object that provides data-to-simulation scale factors for jet pt
   * resolution
   *
   * Only created when running on simulation.
   */
  std::unique_ptr<JME::JetResolutionScaleFactor> jerSFProvider_;

  /// Random number generator
  TabulatedRandomGenerator tabulatedRng_;

  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcMass_;
  mutable TTreeReaderArray<float> srcArea_, srcRawFactor_;
  mutable TTreeReaderArray<float> srcBTag_;
  mutable TTreeReaderArray<int> srcId_;
  mutable TTreeReaderValue<float> puRho_;
  mutable std::unique_ptr<TTreeReaderArray<int>> srcHadronFlavour_;

  // Properties of soft jets, which are used in the type 1 correction of ptmiss
  mutable TTreeReaderArray<float> softRawPt_, softEta_, softPhi_, softArea_;
};

#endif  // JETBUILDER_H_

