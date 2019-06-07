#ifndef JETBUILDER_H_
#define JETBUILDER_H_

#include <initializer_list>
#include <memory>
#include <vector>

#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TRandom.h>

#include <CollectionBuilder.h>
#include <Dataset.h>
#include <GenJetBuilder.h>
#include <Options.h>
#include <PhysicsObjects.h>


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
 * When running over simulation, JER smearing is applied. To follow the standard
 * smearing algorithm, this builder needs to be made aware of generator-level
 * jets via method \ref SetGenJetBuilder. Systematic variations in JEC and JER
 * are supported. The changes in momenta of jets with pt > 15 GeV are aggregated
 * for \ref GetSumMomentumShift.
 */
class JetBuilder : public CollectionBuilder<Jet> {
 public:
  JetBuilder(Dataset &dataset, Options const &options,
             TRandom &randomGenerator);
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
   * \brief Finds matching generator-level jet
   *
   * Returns a nullptr if no match is found within the allowed cone.
   */
  GenJet const *FindGenMatch(Jet const &jet, double ptResolution) const;

  /// Checks whether jet with given index passes PF ID
  bool PassId(unsigned index) const;

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

  /// Reference to common random number generator
  TRandom &randomGenerator_;

  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcE_;
  mutable TTreeReaderArray<float> srcBTagCsvV2_, srcHadronFlavour_;
  mutable TTreeReaderArray<float> srcChf_, srcNhf_, srcCemf_, srcNemf_,
    srcNumConstituents_, srcChargedMult_, srcNeutralMult_;
  mutable TTreeReaderValue<float> puRho_;
};

#endif  // JETBUILDER_H_

