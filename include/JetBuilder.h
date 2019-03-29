#ifndef JETBUILDER_H_
#define JETBUILDER_H_

#include <initializer_list>
#include <memory>
#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TRandom.h>

#include <CollectionBuilder.h>
#include <EventCache.h>
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
 * \brief Builds collection of reconstructed jets
 *
 * Jets are cleaned against physics objects produced by builders given to method
 * \ref EnableCleaning. When running over simulation, JER smearing is applied.
 * To follow the standard smearing algorithm, this builder needs to be made
 * aware of generator-level jets via method \ref SetGenJetBuilder.
 *
 * Systematic variations in JEC and JER are supported.
 */
class JetBuilder : public CollectionBuilder {
 public:
  JetBuilder(TTreeReader &reader, Options const &options,
             TRandom &randomGenerator);
  ~JetBuilder() noexcept;

  /**
   * \brief Enables cleaning with respect to collections produced by given
   * builders
   *
   * Normally, jets should be cleaned against leptons and photons. Given
   * builders should have an appropriate life time.
   */
  void EnableCleaning(
    std::initializer_list<CollectionBuilder const *> builders);
  
  /// Returns collection of jets
  std::vector<Jet> const &Get() const;

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
  void Build() const;

  /**
   * \brief Finds matching generator-level jet
   *
   * Returns a nullptr if no match is found within the allowed cone.
   */
  GenJet const *FindGenMatch(Jet const &jet, double ptResolution) const;

  /// Returns momentum of jet with given index
  TLorentzVector const &GetMomentum(size_t index) const override;

  /// Returns the number of jets
  size_t GetNumMomenta() const override;

  /**
   * \brief Checks if the given jet overlaps with an object in one of the
   * collections for cleaning
   *
   * \param[in] jet  Candidate jet.
   * \return Boolean indicating whether the given candidate jet overlaps with
   *   another object.
   *
   * The matching is done in the (eta, phi) metric.
   */
  bool IsDuplicate(Jet const &jet) const;

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

  /// An object to facilitate caching
  EventCache cache_;

  /// Indicates whether running on simulation or data
  bool isSim_;

  /**
   * \brief Collection of non-owning pointers to objects that produce
   * collections against which jets need to be cleaned.
   */
  std::vector<CollectionBuilder const *> buildersForCleaning_;

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


inline TLorentzVector const &JetBuilder::GetMomentum(size_t index) const {
  return jets_.at(index).p4;
}


inline size_t JetBuilder::GetNumMomenta() const {
  return jets_.size();
}

#endif  // JETBUILDER_H_

