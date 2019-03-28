#ifndef JETBUILDER_H_
#define JETBUILDER_H_

#include <initializer_list>
#include <memory>
#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <CollectionBuilder.h>
#include <EventCache.h>
#include <Options.h>
#include <PhysicsObjects.h>


// Classes from the JME POG that implement jet corrections are hidden from user
// as in the Pimpl idiom
class JetCorrectionUncertainty;


/// Builds collection of reconstructed jets
class JetBuilder : public CollectionBuilder {
 public:
  JetBuilder(TTreeReader &reader, Options const &options);
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

  /// Minimal pt for jets, GeV
  double minPt_;

  /// Maximal |eta| for jets
  double maxAbsEta_;

  /// Collection of jets
  mutable std::vector<Jet> jets_;

  /// An object to facilitate caching
  EventCache cache_;

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

  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcE_;
  mutable TTreeReaderArray<float> srcBTagCsvV2_, srcHadronFlavour_;
  mutable TTreeReaderArray<float> srcChf_, srcNhf_, srcCemf_, srcNemf_,
    srcNumConstituents_, srcChargedMult_, srcNeutralMult_;
};


inline TLorentzVector const &JetBuilder::GetMomentum(size_t index) const {
  return jets_.at(index).p4;
}


inline size_t JetBuilder::GetNumMomenta() const {
  return jets_.size();
}

#endif  // JETBUILDER_H_

