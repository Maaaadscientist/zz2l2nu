#ifndef JETBUILDER_H_
#define JETBUILDER_H_

#include <initializer_list>
#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <CollectionBuilder.h>
#include <EventCache.h>
#include <Options.h>
#include <PhysicsObjects.h>


/// Builds collection of reconstructed jets
class JetBuilder : public CollectionBuilder {
 public:
  JetBuilder(TTreeReader &reader, Options const &);

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

