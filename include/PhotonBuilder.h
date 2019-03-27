#ifndef PHOTONBUILDER_H_
#define PHOTONBUILDER_H_

#include <initializer_list>
#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <CollectionBuilder.h>
#include <EventCache.h>
#include <Options.h>
#include <PhysicsObjects.h>

class ElectronBuilder;


/**
 * \brief Builds a collection of reconstructed photons
 *
 * The collection of photons can be cleaned against collections of physics
 * objects produced by other builders. They need to be given to method
 * \ref EnableCleaning.
 *
 * The photons are subjected to a tight kinematical selection.
 */
class PhotonBuilder : public CollectionBuilder {
 public:
  PhotonBuilder(TTreeReader &reader, Options const &);

  /**
   * \brief Enables cleaning with respect to collections produced by given
   * builders
   *
   * Normally, photons should be cleaned against electrons. Given builders
   * should have an appropriate life time.
   */
  void EnableCleaning(
    std::initializer_list<CollectionBuilder const *> builders);

  /// Returns collection of photons
  std::vector<Photon> const &Get() const;

 private:
  /// Constructs photons for the current event
  void Build() const;

  /// Returns momentum of photon with given index
  TLorentzVector const &GetMomentum(size_t index) const override;

  /// Returns the number of photons
  size_t GetNumMomenta() const override;

  /**
   * \brief Checks if the given photon overlaps with an object in one of the
   * collections for cleaning
   *
   * \param[in] photon  Candidate photon.
   * \return Boolean indicating whether the given candidate photon overlaps with
   *   another object.
   *
   * The matching is done in the (eta, phi) metric.
   */
  bool IsDuplicate(Photon const &photon) const;

  /// Minimal pt for photons to select, GeV
  double minPt_;

  /// Collection of photons
  mutable std::vector<Photon> photons_;

  /// An object to facilitate caching
  EventCache cache_;

  /**
   * \brief Collection of non-owning pointers to objects that produce
   * collections against which photons need to be cleaned.
   */
  std::vector<CollectionBuilder const *> buildersForCleaning_;

  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcEtaSc_;
  mutable TTreeReaderArray<unsigned> srcId_;
  mutable TTreeReaderArray<float> srcSigmaIEtaIEta_, srcSigmaIPhiIPhi_;
  mutable TTreeReaderValue<std::vector<bool>> srcHasPixelSeed_;
};


inline TLorentzVector const &PhotonBuilder::GetMomentum(size_t index) const {
  return photons_.at(index).p4;
}


inline size_t PhotonBuilder::GetNumMomenta() const {
  return photons_.size();
}

#endif  // PHOTONBUILDER_H_

