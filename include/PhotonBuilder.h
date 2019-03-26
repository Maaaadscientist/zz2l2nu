#ifndef PHOTONBUILDER_H_
#define PHOTONBUILDER_H_

#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <EventCache.h>
#include <Options.h>
#include <PhysicsObjects.h>

class ElectronBuilder;


/**
 * \brief Builds a collection of reconstructed photons
 *
 * The collection of photons can be cleaned against tight electrons produced by
 * an instance of ElectronBuilder. To do this, method \ref EnableCleaning must
 * be called.
 *
 * The photons are subjected to a tight kinematical selection.
 */
class PhotonBuilder {
 public:
  PhotonBuilder(TTreeReader &reader, Options const &);

  /**
   * \brief Enables cleaning with respect to electrons
   *
   * The ElectronBuilder object must have an appropriate life time. It is not
   * owned by this.
   */
  void EnableCleaning(ElectronBuilder const *electronBuilder);

  /// Returns collection of photons
  std::vector<Photon> const &Get() const;

 private:
  /// Constructs photons for the current event
  void Build() const;

  /**
   * \brief Checks if the given photon overlaps with an electron
   *
   * \param[in] photon  Candidate photon.
   * \return Boolean indicating whether the given candidate photon overlaps with
   *   an electron.
   *
   * Tight electrons produced by the registered ElectronBuilder are checked. The
   * matching is done in the (eta, phi) metric. If no ElectronBuilder has been
   * registered, always returns false.
   */
  bool IsDuplicate(Photon const &photon) const;

  /// Minimal pt for photons to select, GeV
  double minPt_;

  /// Collection of photons
  mutable std::vector<Photon> photons_;

  /// An object to facilitate caching
  EventCache cache_;

  /**
   * \brief Non-owning pointer to an object that constructs electrons
   *
   * Used to clean photons against the electrons. May be a nullptr; in that case
   * no cleaning is performed.
   */
  ElectronBuilder const *electronBuilder_;

  TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcEtaSc_;
  TTreeReaderArray<unsigned> srcId_;
  TTreeReaderArray<float> srcSigmaIEtaIEta_, srcSigmaIPhiIPhi_;
  mutable TTreeReaderValue<std::vector<bool>> srcHasPixelSeed_;
};

#endif  // PHOTONBUILDER_H_

