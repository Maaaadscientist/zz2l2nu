#ifndef PHOTONBUILDER_H_
#define PHOTONBUILDER_H_

#include <initializer_list>
#include <vector>

#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <CollectionBuilder.h>
#include <Dataset.h>
#include <Options.h>
#include <PhysicsObjects.h>


/**
 * \brief Lazily builds a collection of reconstructed photons
 *
 * The photons are subjected to a tight kinematical selection.
 */
class PhotonBuilder : public CollectionBuilder<Photon> {
 public:
  PhotonBuilder(Dataset &dataset, Options const &);

  /// Returns collection of photons
  std::vector<Photon> const &Get() const override;

 private:
  /// Constructs photons for the current event
  void Build() const override;

  /// Minimal pt for photons to select, GeV
  double minPt_;

  /// Collection of photons
  mutable std::vector<Photon> photons_;

  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_;
  mutable TTreeReaderArray<int> srcId_;
  mutable TTreeReaderArray<float> srcSigmaIEtaIEta_;
  //mutable TTreeReaderArray<float> srcSigmaIPhiIPhi_;
  mutable TTreeReaderArray<bool> srcIsEtaScEb_, srcHasPixelSeed_;
};

#endif  // PHOTONBUILDER_H_

