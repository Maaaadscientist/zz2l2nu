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
  PhotonBuilder(Dataset &dataset);

  /// Returns collection of photons
  std::vector<Photon> const &Get() const override;

 private:
  /// Constructs photons for the current event
  void Build() const override;

  /// Minimal pt for photons to select, GeV
  double minPt_;

  /// Indicates whether running on simulation or data
  bool isSim_;

  /// Collection of photons
  mutable std::vector<Photon> photons_;

  mutable std::unique_ptr<TTreeReaderArray<float>> srcGenPt_, srcGenEta_;
  mutable std::unique_ptr<TTreeReaderArray<float>> srcGenPhi_;
  mutable std::unique_ptr<TTreeReaderArray<int>> srcPhotonGenPartIndex_;
  mutable std::unique_ptr<TTreeReaderArray<UChar_t>> srcFlavour_;
  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_;
  mutable std::unique_ptr<TTreeReaderArray<int>> srcId_;
  mutable TTreeReaderArray<bool> srcIsEtaScEb_;
};

#endif  // PHOTONBUILDER_H_

