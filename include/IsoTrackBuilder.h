#ifndef ISOTRACKBUILDER_H_
#define ISOTRACKBUILDER_H_

#include <initializer_list>
#include <vector>

#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <CollectionBuilder.h>
#include <Dataset.h>
#include <Options.h>
#include <PhysicsObjects.h>


/**
 * \brief Lazily builds a collection of reconstructed IsoTracks
 *
 * The IsoTrack veto is very effective to reject taus
 *
 */
class IsoTrackBuilder : public CollectionBuilder<IsoTrack> {
 public:
  IsoTrackBuilder(Dataset &dataset, Options const &options);

  /// Returns collection of IsoTracks
  std::vector<IsoTrack> const &Get() const override;

 private:
  /// Constructs IsoTracks for the current event
  void Build() const override;

  /// Minimal pt for IsoTracks to select, GeV
  double minLepPt_, minHadPt_;

  /// Collection of IsoTracks
  mutable std::vector<IsoTrack> IsoTracks_;

  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_;
  mutable TTreeReaderArray<int> srcPdgId_;
  mutable TTreeReaderArray<bool> srcIsPFcand_;
  mutable TTreeReaderArray<float> srcDZ_, srcIso_;

};

#endif  // ISOTRACKBUILDER_H_

