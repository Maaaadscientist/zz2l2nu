#ifndef PTMISSBUILDER_H_
#define PTMISSBUILDER_H_

#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <EventCache.h>
#include <PhysicsObjects.h>


/// Builds reconstructed missing pt
class PtMissBuilder {
 public:
  PtMissBuilder(TTreeReader &reader);

  /// Returns missing pt in the current event
  PtMiss const &Get() const;

 private:
  /// Constructs ptmiss in the current event
  void Build() const;

  /// Object representing ptmiss in the current event
  mutable PtMiss ptMiss_;
  
  /// An object to facilitate caching
  EventCache cache_;

  TTreeReaderArray<float> srcPt_, srcPhi_, srcSignificance_;
 };

#endif  // PTMISSBUILDER_H_

