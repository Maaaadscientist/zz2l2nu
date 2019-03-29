#ifndef PTMISSBUILDER_H_
#define PTMISSBUILDER_H_

#include <initializer_list>
#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <CollectionBuilder.h>
#include <EventCache.h>
#include <PhysicsObjects.h>


/**
 * \brief Lazily builds reconstructed missing pt
 *
 * Changes to the (type 1 corrected) missing pt caused by corrections applied to
 * other objects, such as jets, can be included using method
 * \ref PullCalibration.
 */
class PtMissBuilder {
 public:
  PtMissBuilder(TTreeReader &reader);

  /// Returns missing pt in the current event
  PtMiss const &Get() const;

  /**
   * \brief Requests that changes in the total momentum produced by the given
   * builders are propagated into missing pt
   *
   * For example, changes in jet momenta due to variations in their corrections
   * can be included this way.
   */
  void PullCalibration(
    std::initializer_list<CollectionBuilderBase const *> builders);

 private:
  /// Constructs ptmiss in the current event
  void Build() const;

  /// An object to facilitate caching
  EventCache cache_;

  /// Objects registered with \ref PullCalibration
  std::vector<CollectionBuilderBase const *> calibratingBuilders_;

  /// Object representing ptmiss in the current event
  mutable PtMiss ptMiss_;

  TTreeReaderArray<float> srcPt_, srcPhi_, srcSignificance_;
 };

#endif  // PTMISSBUILDER_H_

