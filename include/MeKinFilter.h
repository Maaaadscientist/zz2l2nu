#ifndef MEKINFILTER_H_
#define MEKINFILTER_H_

#include <memory>

#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <Dataset.h>


/**
 * \brief Applies kinematic selection based on matrix element final state
 *
 * This filter is controlled by per-dataset configuration as follows:
 * \code{.yaml}
 * me_kin_filter:
 *   variable: <variable>
 *   range: [<min>, <max>]
 * \endcode
 * The event is accepted if the value of the specified variable lies in the
 * given range. If the per-dataset configuration does not contain section
 * \c me_kin_filter, all events are accepted. The filtering is always disabled
 * in real data.
 */
class MeKinFilter {
 public:
  MeKinFilter(Dataset &dataset);

  /**
   * \brief Evaluates selection for the current event
   *
   * If the filter is disabled for the current dataset, always returns true.
   */
  bool operator()() const;

 private:
  /**
   * \brief Computes parton-level Ht variable
   *
   * Defined as the sum of pt of quarks and gluons in the final state of the
   * matrix element. This is the definition used by MadGraph for cuts
   * \c htjmin and \c htjmax.
   */
  double ComputeHt() const;

  /// Flag indicating whether the filtering is enabled for current dataset
  bool enabled_;

  /// Allowed range for the filtering variable
  double minValue_, maxValue_;

  mutable std::unique_ptr<TTreeReaderValue<UInt_t>> srcNumPart_;
  mutable std::unique_ptr<TTreeReaderArray<Int_t>> srcPdgId_;
  mutable std::unique_ptr<TTreeReaderArray<Int_t>> srcMotherIndex_;
  mutable std::unique_ptr<TTreeReaderArray<Float_t>> srcPartPt_;
};

#endif  // MEKINFILTER_H_

