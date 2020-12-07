#ifndef EVENTNUMBERFILTER_H_
#define EVENTNUMBERFILTER_H_

#include <memory>
#include <unordered_map>
#include <map>

#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <Dataset.h>


/**
 * \brief Applies selection based on run:lumi:event number
 *
 * This filter is controlled by per-dataset configuration as follows:
 * \code{.yaml}
 * \endcode
 * The event is accepted if the value of the specified variable lies in the
 * given range. If the per-dataset configuration does not contain section
 */
class EventNumberFilter {
 public:
  EventNumberFilter(Dataset &dataset, Options const &options);

  /**
   * \brief Evaluates selection for the current event
   *
   * If the filter is disabled for the current dataset, always returns true.
   */
  bool operator()() const;

  using EventMap = std::map<ULong64_t, bool>;
  using RunMap = std::map<UInt_t, EventMap>;

 private:
  /**
   * \brief Computes parton-level Ht variable
   *
   * Defined as the sum of pt of quarks and gluons in the final state of the
   * matrix element. This is the definition used by MadGraph for cuts
   * \c htjmin and \c htjmax.
   */

  /// Flag indicating whether the filtering is enabled for current dataset
  bool enabled_;
  bool isSim_;
  bool quiet_;

  /**
   * \brief Load event list from file
   *
   */
  static RunMap loadEventList(Dataset &dataset, Options const &options);

  RunMap runMap_;
  mutable RunMap::const_iterator eventMap_;

  mutable TTreeReaderValue<UInt_t> run_;
  mutable TTreeReaderValue<UInt_t> lumiBlock_;
  mutable TTreeReaderValue<ULong64_t> event_;

};

#endif  // EVENTNUMBERFILTER_H_
