#ifndef HZZ2L2NU_INCLUDE_TRIGGERFILTER_H_
#define HZZ2L2NU_INCLUDE_TRIGGERFILTER_H_

#include <map>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <TBranch.h>
#include <TChain.h>
#include <TTreeReader.h>

#include <Dataset.h>
#include <EventCache.h>
#include <Options.h>
#include <RunSampler.h>


/**
 * \brief Implements trigger selection
 *
 * Triggers are grouped into one or more channels. Each trigger can be
 * restricted to a specific run range, which is modelled with the help of
 * RunSampler. An event passes a channel of the trigger selection if it is
 * accepted by at least one of the triggers active for the sampled run number.
 *
 * The selection is specified in the section \c trigger_filter of the master
 * configuration. It must define a mapping with one or more channels. Each
 * channel is a sequence of blocks. Each block is a mapping with the following
 * fields:
 * - \c triggers  A sequence of one or more trigger names. The names must not
 *   contain the "HLT_" prefix or version postfix.
 * - \c run_range  A sequence of two run numbers that defines the run range for
 *   this group of triggers. This field is optional. If not given, the triggers
 *   are used for all runs.
 * The same trigger may be included in different blocks (and then will be used
 * in each of the associated run ranges) and different channels.
 *
 * Branch with decisions of a given trigger is not guaranteed to be present in
 * every file of a dataset. Because of this, the branches are read directly,
 * bypassing the TTreeReader. They must not be accessed outside of this class in
 * order not to overwrite the buffers associated with the branches.
 */
class TriggerFilter {
 public:
  TriggerFilter(Dataset &dataset, Options const &options,
                RunSampler const *runSampler);

  /**
   * \brief Evaluates decision for the given channel of trigger selection
   *
   * The result is cached on per-event basis.
   */
  bool operator()(std::string_view channel) const;

 private:
  /// Type representing run number
  using run_t = RunSampler::run_t;

  /// Trigger with associated branch and buffer
  struct Trigger {
    /// Comparison needed to use Triger in an std::set
    struct Compare {
      bool operator()(Trigger const &lhs, Trigger const &rhs) const {
        return lhs.name < rhs.name;
      }
    };

    Trigger(std::string_view name_)
      : name{name_}, branch{nullptr} {}

    /// Name of the trigger, without "HLT_" prefix and version postfix
    std::string name;

    /**
     * \brief Branch that contains trigger decision
     *
     * Set to nullptr if the current tree doesn't contain a branch for this
     * trigger.
     */
    mutable TBranch *branch;

    /// Buffer into which the branch will be read
    mutable Bool_t decision;
  };

  /// Trigger with associated run range
  struct TriggerInPeriod {
    TriggerInPeriod(Trigger const *trigger_, run_t minRun_, run_t maxRun_)
      : trigger{trigger_}, minRun{minRun_}, maxRun{maxRun_} {}

    /**
     * \brief Returns true if the trigger is used in given run and accepts
     * current event
     */
    bool GetDecision(run_t run) const;

    /// Non-owning pointer to a trigger object
    Trigger const *trigger;

    /// Run range. Boundaries are included in the range.
    run_t minRun, maxRun;
  };

  /// Collection of triggers used in a channel
  struct Channel {
    /**
     * \brief Computes combined trigger decision for the current event using
     * the set of triggers registered for given run number
     *
     * The event is accepted if at least one of the selected triggers accepts
     * it. The decision is cached in \ref decision.
     */
    bool Collect(run_t run) const;

    /**
     * \brief Triggers with associated run ranges contributing to this channel
     *
     * The same instance of Trigger may enter in this collection multiple times
     * with different run ranges.
     */
    std::vector<TriggerInPeriod> triggers;

    /// Cached result of Collect
    mutable bool decision;
  };

  /**
   * \brief Updates per-event cache
   *
   * If needed, configures trigger branches. Reads trigger decisions.
   */
  void Build() const;

  /// Reads trigger selection from configuration
  void LoadConfig(YAML::Node const config);

  /// Provides run numbers
  RunSampler const *runSampler_;

  /// Implements per-event caching
  EventCache cache_;

  /// Non-owning pointer to the TTreeReader from the dataset
  TTreeReader const *reader_;

  /// Chain associated with the TTreeReader
  TChain *chain_;

  /// Index of the current tree in the chain
  mutable int treeIndex_;

  /**
   * \brief Used triggers
   *
   * Same trigger may contribute to multiple channels.
   */
  std::set<Trigger, Trigger::Compare> triggers_;

  /**
   * \brief Channels for trigger selection
   *
   * The key of the map is the name of the channel as specified in the
   * configuration file.
   */
  std::map<std::string, Channel, std::less<>> channels_;
};

#endif  // HZZ2L2NU_INCLUDE_TRIGGERFILTER_H_
