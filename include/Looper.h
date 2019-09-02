#ifndef HZZ2L2NU_INCLUDE_LOOPER_H_
#define HZZ2L2NU_INCLUDE_LOOPER_H_

#include <algorithm>
#include <string>

#include <Dataset.h>
#include <Logger.h>
#include <Options.h>


/**
 * \brief Loops over events in input dataset and feeds them to an analysis class
 *
 * \tparam AnalysisClass  The class representing an analysis. Must define a
 *   constructor of the form <tt>AnalysisClass(Options const &, Dataset &)</tt>.
 *   Method <tt>bool AnalysisClass::ProcessEvent()</tt> is called for each event
 *   in the input dataset. It must return true if the event passes the selection
 *   and false otherwise. Method <tt>void AnalysisClass::PostProcessing()</tt>
 *   is called after the event loop.
 */
template<typename AnalysisClass>
class Looper {
 public:
  /// Constructor from a configuration
  Looper(Options const &options);

  /// Runs over the input dataset
  void Run();
 
 private:
  /// Input dataset constructed from the configuration
  Dataset dataset_;

  /// Implementation of the analysis
  AnalysisClass analysis_;

  /// The number of events to read from the input dataset
  int64_t numEvents_;
};


template<typename AnalysisClass>
Looper<AnalysisClass>::Looper(Options const &options)
    : dataset_{DatasetInfo{options.GetAs<std::string>("catalog"), options},
               options.GetAs<int>("skip-files"),
               options.GetAs<int>("max-files")},
      analysis_{options, dataset_} {
  
  auto const maxEvents = options.GetAs<int64_t>("max-events");

  if (maxEvents >= 0)
    numEvents_ = std::min(maxEvents, dataset_.NumEntries());
  else
    numEvents_ = dataset_.NumEntries();
}


template<typename AnalysisClass>
void Looper<AnalysisClass>::Run() {
  LOG_DEBUG << "Will run over " << numEvents_ << " events.";
  int64_t numSelected = 0;

  for (int64_t iEvent = 0; iEvent < numEvents_; ++iEvent) {
    if (iEvent % 10000 == 0) {
      LOG_INFO << Logger::TimeStamp << " Event " << iEvent << " out of "
          << numEvents_;
    }
    dataset_.SetEntry(iEvent);
    if (analysis_.ProcessEvent())
      ++numSelected;
  }

  analysis_.PostProcessing();
  LOG_INFO << Logger::TimeStamp << " Finishing. Total events selected: "
      << numSelected << ".";
}

#endif  // HZZ2L2NU_INCLUDE_LOOPER_H_

