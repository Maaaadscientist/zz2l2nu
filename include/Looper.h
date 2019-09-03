#ifndef HZZ2L2NU_INCLUDE_LOOPER_H_
#define HZZ2L2NU_INCLUDE_LOOPER_H_

#include <algorithm>
#include <string>

#include <boost/program_options.hpp>

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

  /**
   * \brief Constructs descriptions for command line options
   *
   * Options from \c AnalysisClass are also included.
   */
  static boost::program_options::options_description OptionsDescription();

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
boost::program_options::options_description
Looper<AnalysisClass>::OptionsDescription() {
  namespace po = boost::program_options;
  
  po::options_description optionsDescription{"Dataset"};
  optionsDescription.add_options()
    ("catalog",
     po::value<std::string>()->default_value(
       "/user/npostiau/event_files/MC_ewk/Bonzais-catalog_test_ZZTo2L2Nu-ZZ2l2vPruner.txt"),
     "Path to catalog file")
    ("max-events", po::value<int64_t>()->default_value(-1),
     "Maximal number of events to read; -1 means all")
    ("skip-files", po::value<int>()->default_value(0),
     "Number of files to skip at the beginning of the catalog")
    ("max-files", po::value<int>()->default_value(1),
     "Maximal number of files to read");

  optionsDescription.add(AnalysisClass::OptionsDescription());
  return optionsDescription;
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

