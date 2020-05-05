#ifndef METFILTERS_H_
#define METFILTERS_H_

#include <vector>

#include <TTreeReaderValue.h>

#include <Dataset.h>
#include <Options.h>


/**
 * \brief Applies "<a href="https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2">MET filters</a>"
 * to reject pathological events
 */
class MetFilters {
 public:
  MetFilters(Options const &options, Dataset &dataset);

  /// Checks if the current event is good
  bool operator()() const;

 private:
  /// Readers to access decisions of relevant filters
  mutable std::vector<TTreeReaderValue<Bool_t>> flags_;
};

#endif  // METFILTERS_H_

