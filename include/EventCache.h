#ifndef EVENTCACHE_H_
#define EVENTCACHE_H_

#include <TTreeReader.h>


/**
 * \brief Facilitates per-event caching when reading a tree
 *
 * An object of this class tells whether a TTreeReader has moved to a new entry
 * or stayed at the same one as at the time of the previous check.
 */
class EventCache {
 public:
  EventCache(TTreeReader const &reader);

  /// Checks if current entry has been updated since previous invocation
  bool IsUpdated() const;

 private:
  /// Reader object
  TTreeReader const &reader_;

  /// Index of the previously accessed entry in the tree
  mutable long long latestEntry_;
};


inline EventCache::EventCache(TTreeReader const &reader)
    : reader_{reader}, latestEntry_{-1} {}


inline bool EventCache::IsUpdated() const {
  long long curEntry = reader_.GetCurrentEntry();

  if (curEntry == latestEntry_)
    return false;
  else {
    latestEntry_ = curEntry;
    return true;
  }
}


#endif  // EVENTCACHE_H_

