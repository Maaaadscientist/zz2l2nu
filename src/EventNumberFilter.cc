#include <FileInPath.h>
#include <EventNumberFilter.h>

#include <HZZException.h>
#include <Utils.h>

EventNumberFilter::EventNumberFilter(Dataset &dataset, Options const &options)
    : enabled_{true},
      isSim_{dataset.Info().IsSimulation()},
      quiet_{Options::NodeAs<bool>(
          options.GetConfig(), {"photon_filter", "quiet"})},
      runMap_{LoadEventList(dataset, options)},
      eventMap_{runMap_.end()},
      run_{dataset.Reader(), "run"},
      lumiBlock_{dataset.Reader(), "luminosityBlock"},
      event_{dataset.Reader(), "event"}
{}

EventNumberFilter::RunMap EventNumberFilter::LoadEventList(
    Dataset &dataset, Options const &options) {

  EventNumberFilter::RunMap runMap;

  std::string filePath = Options::NodeAs<std::string>(
      options.GetConfig(), {"photon_filter", "file_location"});
  filePath += "/photonFilterList_" + dataset.Info().Name() + ".yaml";
  YAML::Node fileNode = YAML::LoadFile(FileInPath::Resolve(filePath));
  if (!fileNode) {
    HZZException exception;
    exception << "PhotonPrescales::GetTriggers: Cannot find file \""
              << filePath << "\".";
    throw exception;
  }

  // Loading the filter map from the yaml file
  for (auto rnode : fileNode) {
    UInt_t run = rnode.first.as<UInt_t>();
    EventNumberFilter::EventMap eventMap;
    for (auto enode : rnode.second) {
      eventMap[enode.first.as<ULong64_t>()] = bool(enode.second.as<int>());
    }
    runMap[run] = std::move(eventMap);
  }

  return std::move(runMap);
}

bool EventNumberFilter::operator()() const {
  if (not enabled_) return true;

  // In MC run number is alwasy 1, use lumi section in place of the run number
  UInt_t run = (isSim_)? *lumiBlock_ : *run_;

  // For run number, every run shall exisit in the list
  // Cache for the eventMap
  if (eventMap_ == runMap_.end() || run != eventMap_->first)
    eventMap_ = runMap_.find(run);

  if (eventMap_ == runMap_.end()) {
    if (not quiet_)
      LOG_WARN << "[EventNumberFilter] Cannot find run " << *run_
               << " in the even list!" << std::endl;
    return true;
  }

  auto item = (eventMap_->second).find(*event_);
  if (item == (eventMap_->second).end()) {
    if (not quiet_)
      LOG_WARN << "[EventNumberFilter] Cannot find run: " << *run_
               << " event: " << *event_ << " in the event list!" << std::endl;
    return true;
  }

  return item->second;
}
