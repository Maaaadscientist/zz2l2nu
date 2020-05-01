#include <MetFilters.h>


MetFilters::MetFilters(Options const &options, Dataset &dataset) {
  std::vector<std::string> flagsName;
  if (dataset.Info().IsSimulation())
    flagsName = Options::NodeAs<std::vector<std::string>>(
        options.GetConfig(), {"met_filters", "sim"});
  else
    flagsName = Options::NodeAs<std::vector<std::string>>(
        options.GetConfig(), {"met_filters", "data"});

  for (auto const &flagName : flagsName)
    flags_.emplace_back(dataset.Reader(), flagName.c_str());
}


bool MetFilters::operator()() const {
  for (auto &flag: flags_) {
    if (not *flag)
      return false;
  }

  return true;
}

