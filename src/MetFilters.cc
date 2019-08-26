#include <MetFilters.h>


MetFilters::MetFilters(Dataset &dataset) {
  // Select flags to check as recommended in [1]
  // [1] https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2?rev=137#2016_data
  for (auto const &flagName: {
      "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter",
      "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter",
      "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter"})
    flags_.emplace_back(dataset.Reader(), flagName);

  if (not dataset.Info().IsSimulation())
    flags_.emplace_back(dataset.Reader(), "Flag_eeBadScFilter");
}


bool MetFilters::operator()() const {
  for (auto &flag: flags_) {
    if (not *flag)
      return false;
  }

  return true;
}

