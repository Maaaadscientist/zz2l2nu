#include <L1TPrefiringWeight.h>

#include <HZZException.h>
#include <Logger.h>


L1TPrefiringWeight::L1TPrefiringWeight(
    Dataset &dataset, Options const &options) {
  if (dataset.Info().IsSimulation()
      and Options::NodeAs<bool>(options.GetConfig(), {"l1t_prefiring"})) {
    enabled_ = true;
    srcWeightNominal_.emplace(dataset.Reader(), "L1PreFiringWeight_Nom");
    srcWeightUp_.emplace(dataset.Reader(), "L1PreFiringWeight_Up");
    srcWeightDown_.emplace(dataset.Reader(), "L1PreFiringWeight_Dn");
    LOG_DEBUG << "Weights to account for L1T prefiring will be applied.";
  } else {
    enabled_ = false;
  }

  if (enabled_) {
    auto const systLabel = options.GetAs<std::string>("syst");
    if (systLabel == "l1t_prefiring_up")
      defaultWeight_ = &srcWeightUp_.value();
    else if (systLabel == "l1t_prefiring_down")
      defaultWeight_ = &srcWeightDown_.value();
    else
      defaultWeight_ = &srcWeightNominal_.value();
  }
}


double L1TPrefiringWeight::RelWeight(int variation) const {
  switch (variation) {
    case 0:
      return **srcWeightUp_ / **srcWeightNominal_;
    case 1:
      return **srcWeightDown_ / **srcWeightNominal_;
    default:
      throw HZZException{"Illegal index."};
  }
}


std::string_view L1TPrefiringWeight::VariationName(int variation) const {
  switch (variation) {
    case 0:
      return "l1t_prefiring_up";
    case 1:
      return "l1t_prefiring_down";
    default:
      return "";
  }
}
