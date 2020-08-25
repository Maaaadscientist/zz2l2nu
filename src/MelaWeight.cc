#include <MelaWeight.h>

#include <HZZException.h>
#include <Logger.h>


MelaWeight::MelaWeight(Dataset &dataset, Options const &options)
    : enabled_{false} {
  if (options.Exists("mela-weight")) {
    enabled_ = true;
    weightIndex_ = options.GetAs<unsigned>("mela-weight");
  } else if (
      auto const settingsNode = dataset.Info().Parameters()["mela_weight"];
      settingsNode and not settingsNode.IsNull()
  ) {
    enabled_ = true;
    auto const indexNode = settingsNode["index"];

    if (not indexNode) {
      HZZException exception;
      exception << "Node \"mela_weight\" does not contain mandatory parameter "
          "\"index\".";
      throw exception;
    }

    weightIndex_ = indexNode.as<unsigned>();
  }

  if (enabled_) {
    LOG_DEBUG << "Will apply MELA weight with index " << weightIndex_ << ".";
    weights_.reset(new TTreeReaderArray<float>(dataset.Reader(), "GMELA"));
  } else
    LOG_DEBUG << "Will not apply MELA weights.";
}


double MelaWeight::NominalWeight() const {
  if (enabled_)
    return weights_->At(weightIndex_);
  else
    return 1.;
}
