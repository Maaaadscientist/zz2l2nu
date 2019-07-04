#include <Logger.h>
#include <MelaWeight.h>


MelaWeight::MelaWeight(Dataset &dataset, Options const &options) {
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
      std::ostringstream message;
      message << "Node \"mela_weight\" does not contain mandatory parameter "
          "\"index\".";
      throw std::runtime_error(message.str());
    }

    weightIndex_ = indexNode.as<unsigned>();
  }

  if (enabled_) {
    LOG_DEBUG << "Will apply MELA weight with index " << weightIndex_ << ".";
    weights_.reset(new TTreeReaderArray<float>(dataset.Reader(), "GMELA"));
  } else {
    LOG_DEBUG << "Will not apply MELA weights.";
    enabled_ = false;
  }
}


double MelaWeight::operator()() const {
  if (enabled_)
    return weights_->At(weightIndex_);
  else
    return 1.;
}
