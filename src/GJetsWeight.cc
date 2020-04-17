#include <GJetsWeight.h>

#include <algorithm>

#include <Logger.h>

GJetsWeight::GJetsWeight(Dataset &dataset, PhotonBuilder const *photonBuilder)
    : photonBuilder_{photonBuilder} {
  
  auto const settingsNode = dataset.Info().Parameters()["gjets_lo"];

  if (settingsNode and not settingsNode.IsNull()) {
    if (settingsNode.as<bool>())
      enabled_ = true;
  } else
    enabled_ = false;

  if (enabled_) {
    LOG_DEBUG << "Will apply gamma+jets k factor.";
  }
}

double GJetsWeight::NominalWeight() const {
  double kFactor = 1.;
  if (enabled_) {
    for (auto &photon : photonBuilder_->Get()) {
      kFactor *= std::max(1., 1.716910 - 0.001221 * photon.p4.Pt());
    }
  }
  return kFactor;
}
