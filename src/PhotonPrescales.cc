#include <PhotonPrescales.h>

PhotonPrescales::PhotonPrescales(Dataset &dataset, Options const &options)
    : isSim_{dataset.Info().IsSimulation()},
      photonTriggers_{GetTriggers(dataset, options)} {}

std::vector<PhotonTrigger> PhotonPrescales::GetTriggers(Dataset &dataset, Options const &options) {
  std::vector<PhotonTrigger> photonTriggers;
  auto const &parentNode = options.GetConfig()["photon_triggers"];
  for (auto &node : parentNode){
    PhotonTrigger currentTrigger;
    currentTrigger.name = node["name"].as<std::string>();
    currentTrigger.threshold = node["threshold"].as<float>();
    currentTrigger.prescale = node["prescale"].as<float>();
    currentTrigger.decision.reset(new TTreeReaderValue<Bool_t>(dataset.Reader(),node["name"].as<std::string>().c_str()));
    photonTriggers.emplace_back(std::move(currentTrigger));
  }
  std::sort(photonTriggers.begin(), photonTriggers.end());
  return photonTriggers;
}

const std::vector<double> PhotonPrescales::GetThresholdsBinning() {
  std::vector<double> binEdges;
  binEdges.emplace_back(0.);
  for (auto const &trigger : photonTriggers_) {
    binEdges.emplace_back(trigger.threshold);
  }
  if (binEdges.back() < 1500.) {
    binEdges.emplace_back(1500.);
  }
  return binEdges;
}

const double PhotonPrescales::GetWeight(double photonPt) {
  double triggerWeight = 0.;
  double expectedTriggerThreshold = 0.;
  int expectedTriggerNum = -1;

  for (unsigned trigNum = 0 ; trigNum < photonTriggers_.size() ; trigNum++) {
    if (photonTriggers_[trigNum].threshold < photonPt){
      expectedTriggerThreshold = photonTriggers_[trigNum].threshold;
      expectedTriggerNum = trigNum;
    }
    else break;
  }
  if (expectedTriggerNum == -1) {
    LOG_WARN << "No expected photon trigger threshold could be found." << std::endl;
  }
  if (!*(*photonTriggers_[expectedTriggerNum].decision)) {
    return 0;
  }
  if (isSim_) {
    triggerWeight = 1.;
  } else {
    triggerWeight = photonTriggers_[expectedTriggerNum].prescale;
  }

  return triggerWeight;
}
