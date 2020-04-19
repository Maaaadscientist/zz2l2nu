#include <PhotonPrescales.h>
#include <FileInPath.h>

PhotonPrescales::PhotonPrescales(Dataset &dataset, Options const &options)
    : photonTriggers_{GetTriggers(dataset, options)},
      isSim_{dataset.Info().IsSimulation()},
      ptr_run_{new TTreeReaderValue<UInt_t>(dataset.Reader(), "run")},
      ptr_lumiBlock_{new TTreeReaderValue<UInt_t>(dataset.Reader(), "luminosityBlock")} {}

std::vector<double> PhotonPrescales::GetThresholdsBinning() const {
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

const PhotonTrigger* PhotonPrescales::FindTrigger(double photonPt) const {
  int expectedTriggerNum = -1;

  for (unsigned trigNum = 0; trigNum < photonTriggers_.size(); trigNum++) {
    if (photonTriggers_[trigNum].threshold < photonPt) {
      expectedTriggerNum = trigNum;
    }
    else break;
  }
  if (expectedTriggerNum == -1) {
    LOG_WARN << "[PhotonPrescales::FindTrigger] No expected photon trigger threshold could be found." << std::endl;
    return nullptr;
  }

  return &photonTriggers_[expectedTriggerNum];
}

double PhotonPrescales::GetWeight(double photonPt) const {
  if (isSim_) return 1;  // fast return if is not data

  const PhotonTrigger* trigger = FindTrigger(photonPt);
  if (!trigger || !*(*trigger->decision)) {
    return 0.;
  }

  return trigger->prescale;
}

int PhotonPrescales::GetPhotonPrescale(double photonPt) const {
  if (isSim_) return 1;  // fast return if is not data

  // First determine which trigger to use, by photonPt, and retrieve the prescale map
  const PhotonTrigger* trigger = FindTrigger(photonPt);
  if (!trigger || !*(*trigger->decision)) {
    return 0.;
  }
  int prescale = 1;

  std::shared_ptr<std::map<unsigned,std::map<unsigned,int>>> runMap = trigger->prescaleMap;
  if (runMap) {
    // For run number, every run shall exisit in the list
    unsigned int run(*(*ptr_run_).Get()), lumi(*(*ptr_lumiBlock_).Get());
    std::map<unsigned,std::map<unsigned,int>>::const_iterator lumiMap = runMap->find(run);
    if (lumiMap == runMap->end()) {
      LOG_WARN << "[PhotonPrescales::GetPhotonPrescale] Cannot find run " << run << " in the prescale table!" << std::endl;
      return int(trigger->prescale);
    }
    // For lumi number, only the lowest lumi of each prescale is recorded
    auto ilumi = (lumiMap->second).upper_bound(lumi);
    if (ilumi != (lumiMap->second).begin()) ilumi--;
    if (lumi < ilumi->first) {
      LOG_WARN << "[PhotonPrescales::GetPhotonPrescale] Cannot find luminosity block " << lumi << " in run " << run
               << " in the prescale table! Reverting to old method!" << std::endl;
      return int(trigger->prescale);
    }
    prescale = ilumi->second;
  }

  return prescale;
}

std::vector<PhotonTrigger> PhotonPrescales::GetTriggers(Dataset &dataset, Options const &options) {
  std::vector<PhotonTrigger> photonTriggers;
  std::string psfilePath = Options::NodeAs<std::string>(options.GetConfig(), {"photon_prescale_map"});
  YAML::Node psfileNode = YAML::LoadFile(FileInPath::Resolve(psfilePath));
  if (!psfileNode) {
    throw std::invalid_argument("[PhotonPrescales::GetTriggers] Cannot find file "+psfilePath);
  }

  auto const &parentNode = options.GetConfig()["photon_triggers"];
  for (auto &node : parentNode){
    PhotonTrigger currentTrigger;
    currentTrigger.name = node["name"].as<std::string>();
    currentTrigger.threshold = node["threshold"].as<float>();
    currentTrigger.prescale = node["prescale"].as<float>();
    currentTrigger.decision.reset(new TTreeReaderValue<Bool_t>(dataset.Reader(),
      node["name"].as<std::string>().c_str()));

    // Loading the prescale map from the yaml file
    YAML::Node trigNode = psfileNode["tab_prescales_"+currentTrigger.name];
    if (!trigNode) {
      LOG_WARN << "[PhotonPrescales::GetTriggers] Cannot find the prescale map for trigger " << currentTrigger.name << std::endl;
      continue;
    }
    // Create the map object
    std::map<unsigned,std::map<unsigned,int>>* psmap = new std::map<unsigned,std::map<unsigned,int>>;
    for (auto rnode : trigNode) {
      unsigned run = rnode.first.as<unsigned>();
      std::map<unsigned,int> lumiMap;
      for (auto lnode : rnode.second) {
        lumiMap[lnode.first.as<unsigned>()] = lnode.second.as<int>();
      }
      (*psmap)[run] = std::move(lumiMap);
    }
    currentTrigger.prescaleMap.reset(psmap);

    photonTriggers.emplace_back(std::move(currentTrigger));
  }
  std::sort(photonTriggers.begin(), photonTriggers.end());
  return photonTriggers;
}
