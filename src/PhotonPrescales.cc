#include <PhotonPrescales.h>

#include <HZZException.h>
#include <FileInPath.h>


PhotonPrescales::PhotonPrescales(Dataset &dataset, Options const &options)
    : photonTriggers_{GetTriggers(dataset, options)},
      isSim_{dataset.Info().IsSimulation()},
      run_{dataset.Reader(), "run"},
      luminosityBlock_{dataset.Reader(), "luminosityBlock"} {}


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


int PhotonPrescales::GetPhotonPrescale(double photonPt) const {
  // First determine which trigger to use, by photonPt, and retrieve the prescale map
  const PhotonTrigger* trigger = FindTrigger(photonPt);
  //std::cout<<"photon pt="<<photonPt<<std::endl;
  //if (!trigger) std::cout<<"cannot find trigger with this pt"<<std::endl;
  //else std::cout<<"trigger name:"<<trigger->name<<std::endl;
  if (!trigger || !*(*trigger->decision)) {
    //std::cout<<"trigger decision:"<<*(*trigger->decision)<<std::endl;
    return 0.;
  }
  int prescale = 1;

  if (isSim_) return 1;  // fast return if is not data

  // For run number, every run shall exisit in the list
  auto lumiMap = trigger->prescaleMap->find(*run_);
  if (lumiMap == trigger->prescaleMap->end()) {
    LOG_WARN << "[PhotonPrescales::GetPhotonPrescale] Cannot find run " << *run_
             << " in the prescale table!" << std::endl;
    return 1;
  }
  // For lumi number, only the lowest lumi of each prescale is recorded
  auto ilumi = (lumiMap->second).upper_bound(*luminosityBlock_);
  if (ilumi != (lumiMap->second).begin()) ilumi--;
  if (*luminosityBlock_ < ilumi->first) {
    LOG_WARN << "[PhotonPrescales::GetPhotonPrescale] Cannot find luminosity block "
             << *luminosityBlock_ << " in run " << *run_
             << " in the prescale table! Apply default prescale of 1." << std::endl;
    return 1;
  }
  prescale = ilumi->second;

  return prescale;
}


std::vector<PhotonTrigger> PhotonPrescales::GetTriggers(Dataset &dataset, Options const &options) {
  std::vector<PhotonTrigger> photonTriggers;
  std::string psfilePath = Options::NodeAs<std::string>(
      options.GetConfig(), {"photon_triggers", "photon_prescale_map"});
  YAML::Node psfileNode = YAML::LoadFile(FileInPath::Resolve(psfilePath));
  if (!psfileNode) {
    HZZException exception;
    exception << "PhotonPrescales::GetTriggers: Cannot find file \""
        << psfilePath << "\".";
    throw exception;
  }

  auto const &parentNode = options.GetConfig()["photon_triggers"]["triggers"];
  for (auto &node : parentNode){
    PhotonTrigger currentTrigger;
    currentTrigger.name = node["name"].as<std::string>();
    currentTrigger.threshold = node["threshold"].as<float>();
    currentTrigger.decision.reset(new TTreeReaderValue<Bool_t>(dataset.Reader(),
      node["name"].as<std::string>().c_str()));

    // Loading the prescale map from the yaml file
    YAML::Node trigNode = psfileNode[currentTrigger.name];
    if (!trigNode) {
      LOG_WARN << "[PhotonPrescales::GetTriggers] Cannot find the prescale map for trigger "
               << currentTrigger.name << std::endl;
      continue;
    }
    // Create the map object
    auto* psmap = new std::map<unsigned, std::map<unsigned,int>>;
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


const PhotonTrigger* PhotonPrescales::FindTrigger(double photonPt) const {
  int expectedTriggerNum = -1;

  for (unsigned trigNum = 0; trigNum < photonTriggers_.size(); trigNum++) {
    if (photonTriggers_[trigNum].threshold < photonPt) {
      expectedTriggerNum = trigNum;
    }
    else break;
  }
  if (expectedTriggerNum == -1) {
    LOG_WARN << "[PhotonPrescales::FindTrigger] No expected photon trigger threshold could be found."
             << std::endl;
    return nullptr;
  }

  return &photonTriggers_[expectedTriggerNum];
}
