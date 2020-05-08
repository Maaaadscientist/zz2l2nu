#include <PileUpWeight.h>

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include <TFile.h>

#include <Logger.h>


PileUpWeight::PileUpWeight(Dataset &dataset, Options const &options)
    : cache_{dataset.Reader()}, mu_{dataset.Reader(), "Pileup_nTrueInt"} {

  YAML::Node const config = options.GetConfig()["pileup_weight"];
  if (not config)
    throw std::runtime_error(
        "Master configuration does not contain section \"pileup_weight\".");

  // Read pileup profiles in data
  std::filesystem::path const path = Options::NodeAs<std::string>(
      config, {"data_profile"});
  dataProfiles_[0].reset(ReadHistogram(path, "nominal"));
  dataProfiles_[1].reset(ReadHistogram(path, "up"));
  dataProfiles_[2].reset(ReadHistogram(path, "down"));

  // Read pileup profiles in simulation. First look for a dataset-specific
  // profile and if it's not found, fall back to the default one.
  if (config["sim_profiles"]) {
    simProfile_.reset(ReadHistogram(
        config["sim_profiles"].as<std::string>(), dataset.Info().Name(),
        false));
    if (simProfile_) {
      LOG_DEBUG << "Will use dataset-specific pileup profile in simulation.";
    } else {
      if (not config["default_sim_profile"]) {
        std::ostringstream message;
        message << "File \"" << config["sim_profiles"].as<std::string>()
            << "\" does not contain pileup profile for dataset \""
            << dataset.Info().Name() << "\" and no default pileup profile is "
            "provided in the master configuration.";
        throw std::runtime_error(message.str());
      }
      simProfile_.reset(ReadHistogram(
          config["default_sim_profile"].as<std::string>(), "pileup"));
      LOG_DEBUG << "Will use default pileup profile in simulation.";
    }
  } else {
    if (not config["default_sim_profile"])
      throw std::runtime_error(
          "Illegal master configuration. Section \"pileup_weight\" must "
          "contain at least one of keys \"sim_profiles\" and "
          "\"default_sim_profile\".");
    simProfile_.reset(ReadHistogram(
        config["default_sim_profile"].as<std::string>(), "pileup"));
  }

  // Make sure the profiles are normalized to represent probability density
  for (auto &profile : dataProfiles_)
    profile->Scale(1. / profile->Integral(), "width");
  simProfile_->Scale(1. / simProfile_->Integral(), "width");

  // The default weight index is chosen based on the requested systematic
  // variation
  auto const systLabel = options.GetAs<std::string>("syst");
  if (systLabel == "pileup_up")
    defaultWeightIndex_ = 1;
  else if (systLabel == "pileup_down")
    defaultWeightIndex_ = 2;
  else
    defaultWeightIndex_ = 0;
  LOG_DEBUG << "Index of default pileup weight: " << defaultWeightIndex_;
}


void PileUpWeight::Update() const {
  double const probSim = simProfile_->GetBinContent(
      simProfile_->FindFixBin(*mu_));

  if (probSim <= 0.) {
    LOG_WARN << "Got pileup probability in simulation of " << probSim <<
      " for true pileup of " << *mu_ << ". Set pileup weights to 1.";
    std::fill(weights_.begin(), weights_.end(), 1.);
  }

  for (int i = 0; i < int(weights_.size()); ++i) {
    double const probData = dataProfiles_[i]->GetBinContent(
        dataProfiles_[i]->FindFixBin(*mu_));
    weights_[i] = probData / probSim;
  }
  LOG_TRACE << "Pileup weights: " << *mu_ << " -> " << weights_[0] << ", "
      << weights_[1] << ", " << weights_[2];
}


std::string_view PileUpWeight::VariationName(int variation) const {
  switch (variation) {
    case 0:
      return "pileup_up";
    case 1:
      return "pileup_down";
    default:
      return "";
  }
}


TH1 *PileUpWeight::ReadHistogram(
    std::filesystem::path const &path, std::string const &name,
    bool checkMissing /*= true*/) {

  TFile inputFile{path.c_str()};

  if (inputFile.IsZombie()) {
    std::ostringstream message;
    message << "Could not open file " << path << ".";
    throw std::runtime_error(message.str());
  }

  auto hist = dynamic_cast<TH1 *>(inputFile.Get(name.c_str()));
  if (hist)
    hist->SetDirectory(nullptr);
  inputFile.Close();

  if (checkMissing and not hist) {
    std::ostringstream message;
    message << "File " << path << " does not contain required histogram \"" <<
      name << "\".";
    throw std::runtime_error(message.str());
  }

  return hist;
}

