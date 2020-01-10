#include <PileUpWeight.h>

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include <TFile.h>

#include <Logger.h>


PileUpWeight::PileUpWeight(Dataset &dataset, Options const &options)
    : cache_{dataset.Reader()}, mu_{dataset.Reader(), "Pileup_nTrueInt"} {

  // Read pileup profiles in data and simulation
  std::filesystem::path const path = Options::NodeAs<std::string>(
      options.GetConfig(), {"pileup_weight", "data_profile"});
  dataProfiles_[0].reset(ReadHistogram(path, "nominal"));
  dataProfiles_[1].reset(ReadHistogram(path, "up"));
  dataProfiles_[2].reset(ReadHistogram(path, "down"));
  simProfile_.reset(ReadHistogram(
      Options::NodeAs<std::string>(
          options.GetConfig(), {"pileup_weight", "default_sim_profile"}),
      "pileup"));

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


TH1 *PileUpWeight::ReadHistogram(std::filesystem::path const &path,
                                 std::string const &name) {

  TFile inputFile{path.c_str()};

  if (inputFile.IsZombie()) {
    std::ostringstream message;
    message << "Could not open file " << path << ".";
    throw std::runtime_error(message.str());
  }

  auto hist = dynamic_cast<TH1 *>(inputFile.Get(name.c_str()));

  if (not hist) {
    std::ostringstream message;
    message << "File " << path << " does not contain required histogram \"" <<
      name << "\".";
    throw std::runtime_error(message.str());
  }

  hist->SetDirectory(nullptr);
  inputFile.Close();

  return hist;
}

