#include <PileUpWeight.h>

#include <sstream>
#include <stdexcept>

#include <TFile.h>

#include <Logger.h>


PileUpWeight::PileUpWeight(Dataset &dataset, Options const &options)
    : mu_{dataset.Reader(), "Pileup_nTrueInt"} {

  // Read pileup profiles in data and simulation. The one in data is chosen
  // based on the requested systematic variation.
  std::string dataProfileName;
  auto const systLabel = options.GetAs<std::string>("syst");

  if (systLabel == "pileup_up")
    dataProfileName = "up";
  else if (systLabel == "pileup_down")
    dataProfileName = "down";
  else
    dataProfileName = "nominal";

  dataProfile.reset(ReadHistogram(
    Options::NodeAs<std::string>(
      options.GetConfig(), {"pileup_weight", "data_profile"}),
    dataProfileName));
  simProfile.reset(ReadHistogram(
    Options::NodeAs<std::string>(
      options.GetConfig(), {"pileup_weight", "default_sim_profile"}),
    "pileup"));


  // Make sure the profiles are normalized to represent probability density
  dataProfile->Scale(1. / dataProfile->Integral(), "width");
  simProfile->Scale(1. / simProfile->Integral(), "width");
}


double PileUpWeight::operator()() const {

  double const probSim = simProfile->GetBinContent(simProfile->FindFixBin(*mu_));

  if (probSim <= 0.) {
    LOG_WARN << "Got pileup probability in simulation of " << probSim <<
      " for true pileup of " << *mu_ << ". Set pileup weight to 1.";
    return 1.;
  }

  double const probData = dataProfile->GetBinContent(
    dataProfile->FindFixBin(*mu_));
  double const weight = probData / probSim;
  LOG_TRACE << "Pileup weight: " << *mu_ << " -> " << weight;

  return weight;
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

