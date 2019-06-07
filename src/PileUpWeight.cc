#include <PileUpWeight.h>

#include <sstream>
#include <stdexcept>

#include <TFile.h>

#include <Logger.h>


PileUpWeight::PileUpWeight(Dataset &dataset, Options const &options)
    : mu_{dataset.Reader(), "EvtPuCntTruth"} {

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

  auto const config = options.GetConfig()["pileup_weight"];

  dataProfile.reset(ReadHistogram(
    Options::NodeAs<std::string>(config["data_profile"]),
    dataProfileName));
  simProfile.reset(ReadHistogram(
    Options::NodeAs<std::string>(config["default_sim_profile"]),
    "pileup"));


  // Make sure the profiles are normalized to represent probability density
  dataProfile->Scale(1. / dataProfile->Integral(), "width");
  simProfile->Scale(1. / simProfile->Integral(), "width");
}


double PileUpWeight::operator()() const {

  // Because of the truncation in mu_ (see the documentation for this data
  // member), it corresponds to the range [*mu_, *mu_ + 1) in the expected
  // number of pileup interactions. Take the centre of this range to compute the
  // weight.
  double const mu = *mu_ + 0.5;

  double const probSim = simProfile->GetBinContent(simProfile->FindFixBin(mu));

  if (probSim <= 0.) {
    LOG_WARN << "Got pileup probability in simulation of " << probSim <<
      " for true pileup of " << mu << ". Set pileup weight to 1.";
    return 1.;
  }

  double const probData = dataProfile->GetBinContent(
    dataProfile->FindFixBin(mu));
  double const weight = probData / probSim;
  LOG_TRACE << "Pileup weight: " << mu << " -> " << weight;

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

