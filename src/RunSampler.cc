#include <RunSampler.h>

#include <algorithm>
#include <limits>
#include <map>
#include <numeric>

#include <yaml-cpp/yaml.h>

#include <FileInPath.h>
#include <HZZException.h>


RunSampler::RunSampler(Dataset &dataset, Options const &options,
                       TabulatedRngEngine &rngEngine)
    : samplingEnabled_{dataset.Info().IsSimulation()}, cache_{dataset.Reader()},
      tabulatedRng_{rngEngine, 1} {
  if (samplingEnabled_) {
    auto const config = options.GetConfig()["run_sampler"];
    if (not config)
      throw HZZException(
          "Section \"run_sampler\", which is required by RunSampler, is "
          "missing in the master configuration.");
    LoadData(config);
  } else {
    srcRun_.emplace(dataset.Reader(), "run");
  }
}



RunSampler::run_t RunSampler::operator()() const {
  if (cache_.IsUpdated())
    Build();
  return currentRun_;
}


void RunSampler::Build() const {
  if (samplingEnabled_) {
    double const r = tabulatedRng_.Rndm(0);
    auto const res = std::lower_bound(
        cumulProb_.begin(), cumulProb_.end(), r,
        [](auto const &el, auto value){return el.first < value;});
    currentRun_ = res->second;
  } else {
    currentRun_ = **srcRun_;
  }
}


void RunSampler::LoadData(YAML::Node const &config) {
  run_t minRun = std::numeric_limits<run_t>::min();
  run_t maxRun = std::numeric_limits<run_t>::max();
  auto const rangeNode = config["range"];
  if (rangeNode) {
    if (not rangeNode.IsSequence() or rangeNode.size() != 2)
      throw HZZException(
          "Node [run_sampler][range] in the master configuration must be "
          "a sequence of two elements.");
    minRun = rangeNode[0].as<run_t>();
    maxRun = rangeNode[1].as<run_t>();
    if (minRun > maxRun)
      throw HZZException(
          "Incorrect ordering in run range in node [run_sampler][range] in "
          "the master configuration.");
  }

  std::map<run_t, double> luminosity;
  auto const pathNode = config["luminosity"];
  if (not pathNode)
    throw HZZException(
        "Mandatory node [run_sampler][luminosity] is missing in the master "
        "configuration.");
  auto const path = FileInPath::Resolve(pathNode.as<std::string>());
  auto const dataNode = YAML::LoadFile(path);
  for (auto const &runNode : dataNode) {
    auto const run = runNode.first.as<run_t>();
    if (run < minRun or run > maxRun)
      continue;
    luminosity[run] = runNode.second.as<double>();
  }
  if (luminosity.empty())
    throw HZZException("No runs selected in RunSampler.");

  double const norm = std::accumulate(
      luminosity.begin(), luminosity.end(), 0.,
      [](double sum, auto const &el){return sum + el.second;});
  cumulProb_.reserve(luminosity.size());
  double cumulative = 0.;
  for (auto const &[run, lumi] : luminosity) {
    cumulative += lumi;
    cumulProb_.emplace_back(cumulative / norm, run);
  }
}
