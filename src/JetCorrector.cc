#include <JetCorrector.h>

#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include <FileInPath.h>
#include "JERC/FactorizedJetCorrector.h"
#include "JERC/JetCorrectorParameters.h"
#include <Logger.h>


JetCorrector::JetCorrector(Dataset &dataset, Options const &options)
    : run_{dataset.Reader(), "run"},
      rho_{dataset.Reader(), "fixedGridRhoFastjetAll"},
      currentIov_{nullptr}, cachedRun_{0} {
  
  auto const jecConfig = Options::NodeAs<YAML::Node>(
      options.GetConfig(), {"jets", "corrections"});
  auto const jecIovsConfig = Options::NodeAs<YAML::Node>(
      jecConfig, {dataset.Info().IsSimulation() ? "sim" : "data"});
  ReadIovParams(jecIovsConfig);
}


JetCorrector::~JetCorrector() noexcept {
  // The destructor needs to be defined at a point where
  // FactorizedJetCorrector is a complete class so that std::unique_ptr knows
  // how to destroy it
}


double JetCorrector::GetJecFull(TLorentzVector const &rawP4,
                                double area) const {
  jetEnergyCorrector_->setJetEta(rawP4.Eta());
  jetEnergyCorrector_->setJetPt(rawP4.Pt());
  jetEnergyCorrector_->setJetA(area);
  jetEnergyCorrector_->setRho(*rho_);

  return jetEnergyCorrector_->getCorrection();
}


double JetCorrector::GetJecL1(TLorentzVector const &rawP4, double area) const {
  jetEnergyCorrector_->setJetEta(rawP4.Eta());
  jetEnergyCorrector_->setJetPt(rawP4.Pt());
  jetEnergyCorrector_->setJetA(area);
  jetEnergyCorrector_->setRho(*rho_);

  return jetEnergyCorrector_->getSubCorrections()[0];
}


void JetCorrector::LoadJec() const {
  LOG_DEBUG << "Loading JEC from the following files:";
  std::vector<JetCorrectorParameters> jecParameters;
  
  for (auto const &path : currentIov_->jecLevels) {
    LOG_DEBUG << "  " << path;
    jecParameters.emplace_back(path);
  }

  jetEnergyCorrector_.reset(new FactorizedJetCorrector(jecParameters));
}


void JetCorrector::ReadIovParams(YAML::Node const config) {
  if (not config.IsSequence())
    throw std::runtime_error("Configuration for JEC IOVs is not a sequence.");

  for (auto const &iovNode : config) {
    IovParams params;

    auto const runRangeNode = iovNode["run_range"];
    if (not runRangeNode) {
      LOG_TRACE << "Found JEC IOV without a run range. "
          << "Will use a catch-all range.";
      params.runRange[0] = 0;
      params.runRange[1] = std::numeric_limits<run_t>::max();
    } else {
      auto const runRange = runRangeNode.as<std::vector<run_t>>();
      if (runRange.size() != 2 or runRange[0] > runRange[1])
        throw std::runtime_error("Illegal run range specified for a JEC IOV.");
      params.runRange[0] = runRange[0];
      params.runRange[1] = runRange[1];
    }

    auto const jecLevelsNode = Options::NodeAs<YAML::Node>(iovNode, {"levels"});
    if (not jecLevelsNode.IsSequence())
      throw std::runtime_error("JEC levels is not a sequence.");
    
    for (auto const &pathNode : jecLevelsNode)
      params.jecLevels.emplace_back(FileInPath::Resolve(
            pathNode.as<std::string>()));

    if (params.jecLevels.empty())
      throw std::runtime_error("JEC levels specified with an empty sequence.");

    iovs_.emplace_back(params);
  }

  // Sort IOVs by the beginnings of their run ranges
  std::sort(
      iovs_.begin(), iovs_.end(),
      [](auto const &lhs, auto const &rhs)
      {return lhs.runRange[0] < rhs.runRange[0];});

  // Sanity checks
  if (iovs_.empty())
    throw std::runtime_error("Sequence of JEC IOV is empty.");

  for (int i = 0; i < iovs_.size() - 1; ++i) {
    auto const &r1 = iovs_[i].runRange;
    auto const &r2 = iovs_[i + 1].runRange;

    if (r1[1] >= r2[0]) {
      std::ostringstream message;
      message << "Overlapping JEC IOVs found: [" << r1[0] << ", " << r1[1]
          << "] and [" << r2[0] << ", " << r2[1] << "].";
      throw std::runtime_error(message.str());
    }
  }
}


void JetCorrector::UpdateIov() const {
  auto const currentRun = *run_;

  if (currentIov_) {
    if (currentRun == cachedRun_)
      return;

    if (currentIov_->Contains(currentRun)) {
      // Run has changed, but it's still within the same IOV
      cachedRun_ = currentRun;
      return;
    }
  }

  currentIov_ = nullptr;
  for (auto const &iov : iovs_) {
    if (iov.Contains(currentRun)) {
      currentIov_ = &iov;
      cachedRun_ = currentRun;
      LoadJec();
      break;
    }
  }

  if (not currentIov_) {
    std::ostringstream message;
    message << "No JEC IOV found for run " << currentRun << ".";
    throw std::runtime_error(message.str());
  }
}

