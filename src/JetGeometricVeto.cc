#include <JetGeometricVeto.h>

#include <vector>

#include <HZZException.h>
#include <Logger.h>


JetGeometricVeto::JetGeometricVeto(
    Dataset &dataset, Options const &options, JetBuilder const *jetBuilder,
    TabulatedRngEngine &rngEngine)
    : isSim_{dataset.Info().IsSimulation()},
      jetBuilder_{jetBuilder},
      tabulatedRng_{rngEngine},
      srcRun_{dataset.Reader(), "run"} {
  YAML::Node const config = options.GetConfig()["jet_geometric_veto"];
  if (not config) {
    enabled_ = false;
    LOG_DEBUG << "JetGeometricVeto is disabled.";
  } else {
    enabled_ = true;
    LOG_DEBUG << "JetGeometricVeto will be applied.";

    auto const runRange = Options::NodeAs<std::vector<int>>(
        config, {"run_range"});
    if (runRange.size() != 2) {
      HZZException exception;
      exception << "Configuration parameter \"jet_geometric_veto/run_range\" "
          << "must contain exactly two elements but " << runRange.size()
          << " found.";
      throw exception;
    }
    minRun_ = runRange[0];
    maxRun_ = runRange[1];
    if (minRun_ > maxRun_)
      throw HZZException{
          "Wrong ordering in \"jet_geometric_veto/run_range\"."};

    lumiFraction_ = Options::NodeAs<double>(config, {"lumi_fraction"});
    if (lumiFraction_ < 0. or lumiFraction_ > 1.) {
      HZZException exception;
      exception << "Configuration parameter "
          << "\"jet_geometric_veto/lumi_fraction\" must be within range "
          << "[0, 1], but has a value of " << lumiFraction_ << ".";
      throw exception;
    }

    auto const etaRange = Options::NodeAs<std::vector<double>>(
        config, {"eta_range"});
    if (etaRange.size() != 2) {
      HZZException exception;
      exception << "Configuration parameter \"jet_geometric_veto/eta_range\" "
          << "must contain exactly two elements but " << etaRange.size()
          << " found.";
      throw exception;
    }
    minEta_ = etaRange[0];
    maxEta_ = etaRange[1];
    if (minEta_ >= maxEta_)
      throw HZZException{
          "Wrong ordering in \"jet_geometric_veto/eta_range\"."};

    auto const phiRange = Options::NodeAs<std::vector<double>>(
        config, {"phi_range"});
    if (phiRange.size() != 2) {
      HZZException exception;
      exception << "Configuration parameter \"jet_geometric_veto/phi_range\" "
          << "must contain exactly two elements but " << phiRange.size()
          << " found.";
      throw exception;
    }
    minPhi_ = phiRange[0];
    maxPhi_ = phiRange[1];
    if (minPhi_ >= maxPhi_)
      throw HZZException{
          "Wrong ordering in \"jet_geometric_veto/phi_range\"."};
  }

  LOG_TRACE << "JetGeometricVeto configuration:\n"
      << "  run range: [" << minRun_ << ", "
      << maxRun_ << "], lumi fraction: " << lumiFraction_
      << "\n  (eta, phi) window: (" << minEta_ << ", "
      << minPhi_ << ") to (" << maxEta_ << ", " << maxPhi_ << ")";
}


bool JetGeometricVeto::operator()() const {
  if (not enabled_)
    return true;

  if (isSim_) {
    if (tabulatedRng_.Rndm(0) > lumiFraction_)
      return true;
  } else {
    int const run = *srcRun_;
    if (run < minRun_ or run > maxRun_)
      return true;
  }

  for (auto const &jet : jetBuilder_->Get()) {
    double const eta = jet.p4.Eta();
    double const phi = jet.p4.Phi();
    if (eta > minEta_ and eta < maxEta_ and phi > minPhi_ and phi < maxPhi_) {
      LOG_TRACE << "JetGeometricVeto rejects event.";
      return false;
    }
  }

  return true;
}
