#include <JetCorrector.h>

#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include <FileInPath.h>
#include "JERC/FactorizedJetCorrector.h"
#include "JERC/JetCorrectorParameters.h"
#include "JERC/JetCorrectionUncertainty.h"
#include "JERC/JetResolution.h"
#include <Logger.h>


JetCorrector::JetCorrector(Dataset &dataset, Options const &options,
                           TabulatedRngEngine &rngEngine)
    : syst_{Syst::None},
      currentIov_{nullptr}, cachedRun_{0},
      tabulatedRng_{rngEngine, 50},
      run_{dataset.Reader(), "run"},
      rho_{dataset.Reader(), "fixedGridRhoFastjetAll"} {
  
  bool const isSim = dataset.Info().IsSimulation();

  auto const jecConfig = Options::NodeAs<YAML::Node>(
      options.GetConfig(), {"jets", "corrections"});
  auto const jecIovsConfig = Options::NodeAs<YAML::Node>(
      jecConfig, {isSim ? "sim" : "data"});
  ReadIovParams(jecIovsConfig);

  if (isSim) {
    jerProvider_.reset(new JME::JetResolution(FileInPath::Resolve(
        Options::NodeAs<std::string>(
            options.GetConfig(),
            {"jets", "resolution", "sim_resolution"}))));
    jerSFProvider_.reset(new JME::JetResolutionScaleFactor(FileInPath::Resolve(
        Options::NodeAs<std::string>(
            options.GetConfig(),
            {"jets", "resolution", "scale_factors"}))));

    std::string const systLabel{options.GetAs<std::string>("syst")};
    if (systLabel == "jec_up") {
      syst_ = Syst::JEC;
      systDirection_ = SystDirection::Up;
    } else if (systLabel == "jec_down") {
      syst_ = Syst::JEC;
      systDirection_ = SystDirection::Down;
    } else if (systLabel == "jer_up") {
      syst_ = Syst::JER;
      systDirection_ = SystDirection::Up;
    } else if (systLabel == "jer_down") {
      syst_ = Syst::JER;
      systDirection_ = SystDirection::Down;
    }

    if (syst_ == Syst::JEC)
      jecUncProvider_.reset(new JetCorrectionUncertainty(FileInPath::Resolve(
          Options::NodeAs<std::string>(
              options.GetConfig(),
              {"jets", "corrections", "uncertainty"}))));
  }
}


JetCorrector::~JetCorrector() noexcept {
  // The destructor needs to be defined at a point where JME classes to access
  // JER and JER are completely defined so that std::unique_ptr knows how to
  // destroy them
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


double JetCorrector::GetJecUncFactor(TLorentzVector const &corrP4) const {
  if (syst_ != Syst::JEC)
    return 1.;

  jecUncProvider_->setJetEta(corrP4.Eta());
  jecUncProvider_->setJetPt(corrP4.Pt());
  double const uncertainty = jecUncProvider_->getUncertainty(true);

  if (systDirection_ == SystDirection::Up)
    return 1. + uncertainty;
  else
    return 1. - uncertainty;
}


double JetCorrector::GetJerFactor(
    TLorentzVector const &corrP4, GenJet const *genJet,
    double ptResolution, int rngChannel) const {
  // Find data-to-simulation scale factor
  Variation jerDirection;

  if (syst_ == Syst::JER) {
    if (systDirection_ == SystDirection::Up)
      jerDirection = Variation::UP;
    else
      jerDirection = Variation::DOWN;
  } else
    jerDirection = Variation::NOMINAL;

  double const jerSF = jerSFProvider_->getScaleFactor(
      {{JME::Binning::JetEta, corrP4.Eta()}}, jerDirection);

  // Depending on the presence of a matching generator-level jet, perform
  // deterministic or stochastic smearing [1]
  // [1] https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=71#Smearing_procedures
  if (genJet) {
    double const jerFactor = 1.
        + (jerSF - 1.) * (corrP4.Pt() - genJet->p4.Pt()) / corrP4.Pt();
    return jerFactor;
  } else {
    double const jerFactor = 1.
        + tabulatedRng_.Gaus(rngChannel, 0., ptResolution)
        * std::sqrt(std::max(std::pow(jerSF, 2) - 1., 0.));
    return jerFactor;
  }
}


double JetCorrector::GetPtResolution(TLorentzVector const &corrP4) const {
  // Relative jet pt resolution in simulation
  double const ptResolution = jerProvider_->getResolution(
      {{JME::Binning::JetPt, corrP4.Pt()},
       {JME::Binning::JetEta, corrP4.Eta()},
       {JME::Binning::Rho, *rho_}});

  return ptResolution;
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

