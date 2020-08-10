#include <GenWeight.h>

#include <algorithm>
#include <initializer_list>
#include <sstream>
#include <stdexcept>


GenWeight::GenWeight(Dataset &dataset, Options const &options)
  : srcLheNominalWeight_{dataset.Reader(), "LHEWeight_originalXWGTUP"},
    srcGenNominalWeight_{dataset.Reader(), "Generator_weight"},
    srcPdfWeights_{dataset.Reader(), "LHEPdfWeight"},
    srcScaleWeights_{dataset.Reader(), "LHEScaleWeight"} {

  DatasetInfo const &info = dataset.Info();
  datasetWeight_ = info.CrossSection()
      / (info.NumEventsTotal() * info.MeanWeight());

  lheScaleWeightsPresent_ = false;
  if (auto const &weightInfo = info.Parameters()["weights"]; weightInfo) {
    if (auto const &lheScale = weightInfo["lhe_scale"]; lheScale)
      lheScaleWeightsPresent_ = lheScale.as<bool>();
  }

  hasBuggedLheWeights_ = false;
  // Drop WG samples of 2017/2018 as they have one fewer weight.
  auto const WGNode = dataset.Info().Parameters()["wgamma_lnugamma"];
  if (WGNode and not WGNode.IsNull())  // Also drop the LO case for consistency
    hasBuggedLheWeights_ = true;

  auto const systLabel = options.GetAs<std::string>("syst");
  if (systLabel == "me_renorm_up")
    defaultVariationIndex_ = (lheScaleWeightsPresent_) ? 0 : -1;
  else if (systLabel == "me_renorm_down")
    defaultVariationIndex_ = (lheScaleWeightsPresent_) ? 1 : -1;
  else if (systLabel == "factor_up")
    defaultVariationIndex_ = (lheScaleWeightsPresent_) ? 2 : -1;
  else if (systLabel == "factor_down")
    defaultVariationIndex_ = (lheScaleWeightsPresent_) ? 3 : -1;
  else
    defaultVariationIndex_ = -1;

  // Set up the mapping between the ME scale variations and indices in the
  // vector of weights. The convention for the weights is extracted from [1].
  // [1] https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc80X_doc.html#LHEScaleWeight
  meScaleIndices_[{Var::Down, Var::Down}] = 0;
  meScaleIndices_[{Var::Down, Var::Nominal}] = 1;
  meScaleIndices_[{Var::Down, Var::Up}] = 2;
  meScaleIndices_[{Var::Nominal, Var::Down}] = 3;
  meScaleIndices_[{Var::Nominal, Var::Nominal}] = 4;
  meScaleIndices_[{Var::Nominal, Var::Up}] = 5;
  meScaleIndices_[{Var::Up, Var::Down}] = 6;
  meScaleIndices_[{Var::Up, Var::Nominal}] = 7;
  meScaleIndices_[{Var::Up, Var::Up}] = 8;
}


double GenWeight::EnvelopeMEScale(Var direction) const {
  // Find weights for all variations in the two ME scales except for the cases
  // when they go in opposite directions
  std::initializer_list<double> const weights{
    RelWeightMEScale(Var::Nominal, Var::Up),
    RelWeightMEScale(Var::Nominal, Var::Down),
    RelWeightMEScale(Var::Up, Var::Nominal),
    RelWeightMEScale(Var::Down, Var::Nominal),
    RelWeightMEScale(Var::Up, Var::Up),
    RelWeightMEScale(Var::Down, Var::Down)
  };

  if (direction == Var::Up)
    return std::max(weights);
  else if (direction == Var::Down)
    return std::min(weights);
  else
    return 1.;
}


double GenWeight::NominalWeight() const {
  return *srcGenNominalWeight_ * datasetWeight_;
}


double GenWeight::RelWeight(int variation) const {
  switch (variation) {
    case 0:
      return RelWeightMEScale(Var::Up, Var::Nominal);
    case 1:
      return RelWeightMEScale(Var::Down, Var::Nominal);
    case 2:
      return RelWeightMEScale(Var::Nominal, Var::Up);
    case 3:
      return RelWeightMEScale(Var::Nominal, Var::Down);
    default:
      throw std::runtime_error("Illegal index.");
  }
}


double GenWeight::RelWeightAlphaS(Var) const {
  // In POWHEG gg->ZZ data set [1], LHE weights with indices 109 and 110
  // correspond to PDF sets 265000 (NNPDF30_nlo_as_0117) and 266000
  // (NNPDF30_nlo_as_0119), which thus define the alpha_s variation in PDF. In
  // the stored vector, these weights are found at indices 110 and 111. However,
  // it is not guaranteed that these weights have the same meaning for all
  // samples.
  // [1] /ZZTo2L2Nu_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM
  
  /*
  if (srcWeights_.GetSize() < 112) {
    std::ostringstream message;
    message << "Cannot access alpha_s variation in PDF (weights with indices "
      "110 and 111) because only " << srcWeights_.GetSize() << " weights are "
      "available.";
    throw std::runtime_error(message.str());
  }

  if (direction == Var::Up)
    return srcWeights_[111] / srcWeights_[1];
  else if (direction == Var::Down)
    return srcWeights_[110] / srcWeights_[1];
  else
    return 1.;
  */
  return 1.; // No alpha_S variation in NanoAOD.
}


double GenWeight::RelWeightMEScale(Var renorm, Var factor) const {
  if (hasBuggedLheWeights_) return 1.;
  if (srcScaleWeights_.GetSize() < 9) {
    std::ostringstream message;
    message << "Cannot access ME scale variations (weights with indices 0 to 8) "
      "because only " << srcScaleWeights_.GetSize() << " weights are available.";
    throw std::runtime_error(message.str());
  }
  
  return srcScaleWeights_[meScaleIndices_.at({renorm, factor})] / srcScaleWeights_[4];
}


double GenWeight::RelWeightPdf(int replica) const {
  if (srcPdfWeights_.GetSize() < 100) {
    std::ostringstream message;
    message << "Cannot access PDF variations (weights with indices 0 to 99) "
      "because only " << srcPdfWeights_.GetSize() << " weights are available.";
    throw std::runtime_error(message.str());
  }

  return srcPdfWeights_[replica]; // Relative wrt nominal weight
}


std::string_view GenWeight::VariationName(int variation) const {
  switch (variation) {
    case 0:
      return "me_renorm_up";
    case 1:
      return "me_renorm_down";
    case 2:
      return "factor_up";
    case 3:
      return "factor_down";
    default:
      return "";
  }
}

