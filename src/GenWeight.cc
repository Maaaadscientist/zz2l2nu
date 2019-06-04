#include <GenWeight.h>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <sstream>
#include <stdexcept>


GenWeight::GenWeight(TTreeReader &reader)
  : srcWeights_{reader, "EvtWeights"} {

  // Set up the mapping between the ME scale variations and indices in the
  // vector of weights. The meaning of weights in LHEEventProduct::weights() can
  // be found in [1]. The indices in the source vector are offset by +1 compared
  // to LHEEventProduct::weights().
  // [1] https://github.com/andrey-popov/PEC-tuples/issues/86#issuecomment-481698177
  meScaleIndices_[{Var::Nominal, Var::Nominal}] = 1;
  meScaleIndices_[{Var::Nominal, Var::Up}] = 2;
  meScaleIndices_[{Var::Nominal, Var::Down}] = 3;
  meScaleIndices_[{Var::Up, Var::Nominal}] = 4;
  meScaleIndices_[{Var::Up, Var::Up}] = 5;
  meScaleIndices_[{Var::Up, Var::Down}] = 6;
  meScaleIndices_[{Var::Down, Var::Nominal}] = 7;
  meScaleIndices_[{Var::Down, Var::Up}] = 8;
  meScaleIndices_[{Var::Down, Var::Down}] = 9;
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
    // Return weight with the largest absolute value
    return *std::max_element(
      weights.begin(), weights.end(),
      [](double a, double b){return (std::abs(a) < std::abs(b));});
  else if (direction == Var::Down)
    // Return weight with the smallest absolute value
    return *std::min_element(
      weights.begin(), weights.end(),
      [](double a, double b){return (std::abs(a) < std::abs(b));});
  else
    return 1.;
}


double GenWeight::operator()() const {
  if (srcWeights_.GetSize() > 0)
    return srcWeights_[0];
  else
    return 1.;
}


double GenWeight::RelWeightAlphaS(Var direction) const {
  // In POWHEG gg->ZZ data set [1], LHE weights with indices 109 and 110
  // correspond to PDF sets 265000 (NNPDF30_nlo_as_0117) and 266000
  // (NNPDF30_nlo_as_0119), which thus define the alpha_s variation in PDF. In
  // the stored vector, these weights are found at indices 110 and 111. However,
  // it is not guaranteed that these weights have the same meaning for all
  // samples.
  // [1] /ZZTo2L2Nu_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM
  
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
}


double GenWeight::RelWeightMEScale(Var renorm, Var factor) const {
  if (srcWeights_.GetSize() < 10) {
    std::ostringstream message;
    message << "Cannot access ME scale variations (weights with indices 1 to 9) "
      "because only " << srcWeights_.GetSize() << " weights are available.";
    throw std::runtime_error(message.str());
  }
  
  return srcWeights_[meScaleIndices_.at({renorm, factor})] / srcWeights_[1];
}


double GenWeight::RelWeightPdf(int replica) const {
  if (srcWeights_.GetSize() < 110) {
    std::ostringstream message;
    message << "Cannot access PDF variations (weights with indices 10 to 109) "
      "because only " << srcWeights_.GetSize() << " weights are available.";
    throw std::runtime_error(message.str());
  }

  return srcWeights_[10 + replica] / srcWeights_[1];
}

