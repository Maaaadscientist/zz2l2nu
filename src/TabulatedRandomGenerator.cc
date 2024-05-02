#include <TabulatedRandomGenerator.h>

#include <random>
#include <stdexcept>

#include <boost/math/special_functions/erf.hpp>


TabulatedRngEngine::TabulatedRngEngine(Dataset &dataset)
    : numChannelsRegistered_{0},
      event_{dataset.Reader(), "event"} {

  //std::mt19937 gen{kSeed_};
  std::mt19937 gen{static_cast<std::mt19937::result_type>(kSeed_)};
  std::uniform_int_distribution<value_t> distr{0, MaxValue()};

  table_.reserve(kVolume_);
  for (int i = 0; i < kVolume_; ++i)
    table_.emplace_back(distr(gen));
}


TabulatedRngEngine::value_t TabulatedRngEngine::Read(int channel) const {
  return table_[(*event_ + channel) % kVolume_];
}


int TabulatedRngEngine::Register(int numChannels) {
  numChannelsRegistered_ += numChannels;
  return numChannelsRegistered_ - numChannels;
}


TabulatedRandomGenerator::TabulatedRandomGenerator(
    TabulatedRngEngine &engine, int numChannels)
    : engine_{engine}, numChannels_{numChannels} {
  offset_ = engine.Register(numChannels_);    
}


double TabulatedRandomGenerator::Gaus(
    int channel, double mean, double sigma) const {
  double const p = Rndm(channel);

  // Compute quantile of the standard normal distribution
  double const x = std::sqrt(2) * boost::math::erf_inv(2 * p - 1);

  return mean + sigma * x;
}


double TabulatedRandomGenerator::Rndm(int channel) const {
  return double(engine_.Read(offset_ + channel % numChannels_))
      / engine_.MaxValue();
}

