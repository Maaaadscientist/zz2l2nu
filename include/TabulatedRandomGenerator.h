#ifndef HZZ2L2NU_INCLUDE_TABULATEDRANDOMGENERATOR_H_
#define HZZ2L2NU_INCLUDE_TABULATEDRANDOMGENERATOR_H_

#include <limits>
#include <vector>

#include <TTreeReaderValue.h>

#include <Dataset.h>


/**
 * \brief Engine for TabulatedRandomGenerator
 *
 * This class provides fixed random numbers. Multiple numbers per event can be
 * obtained using different channels. The number returned for a given channel in
 * an event with a given ID, is always the same.
 *
 * This class is not expected to be used directly but via
 * TabulatedRandomGenerator instead. Each instance of the latter registers
 * itself as a consumer of random numbers.
 *
 * Internally, a table of pre-generated random integer numbers is maintained.
 * Its elements are indexed by the current event number and given channel
 * number. The table is treated as a circular buffer.
 */
class TabulatedRngEngine {
 public:
  using value_t = uint32_t;

  /// Constructor from a dataset
  TabulatedRngEngine(Dataset &dataset);

  /// Returns maximal possible value of stored random numbers
  constexpr value_t MaxValue() const {
    return std::numeric_limits<value_t>::max();
  }

  /**
   * \brief Reads the random number for the given channel for the current event
   *
   * It is responsibility of the caller to include the offset returned by
   * \ref Register into the provided channel number. No check on the given
   * channel number is done.
   */
  value_t Read(int channel) const;

  /**
   * \brief Registers a new consumer of random numbers
   *
   * The consumer must respect the returned offset to ensure that the set of
   * random numbers it reads does not overlap with the sets read by other
   * consumers.
   *
   * \param[in] numChannels  Number of random numbers that can be used by this
   *     consumer per event.
   * \return Offset that determines the start of the range of channels allocated
   *     for this consumer.
   */
  int Register(int numChannels = 1);

 private:
  /// Total number of channels registered
  int numChannelsRegistered_;

  /// The seed to fill the table of random numbers
  int const kSeed_ = 1439;

  /// Number of random numbers to generate
  int const kVolume_ = 100'000;

  /// Table of pre-generated random numbers
  std::vector<value_t> table_;

  /// Reader to access the event number from event ID
  mutable TTreeReaderValue<ULong64_t> event_;
};


/**
 * \brief Interface to obtain tabulated random numbers
 *
 * Uses TabulatedRngEngine. The number of random numbers consumed per event
 * (number of channels) is a parameter of an object of this class. All methods
 * that return random numbers accept the channel as the first argument; if it is
 * larger or equal than the requested number of channels, it is silently wrapped
 * around.
 */
class TabulatedRandomGenerator {
 public:
  /**
   * \brief Constructor
   *
   * \param[in] engine  Engine to be used. A reference is saved internally.
   * \param[in] numChannels  (Maximal) number of different random numbers
   *     consumed per event.
   */
  TabulatedRandomGenerator(
      TabulatedRngEngine &engine, int numChannels = 1);

  /// Obtain a random number following normal distribution with given parameters
  double Gaus(int channel, double mean = 0., double sigma = 1.) const;

  /// Obtain a random number uniformly distributed on [0, 1]
  double Rndm(int channel) const;

 private:
  /// Underlying engine
  TabulatedRngEngine &engine_;

  /// Number of channels used by this generator
  int numChannels_;

  /// Offset returned by \ref TabulatedRngEngine::Register
  int offset_;
};

#endif  // HZZ2L2NU_INCLUDE_TABULATEDRANDOMGENERATOR_H_

