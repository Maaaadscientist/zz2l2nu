#ifndef HZZ2L2NU_INCLUDE_WEIGHTBASE_H_
#define HZZ2L2NU_INCLUDE_WEIGHTBASE_H_

#include <limits>
#include <string_view>

/**
 * \brief Abstract base class for per-event weights
 *
 * All classes implementing per-event weights should inherit from this one to
 * comply with the standard interface.
 */
class WeightBase {
 public:
  virtual ~WeightBase() noexcept = default;

  /// Returns nominal weight for the current event
  virtual double NominalWeight() const = 0;

  /**
   * \brief Returns number of supported systematic variations
   *
   * Up and down variations are counted separately. The nominal weight is not
   * included.
   */
  virtual int NumVariations() const {
    return 0;
  };

  /**
   * \brief Returns the default weight
   *
   * This might be different from the nominal weight if a specific variation
   * has been requested via command line options.
   */
  virtual double operator()() const {
    return NominalWeight();
  };

  /**
   * \brief Returns relative weight for requested systematic variation
   *
   * The relative weight is computed with respect to the nominal weight. The
   * argument must satisfy <tt>0 <= variation < NumVariations()</tt>.
   */
  virtual double RelWeight(int /*variation*/) const {
    return std::numeric_limits<double>::quiet_NaN();
  };

  /**
   * \brief Returns the label of the systematic variation with the given index
   *
   * If the variation has the meaning of an "up" or "down" one, the label must
   * end with "_up" or "_down" respectively. The argument must satisfy
   * <tt>0 <= variation < NumVariations()</tt>.
   */
  virtual std::string_view VariationName(int /*variation*/) const {
    return "";
  };
};

#endif  // HZZ2L2NU_INCLUDE_WEIGHTBASE_H_

