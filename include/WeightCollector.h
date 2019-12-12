#ifndef HZZ2L2NU_INCLUDE_WEIGHTCOLLECTOR_H_
#define HZZ2L2NU_INCLUDE_WEIGHTCOLLECTOR_H_

#include <initializer_list>
#include <string_view>
#include <tuple>
#include <vector>

#include <WeightBase.h>

/**
 * \brief Collects weights from multiple instances of WeightBase
 *
 * An object of this class keeps a sequence of pointers to weight computers,
 * which inherit from WeightBase. Nominal weights given by each computer are
 * multiplied together. Weights corresponding to systematic variations can be
 * accessed via the global variation index that counts consecutively all
 * variations in all registered weight computers.
 */
class WeightCollector {
 public:
  /**
   * \brief Constructor from a list of weight computers
   *
   * Providing weight computers in the argument of the constructor is equivalent
   * to repeated call to method Add.
   */
  WeightCollector(std::initializer_list<WeightBase const *> computers = {});

  /**
   * \brief Registers a new weight computer with this collector
   *
   * The pointed-to object must exist throughout the lifetime of this collector.
   */
  void Add(WeightBase const *computer);

  /**
   * \brief Computes the combined nominal event weight
   *
   * The value is computed as the product of the nominal weights from all
   * registered computers.
   */
  double NominalWeight() const;

  /**
   * \brief Computes combination of default event weights
   *
   * The value is computed as the product of the default weights from all
   * registered computers.
   */
  double operator()() const;

  /// Total number of systematic variations in all registered computers
  int NumVariations() const;

  /**
   * \brief Returns relative weight for the variation with the given global
   * index
   *
   * This relative weight is computed with respect to the nominal weight. The
   * argument must satisfy <tt>0 <= variation < NumVariations()</tt>.
   */
  double RelWeight(int variation) const;

  /**
   * \brief Returns the label of the variation with the given global index
   *
   * The argument must satisfy <tt>0 <= variation < NumVariations()</tt>.
   */
  std::string_view VariationName(int variation) const;

 private:
  /**
   * \brief Converts a global index of a variation into the index of the
   * computer that provides that variation and the index of the variation in
   * that computer.
   */
  std::tuple<int, int> TranslateIndex(int variation) const;

  /// Non-owning pointers to registered weight computers
  std::vector<WeightBase const *> computers_;
};

#endif  // HZZ2L2NU_INCLUDE_WEIGHTCOLLECTOR_H_

