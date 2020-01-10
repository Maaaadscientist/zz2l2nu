#ifndef HZZ2L2NU_INCLUDE_L1TPREFIRING_H_
#define HZZ2L2NU_INCLUDE_L1TPREFIRING_H_

#include <WeightBase.h>

#include <optional>

#include <TTreeReaderValue.h>

#include <Dataset.h>
#include <Options.h>


/**
 * \brief Implements reweighting to account for self-veto due to ECAL L1T
 * prefiring
 *
 * The prefiring problem is described
 * <a href="https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe">here</a>.
 * The reweighting is enabled only when running on simulation and only if the
 * master configuration contains key \c l1t_prefiring set to true.
 */
class L1TPrefiringWeight : public WeightBase {
 public:
  L1TPrefiringWeight(Dataset &dataset, Options const &options);
  
  virtual double NominalWeight() const override {
    return (enabled_) ? **srcWeightNominal_ : 1.;
  }

  virtual int NumVariations() const override {
    return (enabled_) ? 2 : 0;
  }

  virtual double operator()() const override {
    return (enabled_) ? *defaultWeight_->Get() : 1.;
  }

  virtual double RelWeight(int variation) const override;
  virtual std::string_view VariationName(int variation) const override;

 private:
  bool enabled_;
  mutable std::optional<TTreeReaderValue<Float_t>> srcWeightNominal_,
      srcWeightUp_, srcWeightDown_;

  /**
   * \brief Non-owning pointer to one of the source weights
   *
   * Used to provide the default weight by \ref operator(). Only initialized if
   * the reweighting is enabled.
   */
  TTreeReaderValue<Float_t> *defaultWeight_;
};

#endif  // HZZ2L2NU_INCLUDE_L1TPREFIRING_H_

