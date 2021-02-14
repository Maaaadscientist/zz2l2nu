#ifndef HZZ2L2NU_INCLUDE_TRIGGERWEIGHT_H_
#define HZZ2L2NU_INCLUDE_TRIGGERWEIGHT_H_

#include <WeightBase.h>

#include <Dataset.h>
#include <ElectronBuilder.h>
#include <MuonBuilder.h>
#include <Options.h>
#include <PhysicsObjects.h>

#include <TH2.h>

/**
 * \brief Applies trigger efficiency scale factors
 *
 * The scale factors and efficiencies are read from files specified in the
 * master configuration. 
 *
 * Systematic variations are also provided.
 */
class TriggerWeight : public WeightBase {
 public:
  /// Constructor
  TriggerWeight(Dataset &dataset, Options const &options,
               ElectronBuilder const *electronBuilder,
               MuonBuilder const *muonBuilder, int efficiencyType = 0);

  ~TriggerWeight();

  /// Computes the total trigger efficiency weight for the current event
  double NominalWeight() const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[0];
  }


  double operator()() const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[defaultWeightIndex_];
  }

  double RelWeight(int variation) const override {

    if (cache_.IsUpdated())
      Update();
    return weights_[variation + 1] / weights_[defaultWeightIndex_];
  }
	
  int NumVariations() const override {
    return 6;
  }
  
  double GetEfficiency(Lepton const *l1, Lepton const *l2, int leptonCat, int syst) const;
	
	/// Computes transfer factor for emu reweighting
  double TransferFactor(Lepton const *l1, Lepton const *l2, int syst = 0) const;

  std::string_view VariationName(int variation) const override;
 
 private:
  /// Computes all weights for the current event
  void Update() const;
 

  /**
   * \brief Store efficiency tables in the map
   *
   * Read the root file with given path and store all
   * histograms together with their names.
   */
  void ReadHistogram(
    std::string const &pathWithName,
		std::map<std::string, std::unique_ptr<TH2D>> &th2Map);
  
  /**
   * \brief Cached weights
   *
   * They are given in the order nominal, 
   * triggerEE_up, triggerEE_down, 
   * triggerMuMu_up, triggerMuMu_down, 
   * triggerEMu_up, triggerEMu_down, 
   */
  mutable std::array<double, 7> weights_;

  EventCache cache_;
 
  /**
   * \brief Index of the default weight to be returned by operator()
   *
   * Corresponds to array \ref weights_.
   */
  int defaultWeightIndex_;
  
  enum EfficiencyType { scaleFactor = 0,
                        dataEfficiency = 1,
                        mcEfficiency = 2 };

  enum LeptonCat { kEE = 0,
                   kMuMu = 1,
                   kEMu = 2 };

  int efficiencyType_;

	/// Trigger efficiency table for 3 channels
  std::map<std::string, std::unique_ptr<TH2D>> 
  mumuScaleFactors_, eeScaleFactors_, emuScaleFactors_;

  /// Non-owning pointer to object that provides collection of electrons
  ElectronBuilder const *electronBuilder_;

  /// Non-owning pointer to object that provides collection of muons
  MuonBuilder const *muonBuilder_;
};

#endif  // HZZ2L2NU_INCLUDE_TRIGGERWEIGHT_H_

