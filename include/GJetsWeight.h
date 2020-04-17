#ifndef HZZ2L2NU_INCLUDE_GJETSWEIGHT_H_
#define HZZ2L2NU_INCLUDE_GJETSWEIGHT_H_

#include <WeightBase.h>

#include <PhotonBuilder.h>

/**
 * \brief Computes LO to NLO k factor for Gamma+Jets events
 *
 * We use the same weights than the ones used by JME-17-001. However when the 
 * weight becomes lower than one (at 587.15 GeV) we keep the weight =1. This 
 * looks like the weights we found when comparing our LO samples to our NLO 
 * samples.
 * Notice that the Gamma+Jets samples are used for plotting only, and not for 
 * physics (it is not one of the "genuine MET" processes).
 */
class GJetsWeight : public WeightBase {
 public:
  GJetsWeight(Dataset &dataset, PhotonBuilder const *photonBuilder);

  virtual double NominalWeight() const override;
 
 private:
  bool enabled_;

  /// Non-owning pointer to object that provides collection of photons
  PhotonBuilder const *photonBuilder_;
};

#endif  // HZZ2L2NU_INCLUDE_GJETSWEIGHT_H_
