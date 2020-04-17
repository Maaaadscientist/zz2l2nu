#ifndef HZZ2L2NU_INCLUDE_PHOTONWEIGHT_H_
#define HZZ2L2NU_INCLUDE_PHOTONWEIGHT_H_

#include <WeightBase.h>

#include <filesystem>
#include <memory>
#include <string>

#include <TH2.h>

#include <Dataset.h>
#include <Options.h>
#include <PhotonBuilder.h>
#include <PhysicsObjects.h>

/**
 * \brief Applies photon identification efficiency scale factors
 *
 * The computation is done using the tight photon provided by the builder given 
 * to the constructor. The scale factors applied are read from ROOT files 
 * specified in the master configuration.
 * Systematic uncertainties are not (yet) provided.
 */
class PhotonWeight : public WeightBase {
 public:
  /// Constructor
  PhotonWeight(Dataset &dataset, Options const &options,
               PhotonBuilder const *photonBuilder);
  
  /// Computes the total photon efficiency weight for the current event
  virtual double NominalWeight() const override;
  /**
   * \brief Gives scale factor for a photon
   *
   * This function searches in the scale factor table and returns the scale
   * factor for the photon. If pt is higher than the last bin, it will return
   * the value of last bin.
   */
  double PhotonSF(Photon const &photon) const;
    
 private:
  /**
   * \brief Reads a histogram with given path and name
   *
   * Checks for and reports errors. The returned histogram is owned by the
   * caller.
   */  
  static std::unique_ptr<TH2> ReadHistogram(std::string const &pathsWithNames);
  
  /**
   * \brief Histograms of different type of photon scale factors
   *
   * They are stored in 2D format with eta in X-axis and pt in Y-axis.
   */
  std::vector<std::unique_ptr<TH2>> photonTable_; 

  /// Non-owning pointer to object that provides collection of photons
  PhotonBuilder const *photonBuilder_;
};

#endif  // HZZ2L2NU_INCLUDE_PHOTONWEIGHT_H_
