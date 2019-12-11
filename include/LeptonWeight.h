#ifndef HZZ2L2NU_INCLUDE_LEPTONWEIGHT_H_
#define HZZ2L2NU_INCLUDE_LEPTONWEIGHT_H_

#include <string>

#include <filesystem>
#include <memory>

#include <TH2.h>

#include <Dataset.h>
#include <ElectronBuilder.h>
#include <MuonBuilder.h>
#include <Options.h>
#include <PhysicsObjects.h>

/**
 * \brief Applies lepton identification efficiency scale factors
 *
 * The computation is done using tight leptons provided by the builder given to
 * the constructor. The scale factors applied are read from ROOT files specified
 * in the master configuration. Trigger efficiency scale factors are missing.
 * Systematic uncertainties are not provided.
 */
class LeptonWeight {
 public:
  /// Constructor
  LeptonWeight(Dataset &dataset, Options const &options,
               ElectronBuilder const *electronBuilder,
               MuonBuilder const *muonBuilder);
  
  /// Computes the total lepton efficiency weight for the current event
  double operator()() const;
  
  /**
   * \brief Gives scale factor for an electron
   *
   * This function searches in the scale factor table and returns the scale
   * factor for the electron. If pt is higher than the last bin, it will return
   * the value of last bin.
   */
  double ElectronSF(Electron const &electron) const;
    
  /**
   * \brief Gives scale factor for a muon
   *
   * This function searches in the scale factor table and returns the scale
   * factor for the muon. If pt is higher than the last bin, it will return the
   * value of last bin.
   */
  double MuonSF(Muon const &muon) const;
 
 private:
  /**
   * \brief Reads a histogram with given path and name
   *
   * Checks for and reports errors. The returned histogram is owned by the
   * caller.
   */  
  static std::unique_ptr<TH2> ReadHistogram(std::string const &pathsWithNames);
  
  /**
   * \brief Histograms of different type of lepton scale factors
   *
   * They are stored in 2D format with eta in X-axis and pt in Y-axis.
   */
  std::vector<std::unique_ptr<TH2>> muonTable_, electronTable_; 

  /// Non-owning pointer to object that provides collection of electrons
  ElectronBuilder const *electronBuilder_;

  /// Non-owning pointer to object that provides collection of muons
  MuonBuilder const *muonBuilder_;
};

#endif  // HZZ2L2NU_INCLUDE_LEPTONWEIGHT_H_

