#ifndef HZZ2L2NU_INCLUDE_LEPTONWEIGHT_H_
#define HZZ2L2NU_INCLUDE_LEPTONWEIGHT_H_

#include <string>

#include <filesystem>
#include <memory>

#include <TH2.h>

#include <Dataset.h>
#include <Options.h>
#include <PhysicsObjects.h>

/**
 * \brief Performs implementation of lepton eff. scale factors
 *
 * Paths to ROOT files containing the 2D histograms used to retrieve scale
 * factors. Systematic uncertainty for now is not included.  Trigger efficiency
 * scale factors are still missing.
 */
class LeptonWeight {
 public:
  /// Constructor
  LeptonWeight(Dataset &dataset, Options const &options);
  
  /// Computes the total lepton efficiency weight for the current event
  double operator()(std::vector<Muon> const &muons,
                    std::vector<Electron> const &electrons) const;
  
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
};

#endif  // HZZ2L2NU_INCLUDE_LEPTONWEIGHT_H_

