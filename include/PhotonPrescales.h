#ifndef HZZ2L2NU_INCLUDE_PHOTONPRESCALES_H_
#define HZZ2L2NU_INCLUDE_PHOTONPRESCALES_H_

#include <memory>
#include <vector>

#include <TTreeReaderValue.h>

#include <Dataset.h>
#include <Options.h>
#include <PhysicsObjects.h>

/// Parameters of a given photon trigger used in the analysis
struct PhotonTrigger {
  bool operator<(const PhotonTrigger& photonTrg) const {
    return (threshold < photonTrg.threshold);
  }
  std::string name;
  double threshold;
  double prescale;
  std::unique_ptr<TTreeReaderValue<Bool_t>> decision;
  std::shared_ptr<std::map<unsigned,std::map<unsigned,int>>> prescaleMap;
};

/**
 * \brief Finds and applies trigger prescales on Single Photon data
 *
 * The names of the triggers used, prescales and thresholds are located in the
 * config file.
 *
 * Performs also a cleaning based on the reconstructed photon pT: if the event
 * does not fire the trigger with the pT threshold just below the photon pT,
 * the event is rejected. This cleaning is performed both on data and MC, while
 * the prescales are applied only on data.
 */
class PhotonPrescales {
 public:
  /**
   * \brief Constructor
   *
   * \param[in] dataset  Dataset that will be processed.
   * \param[in] options  Configuration options.
   */
  PhotonPrescales(Dataset &dataset, Options const &options);

  /**
   * \brief Gets vector of pT thresholds bin boundaries for the weights histograms
   *
   * By construction, first bin starts at 0 and last bin ends at 1500 GeV
   * (unless the last threshold is above 1500).
   */
  std::vector<double> GetThresholdsBinning() const;

  /**
   * \brief Gets the prescale for a given event
   *
   * If the event is from a MC file, the function will only return 0 or 1.
   *
   * Returns 0 if the event fails the check based on the trigger threshold.
   *
   * \param[in] photonPt        p_T of the (first) photon in the event.
   */
  double GetWeight(double photonPt) const;

  /**
   * \brief Gets the prescale for a given run and lumi, from the prescale map
   *
   * If the event is from a MC file, the function will only return 0 or 1.
   *
   * Returns 0 if the event fails the check based on the trigger threshold.
   *
   * \param[in] photonPt   p_T of the (first) photon in the event.
   * \param[in] run        the run number of the event
   * \param[in] lumi       the luminosity block of the event
   */
  int GetPhotonPrescale(double photonPt) const;

 private:
  /**
   * \brief Gets the triggers from the config file
   *
   * Also orders the triggers by increasing pT thresholds (mandatory for 
   * \ref GetWeight).
   */
  static std::vector<PhotonTrigger> GetTriggers(
      Dataset &dataset, Options const &options);

  /**
   * \brief Internal, find the trigger object by the photon Pt
   */
  const PhotonTrigger* FindTrigger(double photonPt) const;

  /// Collection of all possible photon triggers
  std::vector<PhotonTrigger> photonTriggers_;

  /// Indicates if the event is simulation or data
  bool isSim_;

  /// For determining the prescale, work around for const object
  std::unique_ptr<TTreeReaderValue<UInt_t>> ptr_run_;
  std::unique_ptr<TTreeReaderValue<UInt_t>> ptr_lumiBlock_;
};

#endif  // HZZ2L2NU_INCLUDE_PHOTONPRESCALES_H

