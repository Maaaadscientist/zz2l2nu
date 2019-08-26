#ifndef PILEUPWEIGHT_H_
#define PILEUPWEIGHT_H_

#include <filesystem>
#include <memory>
#include <string>

#include <TH1.h>
#include <TTreeReaderValue.h>

#include <Dataset.h>
#include <Options.h>


/**
 * \brief Performs reweighting for pileup profile
 *
 * Paths to ROOT files containing the pileup profile used to produce simulation
 * and the target profile in data are read from the analysis configuration file.
 * Systematic uncertainty in pileup is implemented using alternative data
 * profiles, which have been constructed assuming different values for the
 * pileup cross section.
 */
class PileUpWeight {
 public:
  /// Constructor
  PileUpWeight(Dataset &dataset, Options const &options);

  /// Computes the pileup weight for the current event
  double operator()() const;

 private:
  /**
   * \brief Reads a histogram with given name from a ROOT file
   *
   * Checks for and reports errors. The returned histogram is owned by the
   * caller.
   */
  static TH1 *ReadHistogram(std::filesystem::path const &path,
                            std::string const &name);

  /**
   * \brief Histograms representing pileup profiles in data and simulation
   *
   * They are normalized to represent probability density. Under- and overflow
   * bins are empty.
   */
  std::unique_ptr<TH1> dataProfile, simProfile;
  
  /**
   * \brief Interface to read the expected number of pileup interactions
   */
  mutable TTreeReaderValue<float> mu_;
};

#endif  // PILEUPWEIGHT_H_

