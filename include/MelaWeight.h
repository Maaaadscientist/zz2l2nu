#ifndef MELAWEIGHT_H_
#define MELAWEIGHT_H_

#include <TTreeReaderArray.h>

#include <Dataset.h>
#include <Options.h>

/**
 * \brief Applies MELA weights
 *
 * Intended to be used with datasets with a heavy Higgs boson, both gg fusion
 * and VBF. The index of the weight to be applied is read either from
 * command-line option --mela-weight or from node "mela_weight/index" in the
 * dataset definition. The command-line option takes precedence. If no weight
 * index is given, returns a weight of 1 for every event.
 */
class MelaWeight {
 public:
  MelaWeight(Dataset &dataset, Options const &options);

  /// Self operator obtaining MELA weight
  double operator()() const;

 private:
  bool enabled_;

  /// This index is from the option main code
  unsigned weightIndex_;

  /**
   * \brief MELA weights as read from input file
   *
   * Definition of different bits:
   *  - SOnly & Width==5  GeV : 0
   *  - SOnly & Width==10 GeV : 1
   *  - SOnly & Width==100GeV : 2
   *  - BOnly & Width==5  GeV : 3
   *  - BOnly & Width==10 GeV : 4
   *  - BOnly & Width==100GeV : 5
   *  - BSI   & Width==5  GeV : 6
   *  - BSI   & Width==10 GeV : 7
   *  - BSI   & Width==100GeV : 8
   */
  std::unique_ptr<TTreeReaderArray<float>> weights_;
};

#endif  // MELAWEIGHT_H_

