#ifndef HZZ2L2NU_INCLUDE_GENZZBUILDER_H_
#define HZZ2L2NU_INCLUDE_GENZZBUILDER_H_

#include <Dataset.h>

#include <TLorentzVector.h>
#include <TTreeReaderArray.h>


/**
 * \brief Reconstruct the ZZ system form the LHE final state
 *
 * Works by summing up four-momenta of all charged and neutral leptons. Will
 * print a warning if an unexpected number of leptons if found.
 *
 * This class must only be used with simulation.
 */
class GenZZBuilder {
 public:
  GenZZBuilder(Dataset &dataset);

  /// Computes four-momentum of the ZZ system
  TLorentzVector P4ZZ() const;

 private:
  mutable TTreeReaderArray<Int_t> srcPdgId_;
  mutable TTreeReaderArray<Float_t> srcPt_, srcEta_, srcPhi_, srcMass_;
};

#endif  // HZZ2L2NU_INCLUDE_GENZZBUILDER_H_

