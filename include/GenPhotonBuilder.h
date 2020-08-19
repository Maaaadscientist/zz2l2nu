#ifndef HZZ2L2NU_INCLUDE_GENPHOTONBUILDER_H_
#define HZZ2L2NU_INCLUDE_GENPHOTONBUILDER_H_

#include <Dataset.h>

#include <TLorentzVector.h>
#include <TTreeReaderArray.h>


/**
 * \brief Reconstruct the photon from the LHE final state
 *
 * This class must only be used with simulation, for the CR.
 */
class GenPhotonBuilder {
 public:
  GenPhotonBuilder(Dataset &dataset);

  /// Computes four-momentum of the photon
  TLorentzVector P4Gamma() const;

 private:
  mutable TTreeReaderArray<Int_t> srcPdgId_;
  mutable TTreeReaderArray<Float_t> srcPt_, srcEta_, srcPhi_, srcMass_;
};

#endif  // HZZ2L2NU_INCLUDE_GENPHOTONBUILDER_H_

