#include <GenPhotonBuilder.h>

#include <cmath>
#include <cstdlib>

#include <Logger.h>


GenPhotonBuilder::GenPhotonBuilder(Dataset &dataset)
    : srcPdgId_{dataset.Reader(), "LHEPart_pdgId"},
      srcPt_{dataset.Reader(), "LHEPart_pt"},
      srcEta_{dataset.Reader(), "LHEPart_eta"},
      srcPhi_{dataset.Reader(), "LHEPart_phi"},
      srcMass_{dataset.Reader(), "LHEPart_mass"} {}

TLorentzVector GenPhotonBuilder::P4Gamma() const {
  TLorentzVector p4;
  bool photonFound = false;
  for (int i = 0; i < int(srcPdgId_.GetSize()); i++) {
    if (std::abs(srcPdgId_[i]) != 22)
      continue;
    p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], srcMass_[i]);
    if (photonFound) 
      LOG_WARN << "More than 1 photon in this event.";
    photonFound = true;
  }

  if (!photonFound)
    LOG_WARN << "No photon found in this event.";

  return p4;

}
