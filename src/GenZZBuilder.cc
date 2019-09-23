#include <GenZZBuilder.h>

#include <cmath>
#include <cstdlib>

#include <Logger.h>


GenZZBuilder::GenZZBuilder(Dataset &dataset)
    : srcPdgId_{dataset.Reader(), "LHEPart_pdgId"},
      srcPt_{dataset.Reader(), "LHEPart_pt"},
      srcEta_{dataset.Reader(), "LHEPart_eta"},
      srcPhi_{dataset.Reader(), "LHEPart_phi"},
      srcMass_{dataset.Reader(), "LHEPart_mass"} {}


TLorentzVector GenZZBuilder::P4ZZ() const {
  TLorentzVector sumP4;
  int numLeptons = 0;

  for (int i = 0; i < int(srcPdgId_.GetSize()); ++i) {
    int const absPdgId = std::abs(srcPdgId_[i]);

    // Consider only leptons
    if (absPdgId < 11 or absPdgId > 18)
      continue;

    TLorentzVector p4;
    p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], srcMass_[i]);
    sumP4 += p4;
    ++numLeptons;
  }

  if (numLeptons != 4)
    LOG_WARN << "Found " << numLeptons << " in GenZZBuilder, while expect 4.";

  return sumP4;
}

