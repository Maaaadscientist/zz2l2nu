#include <ElectronBuilder.h>

#include <cmath>
#include <cstdlib>
#include <algorithm>


ElectronBuilder::ElectronBuilder(Dataset &dataset, Options const &)
    : CollectionBuilder{dataset.Reader()},
      minPtLoose_{10.}, minPtTight_{25.},
      srcPt_{dataset.Reader(), "Electron_pt"},
      srcEta_{dataset.Reader(), "Electron_eta"},
      srcPhi_{dataset.Reader(), "Electron_phi"},
      srcMass_{dataset.Reader(), "Electron_mass"},
      srcDeltaEtaSc_{dataset.Reader(), "Electron_deltaEtaSC"},
      srcCharge_{dataset.Reader(), "Electron_charge"},
      srcId_{dataset.Reader(), "Electron_cutBased"},
      srcECorr_{dataset.Reader(), "Electron_eCorr"} {}


std::vector<Electron> const &ElectronBuilder::GetLoose() const {
  Update();
  return looseElectrons_;
}


std::vector<Electron> const &ElectronBuilder::GetTight() const {
  Update();
  return tightElectrons_;
}


void ElectronBuilder::Build() const {

  looseElectrons_.clear();
  tightElectrons_.clear();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    double const etaSc = srcDeltaEtaSc_[i] + srcEta_[i];
    double const absEtaSc = std::abs(etaSc);
    bool const passLooseId = (srcId_[i] >= 2);

    if (srcPt_[i] < minPtLoose_ or absEtaSc > 2.5 or not passLooseId)
      continue;

    Electron electron;
    electron.p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], srcMass_[i]);
    electron.charge = srcCharge_[i];
    electron.etaSc = etaSc;

    if (IsDuplicate(electron.p4, 0.1))
      continue;

    looseElectrons_.emplace_back(electron);

    // Propagate corrections to momenta of loose electrons to ptmiss
    TLorentzVector const uncorrP4 = electron.p4 * (1. / srcECorr_[i]);
    AddMomentumShift(uncorrP4, electron.p4);

    bool const passTightId = srcId_[i] >= 4;

    if (srcPt_[i] < minPtTight_ or not passTightId)
      continue;

    if (absEtaSc > 1.4442 and absEtaSc < 1.5660)  // EB-EE gap
      continue;

    tightElectrons_.emplace_back(electron);
  }

  // Make sure the collections are ordered in pt
  std::sort(looseElectrons_.begin(), looseElectrons_.end(), PtOrdered);
  std::sort(tightElectrons_.begin(), tightElectrons_.end(), PtOrdered);
}

