#include <ElectronBuilder.h>

#include <cmath>
#include <cstdlib>
#include <algorithm>


ElectronBuilder::ElectronBuilder(TTreeReader &reader, Options const &)
    : minPtLoose_{10.}, minPtTight_{25.}, cache_{reader},
      srcPt_{reader, "ElPt"}, srcEta_{reader, "ElEta"},
      srcPhi_{reader, "ElPhi"}, srcE_{reader, "ElE"},
      srcEtaSc_{reader, "ElEtaSc"}, srcCharge_{reader, "ElCh"},
      srcId_{reader, "ElId"} {}


std::vector<Electron> const &ElectronBuilder::GetLooseElectrons() const {
  if (cache_.IsUpdated())
    Update();

  return looseElectrons_;
}


std::vector<Electron> const &ElectronBuilder::GetTightElectrons() const {
  if (cache_.IsUpdated())
    Update();

  return tightElectrons_;
}


void ElectronBuilder::Update() const {

  looseElectrons_.clear();
  tightElectrons_.clear();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    double const absEtaSc = std::abs(srcEtaSc_[i]);
    bool const passLooseId = srcId_[i] & (1 << 1);

    if (srcPt_[i] < minPtLoose_ or absEtaSc > 2.5 or not passLooseId)
      continue;

    Electron electron;
    electron.p4.SetPtEtaPhiE(srcPt_[i], srcEta_[i], srcPhi_[i], srcE_[i]);
    electron.charge = srcCharge_[i];
    electron.etaSc = srcEtaSc_[i];

    looseElectrons_.emplace_back(electron);

    bool const passTightId = srcId_[i] & (1 << 3);

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

