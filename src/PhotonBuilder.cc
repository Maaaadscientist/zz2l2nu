#include <PhotonBuilder.h>

#include <cmath>
#include <cstdlib>
#include <algorithm>

#include <Utils.h>


PhotonBuilder::PhotonBuilder(Dataset &dataset, Options const &)
    : CollectionBuilder{dataset.Reader()},
      minPt_{55.},
      srcPt_{dataset.Reader(), "Photon_pt"}, srcEta_{dataset.Reader(), "Photon_eta"},
      srcPhi_{dataset.Reader(), "Photon_phi"},
      srcIsEtaScEb_{dataset.Reader(), "Photon_isScEtaEB"}, // No direct access to photon SC eta in NanoAOD.
      srcId_{dataset.Reader(), "Photon_cutBased"} {}


std::vector<Photon> const &PhotonBuilder::Get() const {
  Update();
  return photons_;
}


void PhotonBuilder::Build() const {
  photons_.clear();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    // Tight ID
    bool const passId = (srcId_[i] >= 3);
    
    if (srcPt_[i] < minPt_ or not passId)
      continue;
    
    // Only consider photons in the barrel
    if (srcIsEtaScEb_[i])
      continue;

    Photon photon;
    photon.p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], 0.);

    // Perform angular cleaning
    if (IsDuplicate(photon.p4, 0.1))
      continue;

    photons_.emplace_back(photon);
  }

  // Make sure the collection is ordered in pt
  std::sort(photons_.begin(), photons_.end(), PtOrdered);
}

