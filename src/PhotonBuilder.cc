#include <PhotonBuilder.h>

#include <cmath>
#include <cstdlib>
#include <algorithm>

#include <Utils.h>


PhotonBuilder::PhotonBuilder(Dataset &dataset, Options const &)
    : CollectionBuilder{dataset.Reader()},
      minPt_{55.},
      srcPt_{dataset.Reader(), "PhotPt"}, srcEta_{dataset.Reader(), "PhotEta"},
      srcPhi_{dataset.Reader(), "PhotPhi"},
      srcEtaSc_{dataset.Reader(), "PhotScEta"},
      srcId_{dataset.Reader(), "PhotId"},
      srcSigmaIEtaIEta_{dataset.Reader(), "PhotSigmaIetaIeta"},
      srcSigmaIPhiIPhi_{dataset.Reader(), "PhotSigmaIphiIphi"},
      srcHasPixelSeed_{dataset.Reader(), "PhotHasPixelSeed"} {}


std::vector<Photon> const &PhotonBuilder::Get() const {
  Update();
  return photons_;
}


void PhotonBuilder::Build() const {
  photons_.clear();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    // Supposedly some tight ID
    bool const passId = srcId_[i] & (1 << 2);
    
    if (srcPt_[i] < minPt_ or not passId)
      continue;
    
    // Only consider photons in the barrel
    if (std::abs(srcEtaSc_[i]) > 1.4442)
      continue;

    // Additional selection to remove "spikes". Motivation is unclear.
    if (srcSigmaIEtaIEta_[i] < 0.001 or srcSigmaIPhiIPhi_[i] < 0.001)
      continue;

    // Drop photons that have pixel seeds. Motivation is unclear.
    if ((*srcHasPixelSeed_)[i])
      continue;

    Photon photon;
    photon.p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], 0.);
    photon.etaSc = srcEtaSc_[i];

    // Perform angular cleaning
    if (IsDuplicate(photon.p4, 0.1))
      continue;

    photons_.emplace_back(photon);
  }

  // Make sure the collection is ordered in pt
  std::sort(photons_.begin(), photons_.end(), PtOrdered);
}

