#include <GenPhotonBuilder.h>

#include <algorithm>


GenPhotonBuilder::GenPhotonBuilder(Dataset &dataset, Options const &)
    : CollectionBuilder{dataset.Reader()},
      srcPhotonGenPartIndex_{dataset.Reader(), "Photon_genPartIdx"},
      srcFlavour_{dataset.Reader(), "Photon_genPartFlav"},
      srcPt_{dataset.Reader(), "GenPart_pt"},
      srcEta_{dataset.Reader(), "GenPart_eta"},
      srcPhi_{dataset.Reader(), "GenPart_phi"},
      srcMass_{dataset.Reader(), "GenPart_mass"} {}


std::vector<GenPhoton> const &GenPhotonBuilder::Get() const {
  Update();
  return photons_;
}


void GenPhotonBuilder::Build() const {
  photons_.clear();

  for (unsigned i = 0 ; i < srcFlavour_.GetSize() ; i++) {
    GenPhoton photon;
    photon.p4.SetPtEtaPhiM(srcPt_.At(srcPhotonGenPartIndex_[i]),
      srcEta_.At(srcPhotonGenPartIndex_[i]),
      srcPhi_.At(srcPhotonGenPartIndex_[i]),
      srcMass_.At(srcPhotonGenPartIndex_[i]));
    photon.flavour = int(srcFlavour_[i]);

    if (IsDuplicate(photon.p4, 0.1))
      continue;

    photons_.emplace_back(photon);
  }

  // Make sure jets are sorted in pt
  std::sort(photons_.begin(), photons_.end(), PtOrdered);
}
