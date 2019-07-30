#include <GenJetBuilder.h>

#include <algorithm>


GenJetBuilder::GenJetBuilder(Dataset &dataset, Options const &)
    : CollectionBuilder{dataset.Reader()},
      srcPt_{dataset.Reader(), "GenJet_pt"},
      srcEta_{dataset.Reader(), "GenJet_eta"},
      srcPhi_{dataset.Reader(), "GenJet_phi"},
      srcMass_{dataset.Reader(), "GenJet_mass"} {}


std::vector<GenJet> const &GenJetBuilder::Get() const {
  Update();
  return jets_;
}


void GenJetBuilder::Build() const {
  jets_.clear();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    GenJet jet;
    jet.p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], srcMass_[i]);

    if (IsDuplicate(jet.p4, 0.4))
      continue;

    jets_.emplace_back(jet);
  }

  // Make sure jets are sorted in pt
  std::sort(jets_.begin(), jets_.end(), PtOrdered);
}

