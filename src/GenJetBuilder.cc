#include <GenJetBuilder.h>

#include <algorithm>


GenJetBuilder::GenJetBuilder(TTreeReader &reader, Options const &)
    : CollectionBuilder{reader},
      srcPt_{reader, "GJetAk04Pt"}, srcEta_{reader, "GJetAk04Eta"},
      srcPhi_{reader, "GJetAk04Phi"}, srcE_{reader, "GJetAk04E"} {}


std::vector<GenJet> const &GenJetBuilder::Get() const {
  Update();
  return jets_;
}


void GenJetBuilder::Build() const {
  jets_.clear();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    GenJet jet;
    jet.p4.SetPtEtaPhiE(srcPt_[i], srcEta_[i], srcPhi_[i], srcE_[i]);

    if (IsDuplicate(jet.p4, 0.4))
      continue;

    jets_.emplace_back(jet);
  }

  // Make sure jets are sorted in pt
  std::sort(jets_.begin(), jets_.end(), PtOrdered);
}

