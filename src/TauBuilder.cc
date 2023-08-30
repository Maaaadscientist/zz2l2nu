#include <TauBuilder.h>

#include <cmath>
#include <cstdlib>
#include <algorithm>

#include <Utils.h>


TauBuilder::TauBuilder(Dataset &dataset, Options const &)
    : CollectionBuilder{dataset.Reader()},
      minLepPt_{18.},
      srcPt_{dataset.Reader(), "Tau_pt"},
      srcEta_{dataset.Reader(), "Tau_eta"},
      srcPhi_{dataset.Reader(), "Tau_phi"},
      srcDecayMode_{dataset.Reader(), "Tau_idDecayModeNewDMs"},
{}


std::vector<Tau> const &TauBuilder::Get() const {
  Update();
  return Tau_;
}


void TauBuilder::Build() const {
  Tau_.clear();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    if (!srcDecayMode_[i]) continue;  // decay modes 5 and 6 are rejected
    if (pt < minLepPt_ ) continue; // lower pt treshold for taus

    if (fabs(srcEta_[i]) > 2.3) continue; 

    Tau tau;
    tau.p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], 0.);

    // Perform angular cleaning
    if (IsDuplicate(tau.p4, 0.4))
      continue;

    Taus_.emplace_back(tau);
  }

  // Make sure the collection is ordered in pt
  std::sort(Taus_.begin(), Taus_.end(), PtOrdered);
}
