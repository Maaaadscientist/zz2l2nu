#include <IsoTrackBuilder.h>

#include <cmath>
#include <cstdlib>
#include <algorithm>

#include <Utils.h>


IsoTrackBuilder::IsoTrackBuilder(Dataset &dataset, Options const &)
    : CollectionBuilder{dataset.Reader()},
      minLepPt_{5.}, minHadPt_{10.},
      srcPt_{dataset.Reader(), "IsoTrack_pt"},
      srcEta_{dataset.Reader(), "IsoTrack_eta"},
      srcPhi_{dataset.Reader(), "IsoTrack_phi"},
      srcPdgId_{dataset.Reader(), "IsoTrack_pdgId"},
      srcIsPFcand_{dataset.Reader(), "IsoTrack_isPFcand"},
      srcDZ_{dataset.Reader(), "IsoTrack_dz"},
      srcIso_{dataset.Reader(), "IsoTrack_pfRelIso03_chg"}
{}


std::vector<IsoTrack> const &IsoTrackBuilder::Get() const {
  Update();
  return IsoTracks_;
}


void IsoTrackBuilder::Build() const {
  IsoTracks_.clear();

  // Selection follows: https://indico.cern.ch/event/885275/contributions/3757314/
  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    //std::cout<<"isotrk No."<<i<< " pt="<<srcPt_[i]<<" eta="<<srcEta_[i]<<" phi="<<srcPhi_[i]<<" pdgID="<<srcPdgId_[i]<<std::endl;
    if (!srcIsPFcand_[i]) continue;  // only consider PF candidates
    int id = abs(srcPdgId_[i]);
    if (id != 11 && id != 13 && id < 100) continue;

    float pt = srcPt_[i];
    if (pt < ((id < 14)? minLepPt_ : minHadPt_)) continue; // lower pt treshold for leptons

    float etaThr = ((id == 11)? 3.0 : 2.4); // larger absEta for electron
    if (fabs(srcEta_[i]) > etaThr) continue; 

    if (fabs(srcDZ_[i]) > 0.1) continue;

    if (pt > 50) {
      if (srcIso_[i] * pt > 5) continue;  // consider absIso
    } else {
      if (srcIso_[i] > 0.1) continue;
    }
    IsoTrack isotrk;
    isotrk.p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], 0.);

    // Perform angular cleaning
    if (IsDuplicate(isotrk.p4, 0.3))
      continue;

    IsoTracks_.emplace_back(isotrk);
  }

  // Make sure the collection is ordered in pt
  std::sort(IsoTracks_.begin(), IsoTracks_.end(), PtOrdered);
}

