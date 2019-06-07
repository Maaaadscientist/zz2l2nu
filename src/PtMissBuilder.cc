#include <PtMissBuilder.h>


PtMissBuilder::PtMissBuilder(Dataset &dataset)
    : cache_{dataset.Reader()},
      srcPt_{dataset.Reader(), "METPtType1XY"},
      srcPhi_{dataset.Reader(), "METPhiType1XY"},
      srcSignificance_{dataset.Reader(), "METsig"} {}


PtMiss const &PtMissBuilder::Get() const {
  if (cache_.IsUpdated())
    Build();

  return ptMiss_;
}


void PtMissBuilder::PullCalibration(
    std::initializer_list<CollectionBuilderBase const *> builders) {
  for (auto const *b : builders)
    calibratingBuilders_.emplace_back(b);
}


void PtMissBuilder::Build() const {
  unsigned const typeIndex = 0;  // Corresponds to PF ptmiss
  ptMiss_.p4.SetPtEtaPhiM(srcPt_[typeIndex], 0., srcPhi_[typeIndex], 0.);
  ptMiss_.significance = srcSignificance_[typeIndex];

  for (auto const *builder : calibratingBuilders_)
    ptMiss_.p4 -= builder->GetSumMomentumShift();
}

