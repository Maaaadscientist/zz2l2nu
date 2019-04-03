#include <PtMissBuilder.h>


PtMissBuilder::PtMissBuilder(TTreeReader &reader)
    : cache_{reader},
      srcPt_{reader, "METPtType1XY"}, srcPhi_{reader, "METPhiType1XY"},
      srcSignificance_{reader, "METsig"} {}


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

