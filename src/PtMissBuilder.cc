#include <PtMissBuilder.h>


PtMissBuilder::PtMissBuilder(Dataset &dataset)
    : cache_{dataset.Reader()},
      //srcSignificance_{dataset.Reader(), "MET_significance"},
      srcPt_{dataset.Reader(), "MET_pt"},
      srcPhi_{dataset.Reader(), "MET_phi"} {}


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
  ptMiss_.p4.SetPtEtaPhiM(*srcPt_, 0., *srcPhi_, 0.);
  //ptMiss_.significance = *srcSignificance_;

  for (auto const *builder : calibratingBuilders_)
    ptMiss_.p4 -= builder->GetSumMomentumShift();
}

