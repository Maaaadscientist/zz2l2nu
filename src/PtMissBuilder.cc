#include <PtMissBuilder.h>

#include <Logger.h>


PtMissBuilder::PtMissBuilder(Dataset &dataset, Options const &options)
    : syst_{Syst::None},
      cache_{dataset.Reader()},
      srcPt_{dataset.Reader(), "RawMET_pt"},
      srcPhi_{dataset.Reader(), "RawMET_phi"},
      srcSignificance_{dataset.Reader(), "MET_significance"} {

  std::string const systLabel{options.GetAs<std::string>("syst")};
  if (systLabel == "metuncl_up")
    syst_ = Syst::UnclEnergyUp;
  else if (systLabel == "metuncl_down")
    syst_ = Syst::UnclEnergyDown;
  else
    syst_ = Syst::None;

  if (syst_ != Syst::None) {
    srcUnclEnergyUpDeltaX_.emplace(
        dataset.Reader(), "MET_MetUnclustEnUpDeltaX");
    srcUnclEnergyUpDeltaY_.emplace(
        dataset.Reader(), "MET_MetUnclustEnUpDeltaY");
    LOG_DEBUG << "Will apply a variation in unclustered momentum in ptmiss.";
  }
}


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
  ptMiss_.significance = *srcSignificance_;

  for (auto const *builder : calibratingBuilders_)
    ptMiss_.p4 -= builder->GetSumMomentumShift();

  if (syst_ == Syst::UnclEnergyUp) {
    ptMiss_.p4.SetPx(ptMiss_.p4.Px() + **srcUnclEnergyUpDeltaX_);
    ptMiss_.p4.SetPy(ptMiss_.p4.Py() + **srcUnclEnergyUpDeltaY_);
  } else if (syst_ == Syst::UnclEnergyDown) {
    ptMiss_.p4.SetPx(ptMiss_.p4.Px() - **srcUnclEnergyUpDeltaX_);
    ptMiss_.p4.SetPy(ptMiss_.p4.Py() - **srcUnclEnergyUpDeltaY_);
  }
}

