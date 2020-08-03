#include <PtMissBuilder.h>

#include <Logger.h>


PtMissBuilder::PtMissBuilder(Dataset &dataset, Options const &options)
    : syst_{Syst::None}, applyEeNoiseMitigation_{false},
      cache_{dataset.Reader()},
      srcPt_{dataset.Reader(), "RawMET_pt"},
      srcPhi_{dataset.Reader(), "RawMET_phi"} {

  auto const config = Options::NodeAs<YAML::Node>(
      options.GetConfig(), {"ptmiss"});
  if (auto const node = config["fix_ee_2017"]; node and node.as<bool>()) {
    applyEeNoiseMitigation_ = true;
    LOG_DEBUG << "Will apply EE noise mitigation in missing pt.";
  }

  if (applyEeNoiseMitigation_) {
    srcDefaultPt_.emplace(dataset.Reader(), "MET_pt");
    srcDefaultPhi_.emplace(dataset.Reader(), "MET_phi");
    srcFixedPt_.emplace(dataset.Reader(), "METFixEE2017_pt");
    srcFixedPhi_.emplace(dataset.Reader(), "METFixEE2017_phi");
    srcSignificance_.emplace(dataset.Reader(), "METFixEE2017_significance");
  } else {
    srcSignificance_.emplace(dataset.Reader(), "MET_significance");
  }

  std::string const systLabel{options.GetAs<std::string>("syst")};
  if (systLabel == "metuncl_up")
    syst_ = Syst::UnclEnergyUp;
  else if (systLabel == "metuncl_down")
    syst_ = Syst::UnclEnergyDown;
  else
    syst_ = Syst::None;

  if (syst_ != Syst::None) {
    std::string const name{(applyEeNoiseMitigation_) ? "METFixEE2017" : "MET"};
    srcUnclEnergyUpDeltaX_.emplace(
        dataset.Reader(), (name + "_MetUnclustEnUpDeltaX").c_str());
    srcUnclEnergyUpDeltaY_.emplace(
        dataset.Reader(), (name + "_MetUnclustEnUpDeltaY").c_str());
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
  ptMiss_.significance = **srcSignificance_;

  if (applyEeNoiseMitigation_) {
    // Add (METFixEE2017 - MET) to the raw missing pt. This difference PF
    // candidates in the EE, both unclustered and also included in soft jets.
    // It also picks up a part of the type 1 correction from MET, which is
    // removed in JetBuilder. See [1] for details.
    // [1] https://hypernews.cern.ch/HyperNews/CMS/get/met/710.html
    TLorentzVector p4;
    p4.SetPtEtaPhiM(**srcFixedPt_, 0., **srcFixedPhi_, 0.);
    ptMiss_.p4 += p4;
    p4.SetPtEtaPhiM(**srcDefaultPt_, 0., **srcDefaultPhi_, 0.);
    ptMiss_.p4 -= p4;
  }

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
