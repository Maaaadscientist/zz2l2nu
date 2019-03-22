#include <MuonBuilder.h>

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>

#include <TRandom.h>

#include <Utils.h>


MuonBuilder::MuonBuilder(TTreeReader &reader, Options const &options,
                         TRandom &randomGenerator)
    : minPtLoose_{10.}, minPtTight_{25.}, cache_{reader},
      isSim_{options.GetAs<bool>("is-mc")},
      randomGenerator_{randomGenerator},
      srcPt_{reader, "MuPt"}, srcEta_{reader, "MuEta"},
      srcPhi_{reader, "MuPhi"}, srcE_{reader, "MuE"},
      srcCharge_{reader, "MuCh"}, srcIsolation_{reader, "MuPfIso"},
      srcId_{reader, "MuId"}, srcIdTight_{reader, "MuIdTight"},
      srcTrackerLayers_{reader, "MuTkLayerCnt"},
      genLeptonId_{reader, "GLepBareId"},
      genLeptonPt_{reader, "GLepBarePt"}, genLeptonEta_{reader, "GLepBareEta"},
      genLeptonPhi_{reader, "GLepBarePhi"} {

  std::string const installPath(std::getenv("HZZ2L2NU_BASE"));
  rochesterCorrection_.reset(new RoccoR(installPath + "/data/rcdata.2016.v3/"));
}


std::vector<Muon> const &MuonBuilder::GetLooseMuons() const {
  if (cache_.IsUpdated())
    Update();

  return looseMuons_;
}


std::vector<Muon> const &MuonBuilder::GetTightMuons() const {
  if (cache_.IsUpdated())
    Update();

  return tightMuons_;
}


void MuonBuilder::ApplyRochesterCorrection(
    Muon *muon, int trackerLayers) const {

  // Apply the correction only in its domain of validity
  if (muon->p4.Pt() > 200. or std::abs(muon->p4.Eta()) > 2.4)
    return;


  double scaleFactor = 1.;

  if (isSim_) {
    auto const genMatch = FindGenMatch(*muon, 0.01);

    if (genMatch)
      scaleFactor = rochesterCorrection_->kScaleFromGenMC(
        muon->charge, muon->p4.Pt(), muon->p4.Eta(), muon->p4.Phi(),
        trackerLayers, genMatch->p4.Pt(), randomGenerator_.Uniform());
    else
      scaleFactor = rochesterCorrection_->kScaleAndSmearMC(
        muon->charge, muon->p4.Pt(), muon->p4.Eta(), muon->p4.Phi(),
        trackerLayers, randomGenerator_.Uniform(), randomGenerator_.Uniform());
  } else
    scaleFactor = rochesterCorrection_->kScaleDT(
      muon->charge, muon->p4.Pt(), muon->p4.Eta(), muon->p4.Phi());
  
  
  muon->uncorrP4 = muon->p4;
  muon->p4.SetPtEtaPhiM(muon->p4.Pt() * scaleFactor, muon->p4.Eta(),
                        muon->p4.Phi(), muon->p4.M());
}


std::optional<GenParticle> MuonBuilder::FindGenMatch(
    Muon const &muon, double maxDR) const {

  unsigned iClosest = -1;
  double minDR2 = std::pow(maxDR, 2);

  for (unsigned i = 0; i < genLeptonId_.GetSize(); ++i) {
    if (std::abs(genLeptonId_[i]) != 13)
      // Only consider muons
      continue;

    double const dR2 = utils::DeltaR2(
      muon.p4.Eta(), muon.p4.Phi(), genLeptonEta_[i], genLeptonPhi_[i]);

    if (dR2 < minDR2) {
      iClosest = i;
      minDR2 = dR2;
    }
  }

  if (iClosest != unsigned(-1)) {
    GenParticle matchedParticle{genLeptonId_[iClosest]};
    matchedParticle.p4.SetPtEtaPhiM(
      genLeptonPt_[iClosest], genLeptonEta_[iClosest], genLeptonPhi_[iClosest],
      0.1057
    );
    return matchedParticle;
  } else
    return {};
} 


void MuonBuilder::Update() const {

  looseMuons_.clear();
  tightMuons_.clear();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    bool const passLooseId = srcId_[i] & (1 << 0);

    if (std::abs(srcEta_[i]) > 2.4 or not passLooseId or
        srcIsolation_[i] > 0.25)
      continue;

    Muon muon;
    muon.p4.SetPtEtaPhiE(srcPt_[i], srcEta_[i], srcPhi_[i], srcE_[i]);
    muon.charge = srcCharge_[i];

    ApplyRochesterCorrection(&muon, srcTrackerLayers_[i]);

    if (muon.p4.Pt() < minPtLoose_)
      continue;

    looseMuons_.emplace_back(muon);


    bool const passTightId = srcIdTight_[i] & (1 << 0);  // W.r.t. vertex #0

    if (muon.p4.Pt() < minPtTight_ or not passTightId or
        srcIsolation_[i] > 0.15)
      continue;

    tightMuons_.emplace_back(muon);
  }


  // Make sure the collections are sorted in pt
  std::sort(looseMuons_.begin(), looseMuons_.end(), PtOrdered);
  std::sort(tightMuons_.begin(), tightMuons_.end(), PtOrdered);
}

