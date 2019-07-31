#include <MuonBuilder.h>

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>

#include <TRandom.h>

#include <FileInPath.h>
#include <Utils.h>


MuonBuilder::MuonBuilder(Dataset &dataset, Options const &options,
                         TRandom &randomGenerator)
    : CollectionBuilder{dataset.Reader()},
      minPtLoose_{10.}, minPtTight_{25.},
      isSim_{dataset.Info().IsSimulation()},
      randomGenerator_{randomGenerator},
      srcPt_{dataset.Reader(), "Muon_pt"}, srcEta_{dataset.Reader(), "Muon_eta"},
      srcPhi_{dataset.Reader(), "Muon_phi"}, srcMass_{dataset.Reader(), "Muon_mass"},
      srcCharge_{dataset.Reader(), "Muon_charge"},
      srcIsolation_{dataset.Reader(), "Muon_pfRelIso04_all"},
      srcIsPfMuon_{dataset.Reader(), "Muon_isPFcand"},
      srcIsGlobalMuon_{dataset.Reader(), "Muon_isGlobal"},
      srcIsTrackerMuon_{dataset.Reader(), "Muon_isTracker"},
      srcIdTight_{dataset.Reader(), "Muon_tightId"},
      srcTrackerLayers_{dataset.Reader(), "Muon_nTrackerLayers"} {

  if(isSim_){
    genPartId_.reset(new TTreeReaderArray<int>(dataset.Reader(), "GenPart_pdgId"));
    genPartPt_.reset(new TTreeReaderArray<float>(dataset.Reader(), "GenPart_pt"));
    genPartEta_.reset(new TTreeReaderArray<float>(dataset.Reader(), "GenPart_eta"));
    genPartPhi_.reset(new TTreeReaderArray<float>(dataset.Reader(), "GenPart_phi"));
  }
  rochesterCorrection_.reset(new RoccoR(FileInPath::Resolve("rcdata.2016.v3")));
}


std::vector<Muon> const &MuonBuilder::GetLoose() const {
  Update();
  return looseMuons_;
}


std::vector<Muon> const &MuonBuilder::GetTight() const {
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
  
  
  muon->p4.SetPtEtaPhiM(muon->p4.Pt() * scaleFactor, muon->p4.Eta(),
                        muon->p4.Phi(), muon->p4.M());
}


void MuonBuilder::Build() const {

  looseMuons_.clear();
  tightMuons_.clear();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    // Loose ID as per https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Loose_Muon
    bool const passLooseId = srcIsPfMuon_[i] && srcIsGlobalMuon_[i] && srcIsTrackerMuon_[i];

    if (std::abs(srcEta_[i]) > 2.4 or not passLooseId or
        srcIsolation_[i] > 0.25)
      continue;

    Muon muon;
    muon.p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], srcMass_[i]);
    muon.uncorrP4 = muon.p4;
    muon.charge = srcCharge_[i];

    ApplyRochesterCorrection(&muon, srcTrackerLayers_[i]);

    if (muon.p4.Pt() < minPtLoose_)
      continue;

    looseMuons_.emplace_back(muon);

    // Propagate changes in momenta of loose muons into ptmiss
    AddMomentumShift(muon.uncorrP4, muon.p4);

    bool const passTightId = srcIdTight_[i];

    if (muon.p4.Pt() < minPtTight_ or not passTightId or
        srcIsolation_[i] > 0.15)
      continue;

    tightMuons_.emplace_back(muon);
  }


  // Make sure the collections are sorted in pt
  std::sort(looseMuons_.begin(), looseMuons_.end(), PtOrdered);
  std::sort(tightMuons_.begin(), tightMuons_.end(), PtOrdered);
}


std::optional<GenParticle> MuonBuilder::FindGenMatch(
    Muon const &muon, double maxDR) const {

  unsigned iClosest = -1;
  double minDR2 = std::pow(maxDR, 2);

  for (unsigned i = 0; i < genPartId_->GetSize(); ++i) {
    if (std::abs(genPartId_->At(i)) != 13)
      // Only consider muons
      continue;

    double const dR2 = utils::DeltaR2(
      muon.p4.Eta(), muon.p4.Phi(), genPartEta_->At(i), genPartPhi_->At(i));

    if (dR2 < minDR2) {
      iClosest = i;
      minDR2 = dR2;
    }
  }

  if (iClosest != unsigned(-1)) {
    GenParticle matchedParticle{genPartId_->At(iClosest)};
    matchedParticle.p4.SetPtEtaPhiM(
      genPartPt_->At(iClosest), genPartEta_->At(iClosest), genPartPhi_->At(iClosest),
      0.1057
    );
    return matchedParticle;
  } else
    return {};
} 

