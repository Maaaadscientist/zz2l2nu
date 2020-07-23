#include <JetBuilder.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <tuple>

#include <Utils.h>


JetBuilder::JetBuilder(
    Dataset &dataset, Options const &options, TabulatedRngEngine &rngEngine,
    PileUpIdFilter const *pileUpIdFilter)
    : CollectionBuilder{dataset.Reader()},
      genJetBuilder_{nullptr}, pileUpIdFilter_{pileUpIdFilter},
      minPtType1Corr_{15.}, applyEeNoiseMitigation_{false},
      isSim_{dataset.Info().IsSimulation()},
      jetCorrector_{dataset, options, rngEngine},
      srcPt_{dataset.Reader(), "Jet_pt"},
      srcEta_{dataset.Reader(), "Jet_eta"},
      srcPhi_{dataset.Reader(), "Jet_phi"},
      srcMass_{dataset.Reader(), "Jet_mass"},
      srcArea_{dataset.Reader(), "Jet_area"},
      srcRawFactor_{dataset.Reader(), "Jet_rawFactor"},
      srcBTag_{dataset.Reader(), (Options::NodeAs<std::string>(
        options.GetConfig(), {"b_tagger", "branch_name"})).c_str()},
      srcId_{dataset.Reader(), "Jet_jetId"},
      srcPileUpId_{dataset.Reader(), "Jet_puId"},
      puRho_{dataset.Reader(), "fixedGridRhoFastjetAll"},
      softRawPt_{dataset.Reader(), "CorrT1METJet_rawPt"},
      softEta_{dataset.Reader(), "CorrT1METJet_eta"},
      softPhi_{dataset.Reader(), "CorrT1METJet_phi"},
      softArea_{dataset.Reader(), "CorrT1METJet_area"} {

  auto const configNode = Options::NodeAs<YAML::Node>(
      options.GetConfig(), {"jets"});
  jetIdBit_ = Options::NodeAs<int>(configNode, {"jet_id_bit"});
  minPt_ = Options::NodeAs<double>(configNode, {"min_pt"});
  maxAbsEta_ = Options::NodeAs<double>(configNode, {"max_abs_eta"});

  if (pileUpIdFilter_)
    LOG_DEBUG << "PileUpIdFilter is registered in JetBuilder.";
  if (pileUpIdFilter_) {
    std::tie(pileUpIdMinPt_, pileUpIdMaxPt_) = pileUpIdFilter_->GetPtRange();
  } else {
    // Defaults from https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID?rev=61#Recommendations_for_13_TeV_data
    pileUpIdMinPt_ = 15.;
    pileUpIdMaxPt_ = 50.;
  }

  if (auto const node = options.GetConfig()["ptmiss_fix_ee_2017"]; node)
    applyEeNoiseMitigation_ = node.as<bool>();

  if (isSim_) {
    srcHadronFlavour_.emplace(dataset.Reader(), "Jet_hadronFlavour");
    srcPartonFlavour_.emplace(dataset.Reader(), "Jet_partonFlavour");
    srcGenJetIdx_.emplace(dataset.Reader(), "Jet_genJetIdx");
  }
}


std::vector<Jet> const &JetBuilder::Get() const {
  Update();
  return jets_;
}


std::vector<Jet> const &JetBuilder::GetRejected() const {
  Update();
  return rejectedJets_;
}


void JetBuilder::SetGenJetBuilder(GenJetBuilder const *genJetBuilder) {
  if (isSim_)
    genJetBuilder_ = genJetBuilder;
}


void JetBuilder::AddType1Correction(
    TLorentzVector const &rawP4, double area,
    double corrFactorOrig, double corrFactorNew) const {
  // Check if this jet should be subjected to the EE noise mitigation procedure
  bool isEeNoise = false;
  if (applyEeNoiseMitigation_) {
    double const absEta = std::abs(rawP4.Eta());
    // Thresholds are from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription?rev=91#Instructions_for_2017_data_with
    isEeNoise = (rawP4.Pt() < 50. and absEta > 2.65 and absEta < 3.139);
  }

  // The normal type 1 correction for jets not affected by the EE noise
  double const jecL1 = jetCorrector_.GetJecL1(rawP4, area);
  if (not isEeNoise and rawP4.Pt() * corrFactorNew > minPtType1Corr_)
    AddMomentumShift(rawP4 * jecL1, rawP4 * corrFactorNew);

  // If the EE noise mitigation is enabled, the computation starts from an
  // adjusted ptmiss instead of the raw one, and some of the contributions in it
  // have to be removed. Technically, the removal procedure is equivalent to
  // applying the type 1 correction with affected jets but using the same JEC as
  // when NanoAOD was produced. See [1] for details.
  // [1] https://hypernews.cern.ch/HyperNews/CMS/get/met/710.html
  if (isEeNoise and rawP4.Pt() * corrFactorOrig > minPtType1Corr_)
    AddMomentumShift(rawP4 * jecL1, rawP4 * corrFactorOrig);
}


void JetBuilder::Build() const {
  jetCorrector_.UpdateIov();
  ProcessJets();

  // Soft jets not included into the main collection contribute to the type 1
  // correction of missing pt nonetheless. Account for them.
  ProcessSoftJets();
}


GenJet const *JetBuilder::FindGenMatch(TLorentzVector const &p4,
                                       double ptResolution) const {
  if (not genJetBuilder_)
    return nullptr;

  // Find the closest generator-level jet in the (eta, phi) metric, with an
  // additional requirement that the difference in pt is loosely compatible with
  // the resolution. The angular distance must not exceed a half of jet radius.
  // These requirements are taken from [1].
  // [1] https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=71#Smearing_procedures
  GenJet const *match = nullptr;
  double const maxDPt = ptResolution * p4.Pt() * 3;  // Use 3 sigma
  double curMinDR2 = std::pow(0.4 / 2, 2);  // Use half of jet radius

  for (auto const &genJet : genJetBuilder_->Get()) {
    double const dR2 = utils::DeltaR2(p4, genJet.p4);

    if (dR2 < curMinDR2 and std::abs(p4.Pt() - genJet.p4.Pt()) < maxDPt) {
      match = &genJet;
      curMinDR2 = dR2;
    }
  }

  return match;
}


double JetBuilder::GetJerFactor(
    TLorentzVector const &corrP4, int rngChannel) const {
  double const ptResolution = jetCorrector_.GetPtResolution(corrP4);
  GenJet const *genJet = FindGenMatch(corrP4, ptResolution);
  return jetCorrector_.GetJerFactor(
      corrP4, genJet, ptResolution, rngChannel);
}


void JetBuilder::ProcessJets() const {
  jets_.clear();
  rejectedJets_.clear();

 for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    if (not (srcId_[i] & 1 << jetIdBit_))
      continue;

    Jet jet;
    jet.p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], srcMass_[i]);

    // Perform angular cleaning
    if (IsDuplicate(jet.p4, 0.4))
      continue;

    double const jecNominal = 1. / (1 - srcRawFactor_[i]);
    TLorentzVector const rawP4 = jet.p4  * (1. / jecNominal);
    double corrFactor = 1.;
    if (isSim_) {
      corrFactor *= jetCorrector_.GetJecUncFactor(jet.p4);
      corrFactor *= GetJerFactor(jet.p4, i);
    }
    jet.p4 *= corrFactor;

    AddType1Correction(rawP4, srcArea_[i], jecNominal, jecNominal * corrFactor);

    // Kinematical cuts for jets to be stored in the collection
    if (jet.p4.Pt() < minPt_ or std::abs(jet.p4.Eta()) > maxAbsEta_)
      continue;

    jet.bTag = srcBTag_[i];
    if (isSim_)
      jet.SetFlavours(srcHadronFlavour_->At(i), srcPartonFlavour_->At(i));

    bool const puIdAccepted = SetPileUpInfo(jet, i);
    if (puIdAccepted)
      jets_.emplace_back(jet);
    else
      rejectedJets_.emplace_back(jet);
  }

  // Make sure jets are sorted in pt
  std::sort(jets_.begin(), jets_.end(), PtOrdered);
}


void JetBuilder::ProcessSoftJets() const {
  for (int i = 0; i < int(softRawPt_.GetSize()); ++i) {
    // Jet energy is not stored, but it's not used for missing pt. Set the mass
    // to 0.
    TLorentzVector rawP4;
    rawP4.SetPtEtaPhiM(softRawPt_[i], softEta_[i], softPhi_[i], 0.);

    if (IsDuplicate(rawP4, 0.4))
      continue;

    double const area = softArea_[i];
    double const jecNominal = jetCorrector_.GetJecFull(rawP4, area);
    double corrFactorFull = jecNominal;

    if (isSim_) {
      TLorentzVector const corrP4{rawP4 * jecNominal};
      corrFactorFull *= jetCorrector_.GetJecUncFactor(corrP4);
      corrFactorFull *= GetJerFactor(corrP4, srcPt_.GetSize() + i);
    }

    AddType1Correction(rawP4, area, jecNominal, corrFactorFull);
  }
}


bool JetBuilder::SetPileUpInfo(Jet &jet, int index) const {
  if (jet.p4.Pt() < pileUpIdMinPt_ or jet.p4.Pt() > pileUpIdMaxPt_) {
    jet.pileUpId = Jet::PileUpId::PassThrough;
  } else {
    int const id = srcPileUpId_[index];
    if (id & 1)
      jet.pileUpId = Jet::PileUpId::Tight;
    else if (id & 1 << 1)
      jet.pileUpId = Jet::PileUpId::Medium;
    else if (id & 1 << 2)
      jet.pileUpId = Jet::PileUpId::Loose;
    else
      jet.pileUpId = Jet::PileUpId::None;
  }

  if (isSim_) {
    // The standard matching [1] is done based on dR, with the cut-off dR < 0.4.
    // However, each particle-level jet can be matched to one reconstructed jet
    // at maximum [2] to guarantee an unambiguous matching.
    // [1] https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_22/PhysicsTools/PatAlgos/python/mcMatchLayer0/jetMatch_cfi.py#L18-L28
    // [2] https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_22/CommonTools/UtilAlgos/interface/PhysObjectMatcher.h#L150-L159
    jet.isPileUp = (srcGenJetIdx_->At(index) == -1);
  }

  if (pileUpIdFilter_ and not (*pileUpIdFilter_)(jet)) {
    LOG_TRACE << "Pileup ID filter rejects jets with pt " << jet.p4.Pt()
        << " GeV, eta " << jet.p4.Eta() << ", and pileup ID WP "
        << int(jet.pileUpId) << ".";
    return false;
  }

  return true;
}
