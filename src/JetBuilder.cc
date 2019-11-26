#include <JetBuilder.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <Utils.h>


JetBuilder::JetBuilder(Dataset &dataset, Options const &options,
                       TabulatedRngEngine &rngEngine)
    : CollectionBuilder{dataset.Reader()}, genJetBuilder_{nullptr},
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
      puRho_{dataset.Reader(), "fixedGridRhoFastjetAll"},
      softRawPt_{dataset.Reader(), "CorrT1METJet_rawPt"},
      softEta_{dataset.Reader(), "CorrT1METJet_eta"},
      softPhi_{dataset.Reader(), "CorrT1METJet_phi"},
      softArea_{dataset.Reader(), "CorrT1METJet_area"} {

  auto const configNode = Options::NodeAs<YAML::Node>(
      options.GetConfig(), {"jets"});
  minPt_ = Options::NodeAs<double>(configNode, {"min_pt"});
  maxAbsEta_ = Options::NodeAs<double>(configNode, {"max_abs_eta"});

  if (isSim_)
    srcHadronFlavour_.reset(new  TTreeReaderArray<int>(
        dataset.Reader(), "Jet_hadronFlavour"));
}


std::vector<Jet> const &JetBuilder::Get() const {
  Update();
  return jets_;
}


void JetBuilder::SetGenJetBuilder(GenJetBuilder const *genJetBuilder) {
  if (isSim_)
    genJetBuilder_ = genJetBuilder;
}


void JetBuilder::Build() const {
  jets_.clear();
  jetCorrector_.UpdateIov();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    if (not srcId_[i] & (1 << 0))
      continue;

    Jet jet;
    jet.p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], srcMass_[i]);
    jet.bTag = srcBTag_[i];
    if(isSim_)
      jet.hadronFlavour = srcHadronFlavour_->At(i);
    else
      jet.hadronFlavour = 0;

    // Perform angular cleaning
    if (IsDuplicate(jet.p4, 0.4))
      continue;

    TLorentzVector const rawP4 = jet.p4 * (1 - srcRawFactor_[i]);
    double corrFactor = 1.;

    if (isSim_) {
      // Evaluate JEC uncertainty
      corrFactor *= jetCorrector_.GetJecUncFactor(jet.p4);

      // Perform JER smearing using jet four-momentum with nominal JEC applied
      double const ptResolution = jetCorrector_.GetPtResolution(jet.p4);
      GenJet const *genJet = FindGenMatch(jet.p4, ptResolution);
      corrFactor *= jetCorrector_.GetJerFactor(jet.p4, genJet, ptResolution, i);
    }

    jet.p4 *= corrFactor;

    // Type 1 correction to missing pt following the full - L1 scheme
    if (jet.p4.Pt() > 15.) {
      double const corrFactorL1 = jetCorrector_.GetJecL1(rawP4, srcArea_[i]);
      AddMomentumShift(rawP4 * corrFactorL1, jet.p4);
    }

    // Kinematical cuts for jets to be stored in the collection
    if (jet.p4.Pt() < minPt_ or std::abs(jet.p4.Eta()) > maxAbsEta_)
      continue;

    jets_.emplace_back(jet);
  }

  // Make sure jets are sorted in pt
  std::sort(jets_.begin(), jets_.end(), PtOrdered);


  // Soft jets not included into the main collection contribute to the type 1
  // correction of missing pt nonetheless. Account for them.
  for (int i = 0; i < softRawPt_.GetSize(); ++i) {
    // Jet energy is not stored, but it's not used for missing pt. Set the mass
    // to 0.
    TLorentzVector rawP4;
    rawP4.SetPtEtaPhiM(softRawPt_[i], softEta_[i], softPhi_[i], 0.);

    if (IsDuplicate(rawP4, 0.4))
      continue;

    double const area = softArea_[i];
    double const jecL1 = jetCorrector_.GetJecL1(rawP4, area);
    double const jecNominal = jetCorrector_.GetJecFull(rawP4, area);
    double corrFactorFull = jecNominal;

    if (isSim_) {
      TLorentzVector const corrP4{rawP4 * jecNominal};
      corrFactorFull *= jetCorrector_.GetJecUncFactor(corrP4);

      double const ptResolution = jetCorrector_.GetPtResolution(corrP4);
      GenJet const *genJet = FindGenMatch(corrP4, ptResolution);
      corrFactorFull *= jetCorrector_.GetJerFactor(
          corrP4, genJet, ptResolution, srcPt_.GetSize() + i);
    }

    if (rawP4.Pt() * corrFactorFull > 15.)
      AddMomentumShift(rawP4 * jecL1, rawP4 * corrFactorFull);
  }
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
