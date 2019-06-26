#include <JetBuilder.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <FileInPath.h>
#include "JERC/JetCorrectionUncertainty.h"
#include "JERC/JetResolution.h"
#include <Utils.h>


JetBuilder::JetBuilder(Dataset &dataset, Options const &options,
                       TRandom &randomGenerator)
    : CollectionBuilder{dataset.Reader()}, genJetBuilder_{nullptr},
      minPt_{30.}, maxAbsEta_{4.7}, isSim_{dataset.Info().IsSimulation()},
      syst_{Syst::None}, randomGenerator_{randomGenerator},
      srcPt_{dataset.Reader(), "JetAk04Pt"},
      srcEta_{dataset.Reader(), "JetAk04Eta"},
      srcPhi_{dataset.Reader(), "JetAk04Phi"},
      srcE_{dataset.Reader(), "JetAk04E"},
      srcBTagCsvV2_{dataset.Reader(), "JetAk04BDiscCisvV2"},
      srcHadronFlavour_{dataset.Reader(), "JetAk04HadFlav"},
      srcChf_{dataset.Reader(), "JetAk04ChHadFrac"},
      srcNhf_{dataset.Reader(), "JetAk04NeutralHadAndHfFrac"},
      srcCemf_{dataset.Reader(), "JetAk04ChEmFrac"},
      srcNemf_{dataset.Reader(), "JetAk04NeutralEmFrac"},
      srcNumConstituents_{dataset.Reader(), "JetAk04ConstCnt"},
      srcChargedMult_{dataset.Reader(), "JetAk04ChMult"},
      srcNeutralMult_{dataset.Reader(), "JetAk04NeutMult"},
      puRho_{dataset.Reader(), "EvtFastJetRho"} {

  if (isSim_) {
    jerProvider_.reset(new JME::JetResolution(FileInPath::Resolve(
      "JERC/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt")));
    jerSFProvider_.reset(new JME::JetResolutionScaleFactor(FileInPath::Resolve(
      "JERC/Summer16_25nsV1_MC_SF_AK4PFchs.txt")));

    std::string const systLabel{options.GetAs<std::string>("syst")};

    if (systLabel == "jec_up") {
      syst_ = Syst::JEC;
      systDirection_ = SystDirection::Up;
    } else if (systLabel == "jec_down") {
      syst_ = Syst::JEC;
      systDirection_ = SystDirection::Down;
    } else if (systLabel == "jer_up") {
      syst_ = Syst::JER;
      systDirection_ = SystDirection::Up;
    } else if (systLabel == "jer_down") {
      syst_ = Syst::JER;
      systDirection_ = SystDirection::Down;
    }

    if (syst_ == Syst::JEC)
      jecUncProvider_.reset(new JetCorrectionUncertainty(FileInPath::Resolve(
        "JERC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt")));
  }
}


JetBuilder::~JetBuilder() noexcept {
  // The destructor needs to be defined at a point where
  // JetCorrectionUncertainty is a complete class so that std::unique_ptr knows
  // how to destroy it
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

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    if (not PassId(i))
      continue;

    Jet jet;
    jet.p4.SetPtEtaPhiE(srcPt_[i], srcEta_[i], srcPhi_[i], srcE_[i]);
    jet.bTagCsvV2 = srcBTagCsvV2_[i];
    jet.hadronFlavour = srcHadronFlavour_[i];

    // Perform angular cleaning
    if (IsDuplicate(jet.p4, 0.4))
      continue;

    double corrFactor = 1.;
    double const corrPt = jet.p4.Pt();  // pt with nominal JEC

    // Evaluate JEC uncertainty
    if (syst_ == Syst::JEC) {
      jecUncProvider_->setJetEta(jet.p4.Eta());
      jecUncProvider_->setJetPt(corrPt);
      double const uncertainty = jecUncProvider_->getUncertainty(true);

      if (systDirection_ == SystDirection::Up)
        corrFactor *= (1. + uncertainty);
      else
        corrFactor *= (1. - uncertainty);
    }

    // Apply JER smearing. Corresponding correction factor is always evaluated
    // with nominal JEC applied, even if a JEC variation has been requested.
    // This aligns with how the JER smearing is usually applied in CMSSW.
    if (isSim_) {
      // Relative jet pt resolution in simulation
      double const ptResolution = jerProvider_->getResolution(
        {{JME::Binning::JetPt, corrPt}, {JME::Binning::JetEta, jet.p4.Eta()},
         {JME::Binning::Rho, *puRho_}});

      // Find data-to-simulation scale factor
      Variation jerDirection;

      if (syst_ == Syst::JER) {
        if (systDirection_ == SystDirection::Up)
          jerDirection = Variation::UP;
        else
          jerDirection = Variation::DOWN;
      } else
        jerDirection = Variation::NOMINAL;

      double const jerSF = jerSFProvider_->getScaleFactor(
        {{JME::Binning::JetEta, jet.p4.Eta()}}, jerDirection);


      // Depending on the presence of a matching generator-level jet, perform
      // deterministic or stochastic smearing [1]
      // [1] https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=71#Smearing_procedures
      GenJet const *genJet = FindGenMatch(jet, ptResolution);

      if (genJet) {
        double const jerFactor = 1. + 
          (jerSF - 1.) * (corrPt - genJet->p4.Pt()) / corrPt;
        corrFactor *= jerFactor;
      } else {
        double const jerFactor = 1. + randomGenerator_.Gaus(0., ptResolution) *
          std::sqrt(std::max(std::pow(jerSF, 2) - 1., 0.));
        corrFactor *= jerFactor;
      }
    }

    TLorentzVector const originalP4 = jet.p4;
    jet.p4 *= corrFactor;

    // Propagate the change in jet momentum for the use in ptmiss. The type 1
    // correction to ptmiss has been applied using jets with originalP4.Pt()
    // above 15 GeV. After the additional corrections applied above, the set of
    // jets with pt > 15 GeV has changed. However, without the access to raw
    // momenta of jets, it is not possible to undo the contribution to the
    // type 1 correction from jets whose pt has changed from above to below
    // 15 GeV. For simplicity, use the same set of jets as in the original
    // type 1 correction (modulus the different jet ID).
    if (originalP4.Pt() > 15.)
      AddMomentumShift(originalP4, jet.p4);

    // Kinematical cuts for jets to be stored in the collection
    if (jet.p4.Pt() < minPt_ or std::abs(jet.p4.Eta()) > maxAbsEta_)
      continue;

    jets_.emplace_back(jet);
  }

  // Make sure jets are sorted in pt
  std::sort(jets_.begin(), jets_.end(), PtOrdered);
}


GenJet const *JetBuilder::FindGenMatch(Jet const &jet,
                                       double ptResolution) const {
  if (not genJetBuilder_)
    return nullptr;

  // Find the closest generator-level jet in the (eta, phi) metric, with an
  // additional requirement that the difference in pt is loosely compatible with
  // the resolution. The angular distance must not exceed a half of jet radius.
  // These requirements are taken from [1].
  // [1] https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=71#Smearing_procedures
  GenJet const *match = nullptr;
  double const maxDPt = ptResolution * jet.p4.Pt() * 3;  // Use 3 sigma
  double curMinDR2 = std::pow(0.4 / 2, 2);  // Use half of jet radius

  for (auto const &genJet : genJetBuilder_->Get()) {
    double const dR2 = utils::DeltaR2(jet.p4, genJet.p4);

    if (dR2 < curMinDR2 and std::abs(jet.p4.Pt() - genJet.p4.Pt()) < maxDPt) {
      match = &genJet;
      curMinDR2 = dR2;
    }
  }

  return match;
}


bool JetBuilder::PassId(unsigned i) const {
  // Loose PF ID for 2016 [1] is implemented below
  // [1] https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2016?rev=11#Recommendations_for_the_13_TeV_d
  bool passId = false;
  double const absEta = std::abs(srcEta_[i]);

  // Multiplicities are encoded with floating-point numbers. Convert them to
  // integers safely in order to avoid rounding errors.
  int const numConstituents = int(std::round(srcNumConstituents_[i]));
  int const chargedMult = int(std::round(srcChargedMult_[i]));
  int const neutralMult = int(std::round(srcNeutralMult_[i]));

  if (absEta <= 2.7) {
    bool const commonCriteria = (srcNhf_[i] < 0.99 and srcNemf_[i] < 0.99 and
      numConstituents > 1);

    if (absEta <= 2.4)
      passId = commonCriteria and srcChf_[i] > 0. and
        chargedMult > 0 and srcCemf_[i] < 0.99;
    else
      passId = commonCriteria;
  } else if (absEta <= 3.)
    passId = neutralMult > 2 and srcNhf_[i] < 0.98 and
      srcNemf_[i] > 0.01;
  else
    passId = neutralMult > 10 and srcNemf_[i] < 0.9;

  return passId;
}

