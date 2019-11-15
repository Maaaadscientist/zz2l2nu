#include <JetBuilder.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <FileInPath.h>
#include "JERC/JetCorrectionUncertainty.h"
#include "JERC/JetResolution.h"
#include <Utils.h>


JetBuilder::JetBuilder(Dataset &dataset, Options const &options,
                       TabulatedRngEngine &rngEngine)
    : CollectionBuilder{dataset.Reader()}, genJetBuilder_{nullptr},
      isSim_{dataset.Info().IsSimulation()},
      syst_{Syst::None},
      jetCorrector_{dataset, options},
      tabulatedRng_{rngEngine, 20},  // Book 20 channels
      srcPt_{dataset.Reader(), "Jet_pt"},
      srcEta_{dataset.Reader(), "Jet_eta"},
      srcPhi_{dataset.Reader(), "Jet_phi"},
      srcMass_{dataset.Reader(), "Jet_mass"},
      srcArea_{dataset.Reader(), "Jet_area"},
      srcRawFactor_{dataset.Reader(), "Jet_rawFactor"},
      srcBTag_{dataset.Reader(), (Options::NodeAs<std::string>(
        options.GetConfig(), {"b_tagger", "branch_name"})).c_str()},
      srcId_{dataset.Reader(), "Jet_jetId"},
      puRho_{dataset.Reader(), "fixedGridRhoFastjetAll"} {

  auto const configNode = Options::NodeAs<YAML::Node>(
      options.GetConfig(), {"jets"});
  minPt_ = Options::NodeAs<double>(configNode, {"min_pt"});
  maxAbsEta_ = Options::NodeAs<double>(configNode, {"max_abs_eta"});

  if (isSim_) {
    srcHadronFlavour_.reset(new  TTreeReaderArray<int>(
        dataset.Reader(), "Jet_hadronFlavour"));
    jerProvider_.reset(new JME::JetResolution(FileInPath::Resolve(
        Options::NodeAs<std::string>(configNode,
                                     {"resolution", "sim_resolution"}))));
    jerSFProvider_.reset(new JME::JetResolutionScaleFactor(FileInPath::Resolve(
        Options::NodeAs<std::string>(configNode,
                                     {"resolution", "scale_factors"}))));

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
          Options::NodeAs<std::string>(configNode,
                                       {"corrections", "uncertainty"}))));
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
      corrFactor *= ComputeJecUncFactor(jet.p4);

      // Perform JER smearing using jet four-momentum with nominal JEC applied
      corrFactor *= ComputeJerFactor(jet.p4, *puRho_, i);
    }

    jet.p4 *= corrFactor;

    // Type 1 correction to missing pt following the full - L1 scheme
    if (jet.p4.Pt() > 15.) {
      double const jecL1 = jetCorrector_.GetJecL1(rawP4, srcArea_[i]);
      AddMomentumShift(rawP4 * jecL1, jet.p4);
    }

    // Kinematical cuts for jets to be stored in the collection
    if (jet.p4.Pt() < minPt_ or std::abs(jet.p4.Eta()) > maxAbsEta_)
      continue;

    jets_.emplace_back(jet);
  }

  // Make sure jets are sorted in pt
  std::sort(jets_.begin(), jets_.end(), PtOrdered);
}


double JetBuilder::ComputeJecUncFactor(TLorentzVector const &corrP4) const {
  if (syst_ != Syst::JEC)
    return 1.;

  jecUncProvider_->setJetEta(corrP4.Eta());
  jecUncProvider_->setJetPt(corrP4.Pt());
  double const uncertainty = jecUncProvider_->getUncertainty(true);

  if (systDirection_ == SystDirection::Up)
    return 1. + uncertainty;
  else
    return 1. - uncertainty;
}


double JetBuilder::ComputeJerFactor(TLorentzVector const &corrP4, double rho,
                                    int rngChannel) const {
    // Relative jet pt resolution in simulation
    double const ptResolution = jerProvider_->getResolution(
        {{JME::Binning::JetPt, corrP4.Pt()},
         {JME::Binning::JetEta, corrP4.Eta()},
         {JME::Binning::Rho, rho}});

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
        {{JME::Binning::JetEta, corrP4.Eta()}}, jerDirection);

    // Depending on the presence of a matching generator-level jet, perform
    // deterministic or stochastic smearing [1]
    // [1] https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=71#Smearing_procedures
    GenJet const *genJet = FindGenMatch(corrP4, ptResolution);

    if (genJet) {
      double const jerFactor = 1.
          + (jerSF - 1.) * (corrP4.Pt() - genJet->p4.Pt()) / corrP4.Pt();
      return jerFactor;
    } else {
      double const jerFactor = 1.
          + tabulatedRng_.Gaus(rngChannel, 0., ptResolution)
          * std::sqrt(std::max(std::pow(jerSF, 2) - 1., 0.));
      return jerFactor;
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
