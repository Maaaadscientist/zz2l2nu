#include <JetBuilder.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <FileInPath.h>
#include "JERC/JetCorrectionUncertainty.h"


JetBuilder::JetBuilder(TTreeReader &reader, Options const &options)
    : minPt_{30.}, maxAbsEta_{4.7}, cache_{reader},
      syst_{Syst::None},
      srcPt_{reader, "JetAk04Pt"}, srcEta_{reader, "JetAk04Eta"},
      srcPhi_{reader, "JetAk04Phi"}, srcE_{reader, "JetAk04E"},
      srcBTagCsvV2_{reader, "JetAk04BDiscCisvV2"},
      srcHadronFlavour_{reader, "JetAk04HadFlav"},
      srcChf_{reader, "JetAk04ChHadFrac"},
      srcNhf_{reader, "JetAk04NeutralHadAndHfFrac"},
      srcCemf_{reader, "JetAk04ChEmFrac"},
      srcNemf_{reader, "JetAk04NeutralEmFrac"},
      srcNumConstituents_{reader, "JetAk04ConstCnt"},
      srcChargedMult_{reader, "JetAk04ChMult"},
      srcNeutralMult_{reader, "JetAk04NeutMult"} {

  if (options.GetAs<bool>("is-mc")) {
    std::string const systLabel{options.GetAs<std::string>("syst")};

    if (systLabel == "jec_up") {
      syst_ = Syst::JEC;
      systDirection_ = SystDirection::Up;
    } else if (systLabel == "jec_down") {
      syst_ = Syst::JEC;
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


void JetBuilder::EnableCleaning(
    std::initializer_list<CollectionBuilder const *> builders) {
  for (auto *b : builders)
    buildersForCleaning_.emplace_back(b);
}


std::vector<Jet> const &JetBuilder::Get() const {
  if (cache_.IsUpdated())
    Build();

  return jets_;
}


void JetBuilder::Build() const {
  jets_.clear();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    if (std::abs(srcEta_[i]) > maxAbsEta_ or not PassId(i))
      continue;

    Jet jet;
    jet.p4.SetPtEtaPhiE(srcPt_[i], srcEta_[i], srcPhi_[i], srcE_[i]);
    jet.bTagCsvV2 = srcBTagCsvV2_[i];
    jet.hadronFlavour = srcHadronFlavour_[i];

    // Perform angular cleaning
    if (IsDuplicate(jet))
      continue;

    double corrFactor = 1.;

    if (syst_ == Syst::JEC) {
      jecUncProvider_->setJetEta(jet.p4.Eta());
      jecUncProvider_->setJetPt(jet.p4.Pt());
      double const uncertainty = jecUncProvider_->getUncertainty(true);

      if (systDirection_ == SystDirection::Up)
        corrFactor *= (1. + uncertainty);
      else
        corrFactor *= (1. - uncertainty);
    }

    jet.p4 *= corrFactor;

    if (jet.p4.Pt() < minPt_)
      continue;

    jets_.emplace_back(jet);
  }

  // Make sure jets are sorted in pt
  std::sort(jets_.begin(), jets_.end(), PtOrdered);
}



bool JetBuilder::IsDuplicate(Jet const &jet) const {
  for (auto const builder : buildersForCleaning_) {
    if (builder->GetMomenta().HasOverlap(jet.p4, 0.4))
      return true;
  }

  return false;
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

