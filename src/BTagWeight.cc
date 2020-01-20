#include <BTagWeight.h>
#include <FileInPath.h>

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

#include <TFile.h>

#include "BTag/BTagCalibrationStandalone.h"

using namespace std::string_literals;


BTagWeight::BTagWeight(Dataset &dataset, Options const &options,
                       BTagger const *bTagger, JetBuilder const *jetBuilder)
    : bTagger_{bTagger}, jetBuilder_{jetBuilder},
      effProfilesPath_{FileInPath::Resolve(
        Options::NodeAs<std::string>(
          options.GetConfig(), {"b_tag_weight", "efficiency_tables"}))},
      scaleFactorReader_{new BTagCalibrationReader{
        BTagEntry::OP_LOOSE, "central", {"up", "down"}}},
      cache_{dataset.Reader()} {

  std::string const scaleFactorsPath{FileInPath::Resolve(
    Options::NodeAs<std::string>(
      options.GetConfig(), {"b_tag_weight", "scale_factors"}))};

  BTagCalibration calibration{"", scaleFactorsPath};

  scaleFactorReader_->load(calibration, BTagEntry::FLAV_B, "mujets");
  scaleFactorReader_->load(calibration, BTagEntry::FLAV_C, "mujets");
  scaleFactorReader_->load(calibration, BTagEntry::FLAV_UDSG, "incl");

  LoadEffProfiles();

  auto const systLabel = options.GetAs<std::string>("syst");
  if (systLabel == "btag_up")
    defaultVariation_ = Variation::kTagUp;
  else if (systLabel == "btag_down")
    defaultVariation_ = Variation::kTagDown;
  else if (systLabel == "bmistag_up")
    defaultVariation_ = Variation::kMistagUp;
  else if (systLabel == "bmistag_down")
    defaultVariation_ = Variation::kMistagDown;
  else
    defaultVariation_ = Variation::kNominal;
}


BTagWeight::~BTagWeight() noexcept {}


std::string_view BTagWeight::VariationName(int variation) const {
  switch (variation) {
    case 0:
      return "btag_up";
    case 1:
      return "btag_down";
    case 2:
      return "bmistag_up";
    case 3:
      return "bmistag_down";
    default:
      return "";
  }
}


double BTagWeight::ComputeWeight(Variation variation) const {
  double weight = 1.;

  for (auto const &jet : jetBuilder_->Get()) {
    if (not bTagger_->IsTaggable(jet))
      continue;

    double const sf = GetScaleFactor(
      jet.p4.Pt(), jet.p4.Eta(), jet.hadronFlavour, variation);

    if ((*bTagger_)(jet)) {
      weight *= sf;
    } else {
      double const eff = GetEfficiency(
          jet.p4.Pt(), jet.p4.Eta(), jet.hadronFlavour);
      weight *= (1 - sf * eff) / (1 - eff);
    }
  }

  return weight;
}


double BTagWeight::GetEfficiency(double pt, double eta, int flavour) const {

  std::string flavourLabel;
  eta = std::fabs(eta);

  switch (std::abs(flavour)) {
    case 5:
      flavourLabel = "b";
      break;

    case 4:
      flavourLabel = "c";
      break;

    default:
      flavourLabel = "udsg";
  }
  int const globalBin = effProfiles_.at(flavourLabel)->FindFixBin(pt, eta);
  return effProfiles_.at(flavourLabel)->GetEfficiency(globalBin);
}


void BTagWeight::LoadEffProfiles() {

  TFile inputFile{effProfilesPath_.c_str()};

  if (inputFile.IsZombie()) {
    std::ostringstream message;
    message << "Could not open file " << effProfilesPath_ << ".";
    throw std::runtime_error(message.str());
  }

  for(std::string const &flavor : {"b", "c", "udsg"}) { 
    effProfiles_[flavor].reset((TEfficiency *) inputFile.Get(flavor.c_str()));

    if (not effProfiles_[flavor]) {
      std::ostringstream message;
      message << "File " << effProfilesPath_ <<
        " does not contain required histogram \"" << flavor << "\".";
      throw std::runtime_error(message.str());
    }

    effProfiles_[flavor]->SetDirectory(nullptr);
  }
  inputFile.Close();
}


double BTagWeight::GetScaleFactor(double pt, double eta, int flavour,
                                  Variation variation) const {
  BTagEntry::JetFlavor translatedFlavour;
  switch (std::abs(flavour)) {
    case 5:
      translatedFlavour = BTagEntry::FLAV_B;
      break;
    case 4:
      translatedFlavour = BTagEntry::FLAV_C;
      break;
    default:
      translatedFlavour = BTagEntry::FLAV_UDSG;
  }

  std::string version{"central"};
  if (translatedFlavour != BTagEntry::FLAV_UDSG) {
    if (variation == Variation::kTagUp)
      version = "up";
    else if (variation == Variation::kTagDown)
      version = "down";
  } else {
    if (variation == Variation::kMistagUp)
      version = "up";
    else if (variation == Variation::kMistagDown)
      version = "down";
  }

  return scaleFactorReader_->eval_auto_bounds(
    version, translatedFlavour, eta, pt);
}


void BTagWeight::Update() const {
  for (int i = 0; i < int(weights_.size()); ++i)
    weights_[i] = ComputeWeight(Variation(i));
}

