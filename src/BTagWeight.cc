#include <BTagWeight.h>
#include <FileInPath.h>

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

#include "BTag/BTagCalibrationStandalone.h"

using namespace std::string_literals;


BTagWeight::BTagWeight(Dataset &dataset, Options const &options,
                       BTagger const *bTagger, JetBuilder const *jetBuilder)
    : bTagger_{bTagger}, jetBuilder_{jetBuilder},
      effTablePath_{Options::NodeAs<std::string>(
        options.GetConfig(), {"b_tag_weight", "efficiency_tables"})},
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

  LoadEffTables();

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

  auto const &table = efficiencyTables_.at(
    flavourLabel);
  return table.getEfficiency(pt, eta);
}


void BTagWeight::LoadEffTables()
{
  size_t pos = effTablePath_.find("{");
  size_t len = effTablePath_.find("}") - pos + 1;
  
  for(auto const &flavour : {"b", "c", "udsg"}) {
    std::string effTablePath = effTablePath_;
    effTablePath.replace(pos, len, flavour);
    efficiencyTables_.insert(
      std::pair<std::string, utils::table>(
        flavour,
        utils::table(FileInPath::Resolve(effTablePath))));
  }
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

