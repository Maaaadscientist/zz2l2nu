#include <BTagWeight.h>

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>


using namespace std::string_literals;


BTagWeight::BTagWeight(Options const &options)
    : scaleFactorReader_{bTagWorkingPoint_, "central", {"up", "down"}},
      syst_{options.GetAs<std::string>("syst")} {

  std::string const scaleFactorsPath{std::getenv("HZZ2L2NU_BASE") +
    "/data/efficiencyTables/CSVv2_Moriond17_B_H.csv"s};
  BTagCalibration calibration{bTagAlgorithm_, scaleFactorsPath};

  scaleFactorReader_.load(calibration, BTagEntry::FLAV_B, "mujets");
  scaleFactorReader_.load(calibration, BTagEntry::FLAV_C, "mujets");
  scaleFactorReader_.load(calibration, BTagEntry::FLAV_UDSG, "incl");

  LoadEffTables();
}


double BTagWeight::operator()(std::vector<Jet> const &jets) const {

  double weight = 1.;

  for (auto const &jet : jets) {
    if (std::abs(jet.p4.Eta()) > 2.5)
      continue;

    double const sf = GetScaleFactor(
      jet.p4.Pt(), jet.p4.Eta(), jet.hadronFlavour);

    if (jet.bTagCsvV2 > bTagCut_) {
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

  std::string const wpLabel{WPToText(bTagWorkingPoint_)};

  auto const &table = efficiencyTables_.at(wpLabel + "-" + flavourLabel);
  return table.getEfficiency(pt, eta);
}


void BTagWeight::LoadEffTables()
{
  std::string const effTablePath = std::getenv("HZZ2L2NU_BASE") +
    "/data/efficiencyTables/"s;
  
  for (auto const &wp : {WPToText(bTagWorkingPoint_)}) {
    for(auto const &flavour : {"b", "c", "udsg"}) {
      efficiencyTables_.insert(
        std::pair<std::string, utils::table>(
          std::string{wp} + "-" + flavour,
          utils::table(effTablePath + "btag-" + wp + "-" + flavour + ".txt")));
    }
  }
}


double BTagWeight::GetScaleFactor(double pt, double eta, int flavour) const {

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

  if (syst_ == "btag_up")
    version = "up";
  else if (syst_ == "btag_down")
    version = "down";

  return scaleFactorReader_.eval_auto_bounds(
    version, translatedFlavour, std::abs(eta), pt);
}


std::string BTagWeight::WPToText(BTagEntry::OperatingPoint wp) {
  
  switch (wp) {
    case BTagEntry::OP_LOOSE:
      return "loose";

    case BTagEntry::OP_MEDIUM:
      return "medium";

    case BTagEntry::OP_TIGHT:
      return "tight";

    default: {
      std::ostringstream message;
      message << "Unknown b tag working point " << wp << "\"";
      throw std::runtime_error(message.str());
    }
  }
}

