#include <PileUpIdWeight.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <map>
#include <string>

#include <TFile.h>

#include <FileInPath.h>


/**
 * \brief Auxiliary class to simplify reading of histograms in
 * PileUpIdWeight::LoadScaleFactors
 *
 * It reads 2D histograms for a ROOT file, wrapping them in std::shared_ptr. If
 * a histogram with the same name is requested again, the already existing
 * std::shared_ptr is returned.
 */
class HistReader {
 public:
  HistReader(std::filesystem::path const &path);
  ~HistReader();
  std::shared_ptr<TH2> operator()(std::string const &name);

 private:
  TFile file_;
  std::map<std::string, std::shared_ptr<TH2>> readHistograms_;
};

HistReader::HistReader(std::filesystem::path const &path)
    : file_{path.c_str()} {}

HistReader::~HistReader() {
  file_.Close();
}

std::shared_ptr<TH2> HistReader::operator()(std::string const &name) {
  auto const res = readHistograms_.find(name);
  if (res != readHistograms_.end())
    return res->second;
  auto hist = std::shared_ptr<TH2>(file_.Get<TH2>(name.c_str()));
  hist->SetDirectory(nullptr);
  readHistograms_[name] = hist;
  return hist;
}


PileUpIdWeight::PileUpIdWeight(
    Dataset &dataset, Options const &options,
    PileUpIdFilter const *pileUpIdFilter,
    JetBuilder const *jetBuilder)
    : pileUpIdFilter_{pileUpIdFilter}, jetBuilder_{jetBuilder},
      absEtaEdges_{pileUpIdFilter_->GetAbsEtaEdges()},
      expPileUp_{dataset.Reader(), "Pileup_nTrueInt"},
      cache_{dataset.Reader()} {
  for (auto const &wp : pileUpIdFilter_->GetWorkingPoints())
    contexts_.emplace_back(wp);

  int const year = Options::NodeAs<int>(options.GetConfig(), {"period"});
  LoadScaleFactors(options.GetConfig(), year);

  auto const effModelPath = FileInPath::Resolve(
      Options::NodeAs<std::string>(
          options.GetConfig(), {"pileup_id", "efficiency"}));
  effCalc_.emplace(effModelPath, effFeatures_.size());

  // The year won't change, so the corresponding features can be set already
  // now
  effFeatures_[4] = (year == 2016) ? 1 : 0;
  effFeatures_[5] = (year == 2017) ? 1 : 0;
  effFeatures_[6] = (year == 2018) ? 1 : 0;

  auto const systLabel = options.GetAs<std::string>("syst");
  if (systLabel == "puid_tag_up")
    defaultVariation_ = Variation::kTagUp;
  else if (systLabel == "puid_tag_down")
    defaultVariation_ = Variation::kTagDown;
  else if (systLabel == "puid_mistag_up")
    defaultVariation_ = Variation::kMistagUp;
  else if (systLabel == "puid_mistag_down")
    defaultVariation_ = Variation::kMistagDown;
  else
    defaultVariation_ = Variation::kNominal;
}


std::string_view PileUpIdWeight::VariationName(int variation) const {
  switch (variation) {
    case 0:
      return "puid_tag_up";
    case 1:
      return "puid_tag_down";
    case 2:
      return "puid_mistag_up";
    case 3:
      return "puid_mistag_down";
    default:
      return "";
  }
}


PileUpIdWeight::Context const &PileUpIdWeight::FindContext(
    Jet const &jet) const {
  int const bin = std::upper_bound(
      absEtaEdges_.begin(), absEtaEdges_.end(), std::abs(jet.p4.Eta()))
      - absEtaEdges_.begin();
  return contexts_[bin];
}


double PileUpIdWeight::GetEfficiency(
    Context const &context, Jet const &jet) const {
  int const f = jet.combFlavour;
  effFeatures_[0] = jet.p4.Pt();
  effFeatures_[1] = jet.p4.Eta();
  effFeatures_[2] = *expPileUp_;
  effFeatures_[3] = int(context.workingPoint);
  // Indices 4 to 6 correspond to year and has been already set
  effFeatures_[7] = (jet.isPileUp) ? 1 : 0;
  effFeatures_[8] = (f == 21 or f == 0) ? 1 : 0;
  effFeatures_[9] = (f == 1 or f == 2 or f == 3) ? 1 : 0;
  effFeatures_[10] = (f == 4) ? 1 : 0;
  effFeatures_[11] = (f == 5) ? 1 : 0;
  effFeatures_[12] = std::abs(jet.p4.Eta());
  return effCalc_->Predict(effFeatures_.data());
}


double PileUpIdWeight::GetScaleFactor(
    Context const &context, Jet const &jet, Variation variation) const {
  std::shared_ptr<TH2> histValue, histUnc;
  if (jet.isPileUp) {
    histValue = context.sfPileUp;
    histUnc = context.sfUncPileUp;
  } else {
    histValue = context.sfMatched;
    histUnc = context.sfUncMatched;
  }

  int const bin = histValue->FindFixBin(jet.p4.Pt(), jet.p4.Eta());
  double const sfNominal = histValue->GetBinContent(bin);

  int shift = 0;
  if (jet.isPileUp) {
    if (variation == Variation::kMistagUp)
      shift = +1;
    else if (variation == Variation::kMistagDown)
      shift = -1;
  } else {
    if (variation == Variation::kTagUp)
      shift = +1;
    else if (variation == Variation::kTagDown)
      shift = -1;
  }

  if (shift == 0) {
    return sfNominal;
  } else {
    int const bin = histUnc->FindFixBin(jet.p4.Pt(), jet.p4.Eta());
    double const sfUnc = histUnc->GetBinContent(bin);
    double sf = sfNominal + shift * sfUnc;
    if (sf < 0.)
      sf = 0.;
    else if (sf > 5.)
      sf = 5.;
    return sf;
  }
}


void PileUpIdWeight::LoadScaleFactors(YAML::Node const config, int year) {
  std::filesystem::path const path = FileInPath::Resolve(
      Options::NodeAs<std::string>(config, {"pileup_id", "scale_factors"}));
  HistReader histReader{path};
  for (auto &context : contexts_) {
    if (context.workingPoint == Jet::PileUpId::None)
      continue;

    std::string wpLabel;
    switch (context.workingPoint) {
      case Jet::PileUpId::Loose:
        wpLabel = "L";
        break;
      case Jet::PileUpId::Medium:
        wpLabel = "M";
        break;
      case Jet::PileUpId::Tight:
        wpLabel = "T";
        break;
      default:
        wpLabel = "";
    }
    std::string const nameFragment = std::to_string(year) + "_" + wpLabel;

    context.sfMatched = histReader("h2_eff_sf" + nameFragment);
    context.sfUncMatched = histReader(
        "h2_eff_sf" + nameFragment + "_Systuncty");
    context.sfPileUp = histReader("h2_mistag_sf" + nameFragment);
    context.sfUncPileUp = histReader(
        "h2_mistag_sf" + nameFragment + "_Systuncty");
  }
}


void PileUpIdWeight::Update() const {
  std::fill(weights_.begin(), weights_.end(), 1.);

  // Jets that pass pileup ID
  for (auto const &jet : jetBuilder_->Get()) {
    if (not pileUpIdFilter_->IsTaggable(jet))
      continue;
    auto const &context = FindContext(jet);
    if (context.workingPoint == Jet::PileUpId::None)
      continue;

    double const eff = GetEfficiency(context, jet);
    for (int iVar = 0; iVar < 5; ++iVar) {
      double const sf = GetScaleFactor(context, jet, Variation(iVar));
      weights_[iVar] *= std::min(sf * eff, 1.) / eff;
    }
  }

  // Jets that fail pileup ID
  for (auto const &jet : jetBuilder_->GetRejected()) {
    if (not pileUpIdFilter_->IsTaggable(jet))
      continue;
    auto const &context = FindContext(jet);
    if (context.workingPoint == Jet::PileUpId::None)
      continue;

    double const eff = GetEfficiency(context, jet);
    for (int iVar = 0; iVar < 5; ++iVar) {
      double const sf = GetScaleFactor(context, jet, Variation(iVar));
      weights_[iVar] *= std::max(1. - sf * eff, 0.) / (1. - eff);
    }
  }

  LOG_TRACE << "Nominal weight in PileUpIdWeight: " << weights_[0] << ".";
}

