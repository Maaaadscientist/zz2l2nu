#include <GenWeight.h>

#include <algorithm>
#include <initializer_list>
#include <regex>
#include <string>

#include <yaml-cpp/yaml.h>

#include <FileInPath.h>
#include <HZZException.h>
#include <Logger.h>


GenWeight::GenWeight(Dataset &dataset, Options const &options)
  : srcLheNominalWeight_{dataset.Reader(), "LHEWeight_originalXWGTUP"},
    srcGenNominalWeight_{dataset.Reader(), "Generator_weight"},
    srcScaleWeights_{dataset.Reader(), "LHEScaleWeight"},
    srcPdfWeights_{dataset.Reader(), "LHEPdfWeight"} {

  DatasetInfo const &info = dataset.Info();
  datasetWeight_ = info.CrossSection()
      / (info.NumEventsTotal() * info.MeanWeight());

  InitializeLheScale(dataset);
  InitializePdf(dataset);

  if (lheScaleWeightsPresent_)
    for (auto const &name : {
        "me_renorm_up", "me_renorm_down", "factor_up", "factor_down"})
      availableVariations_.emplace_back(name);
  if (pdfWeightsPresent_) {
    availableVariations_.emplace_back("pdf_up");
    availableVariations_.emplace_back("pdf_down");
  }
  if (alphaSWeightsPresent_) {
    availableVariations_.emplace_back("alphaS_up");
    availableVariations_.emplace_back("alphaS_down");
  }

  auto const systLabel = options.GetAs<std::string>("syst");
  auto const r = std::find(
      availableVariations_.begin(), availableVariations_.end(), systLabel);
  if (r == availableVariations_.end())
    defaultVariationIndex_ = -1;
  else
    defaultVariationIndex_ = r - availableVariations_.begin();
}


double GenWeight::EnvelopeMEScale(Var direction) const {
  // Find weights for all variations in the two ME scales except for the cases
  // when they go in opposite directions
  std::initializer_list<double> const weights{
    RelWeightMEScale(Var::Nominal, Var::Up),
    RelWeightMEScale(Var::Nominal, Var::Down),
    RelWeightMEScale(Var::Up, Var::Nominal),
    RelWeightMEScale(Var::Down, Var::Nominal),
    RelWeightMEScale(Var::Up, Var::Up),
    RelWeightMEScale(Var::Down, Var::Down)
  };

  if (direction == Var::Up)
    return std::max(weights);
  else if (direction == Var::Down)
    return std::min(weights);
  else
    return 1.;
}


double GenWeight::NominalWeight() const {
  return *srcGenNominalWeight_ * datasetWeight_;
}


double GenWeight::RelWeight(int variation) const {
  enum class Group {
    None,
    MEScale,
    PDF,
    AlphaS
  };

  Group group = Group::None;
  int index = variation;
  for (auto const &[g, offset, present] : {
      std::make_tuple(Group::MEScale, 4, lheScaleWeightsPresent_),
      std::make_tuple(Group::PDF, 2, pdfWeightsPresent_),
      std::make_tuple(Group::AlphaS, 2, alphaSWeightsPresent_)}) {
    if (not present)
      continue;
    if (index < offset) {
      group = g;
      break;
    } else {
      index -= offset;
    }
  }
  if (group == Group::None)
    throw HZZException{"Illegal variation index."};

  if (group == Group::MEScale) {
    switch (index) {
      case 0:
        return RelWeightMEScale(Var::Up, Var::Nominal);
      case 1:
        return RelWeightMEScale(Var::Down, Var::Nominal);
      case 2:
        return RelWeightMEScale(Var::Nominal, Var::Up);
      case 3:
        return RelWeightMEScale(Var::Nominal, Var::Down);
    }
  } else if (group == Group::PDF) {
    switch (index) {
      case 0:
        return RelWeightPdf(Var::Up);
      case 1:
        return RelWeightPdf(Var::Down);
    }
  } else if (group == Group::AlphaS) {
    switch (index) {
      case 0:
        return RelWeightAlphaS(Var::Up);
      case 1:
        return RelWeightAlphaS(Var::Down);
    }
  }

  // Should never reach this point
  throw HZZException{"Illegal variation index."};
}


double GenWeight::RelWeightAlphaS(Var direction) const {
  if (not alphaSWeightsPresent_ or direction == Var::Nominal)
    return 1.;

  double weight;
  if (direction == Var::Up)
    weight = srcPdfWeights_[alphaSWeightsIndices_[0]];
  else
    weight = srcPdfWeights_[alphaSWeightsIndices_[1]];
  double const scaledWeight = 1. + (weight - 1.) * alphaSVarScaleFactor_;
  return scaledWeight;
}


double GenWeight::RelWeightMEScale(Var renorm, Var factor) const {
  if (not lheScaleWeightsPresent_)
    return 1.;
  if (srcScaleWeights_.GetSize() < 9) {
    HZZException exception;
    exception << "Cannot access ME scale variations (weights with indices 0 "
        << "to 8) because only " << srcScaleWeights_.GetSize()
        << " weights are available.";
    throw exception;
  }
  double const weight = srcScaleWeights_[meScaleIndices_.at({renorm, factor})];
  return weight / srcScaleWeights_[4];
}


double GenWeight::RelWeightPdf(Var direction) const {
  if (not pdfWeightsPresent_ or direction == Var::Nominal)
    return 1.;

  double sumW = 0., sumW2 = 0.;
  for (int i = pdfWeightsIndices_.first; i < pdfWeightsIndices_.second; ++i) {
    double const w = srcPdfWeights_[i];
    sumW += w;
    sumW2 += w * w;
  }
  int const n = pdfWeightsIndices_.second - pdfWeightsIndices_.first;
  double const meanWeight = sumW / n;
  double const sumDiff2 = sumW2 - n * std::pow(meanWeight, 2);

  double deltaWeight;
  if (pdfVarType_ == PdfVarType::Hessian)
    // Deviations from mean weight are summed up in quadrature
    deltaWeight = std::sqrt(sumDiff2);
  else
    // Standard deviation of the weights
    deltaWeight = std::sqrt(sumDiff2 / (n - 1));

  // When applying the variation, assume that the mean weight is 1. In practice,
  // it's very close to 1.
  if (direction == Var::Up)
    return 1. + deltaWeight;
  else
    return 1. - deltaWeight;
}


void GenWeight::InitializeLheScale(Dataset &dataset) {
  lheScaleWeightsPresent_ = false;
  if (auto const &weightInfo = dataset.Info().Parameters()["weights"];
      weightInfo) {
    if (auto const &lheScale = weightInfo["lhe_scale"]; lheScale)
      lheScaleWeightsPresent_ = lheScale.as<bool>();
  }

  // Drop WG samples of 2017/2018 as they have one fewer weight. See
  // https://github.com/cms-nanoAOD/cmssw/issues/520
  auto const WGNode = dataset.Info().Parameters()["wgamma_lnugamma"];
  if (WGNode and not WGNode.IsNull())  // Also drop the LO case for consistency
    lheScaleWeightsPresent_ = false;

  // Set up the mapping between the ME scale variations and indices in the
  // vector of weights. The convention for the weights is extracted from [1].
  // [1] https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc80X_doc.html#LHEScaleWeight
  meScaleIndices_[{Var::Down, Var::Down}] = 0;
  meScaleIndices_[{Var::Down, Var::Nominal}] = 1;
  meScaleIndices_[{Var::Down, Var::Up}] = 2;
  meScaleIndices_[{Var::Nominal, Var::Down}] = 3;
  meScaleIndices_[{Var::Nominal, Var::Nominal}] = 4;
  meScaleIndices_[{Var::Nominal, Var::Up}] = 5;
  meScaleIndices_[{Var::Up, Var::Down}] = 6;
  meScaleIndices_[{Var::Up, Var::Nominal}] = 7;
  meScaleIndices_[{Var::Up, Var::Up}] = 8;
}


void GenWeight::InitializePdf(Dataset &dataset) {
  std::string const pdfBranchTitle{
      dataset.Reader().GetTree()->GetBranch("LHEPdfWeight")->GetTitle()};
  if (pdfBranchTitle.empty()) {
    LOG_WARN << "Weights for PDF variations are not found.";
    pdfWeightsPresent_ = false;
    alphaSWeightsPresent_ = false;
    return;
  }
  pdfWeightsPresent_ = true;

  // Extract the range of LHAPDF IDs of stored PDF weights from the title of
  // the branch
  std::regex pattern{R"(LHA\s+IDs\s+(\d+)\s*-\s*(\d+)\b)"};
  std::smatch match;
  if (not std::regex_search(pdfBranchTitle, match, pattern)) {
    HZZException exception;
    exception << "Failed to parse title of branch with PDF weights \""
        << pdfBranchTitle << "\".";
    throw exception;
  }
  int const lhapdfIdFirst = std::stoi(match[1]);
  int const lhapdfIdLast = std::stoi(match[2]);

  int lhapdfId;  // ID of the PDF set as a whole
  bool nominalWeightPresent;
  if (lhapdfIdFirst % 100 == 1) {
    nominalWeightPresent = false;
    lhapdfId = lhapdfIdFirst - 1;
  } else {
    nominalWeightPresent = true;
    lhapdfId = lhapdfIdFirst;
  }

  auto const &[pdfVarType, alphaSVarIndices, alphaSVarSize]
      = LookUpPdfSet(lhapdfId);
  alphaSWeightsPresent_ = bool(alphaSVarIndices);

  pdfWeightsIndices_.first = (nominalWeightPresent) ? 1 : 0;
  pdfWeightsIndices_.second = lhapdfIdLast - lhapdfIdFirst + 1;
  if (alphaSWeightsPresent_)
    pdfWeightsIndices_.second -= 2;  // Last two weights are for alpha_s
  pdfVarType_ = pdfVarType;

  if (alphaSWeightsPresent_) {
    alphaSWeightsIndices_ = alphaSVarIndices.value();
    // PDF4LHC15 recommends a variation of 1.5e-3
    alphaSVarScaleFactor_ = 1.5e-3 / alphaSVarSize;
  }

  LOG_DEBUG << "Found PDF weights for LHAPDF IDs " << lhapdfIdFirst << " to "
      << lhapdfIdLast << ", which " << ((nominalWeightPresent) ? "" : "don't ")
      << "include nominal weight.";
  LOG_DEBUG << "Type of PDF variation: \""
      << ((pdfVarType_ == PdfVarType::Hessian) ? "Hessian" : "MC") << "\".";
  if (alphaSWeightsPresent_)
    LOG_DEBUG << "Variation in alpha_s is present. It will be scaled by "
        << "factor " << alphaSVarScaleFactor_ << ".";
  else
    LOG_DEBUG << "Variation in alpha_s is missing.";
}


std::tuple<GenWeight::PdfVarType, std::optional<std::array<int, 2>>, double>
GenWeight::LookUpPdfSet(int lhapdf) {
  YAML::Node pdfInfos = YAML::LoadFile(FileInPath::Resolve("pdf.yaml"));
  auto pdfInfo = pdfInfos[lhapdf];
  if (not pdfInfo) {
    HZZException exception;
    exception << "Found no metadata for PDF set with LHAPDF ID "
        << lhapdf << ".";
    throw exception;
  }

  for (auto const &field : {
      "pdf_variations_type", "alpha_s_indices", "alpha_s_variation"})
    if (auto node = pdfInfo[field]; not node) {
      HZZException exception;
      exception << "Entry with metadata for PDF set " << lhapdf
          << " doesn't contain mandatory field \"" << field << "\".";
      throw exception;
    }

  PdfVarType pdfVarType;
  auto const s = pdfInfo["pdf_variations_type"].as<std::string>();
  if (s == "Hessian") {
    pdfVarType = PdfVarType::Hessian;
  } else if (s == "MC") {
    pdfVarType = PdfVarType::MC;
  } else {
    HZZException exception;
    exception << "Value \"" << s << "\" for \"pdf_variations_type\" in "
        << "metadata for PDF set " << lhapdf << " is not supported.";
    throw exception;
  }

  std::optional<std::array<int, 2>> alphaSVarIndices;
  if (auto node = pdfInfo["alpha_s_indices"]; not node.IsNull()) {
    if (not node.IsSequence() or node.size() != 2) {
      HZZException exception;
      exception << "Value for \"alpha_s_indices\" in metadata for PDF set "
          << lhapdf << " is not an array of size 2.";
      throw exception;
    }
    alphaSVarIndices.emplace();
    (*alphaSVarIndices)[0] = node[0].as<int>();
    (*alphaSVarIndices)[1] = node[1].as<int>();
  }

  double alphaSVarSize = 0.;
  if (auto node = pdfInfo["alpha_s_variation"]; not node.IsNull())
    alphaSVarSize = node.as<double>();

  return {pdfVarType, alphaSVarIndices, alphaSVarSize};
}
