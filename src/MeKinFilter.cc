#include <MeKinFilter.h>

#include <HZZException.h>
#include <Logger.h>


MeKinFilter::MeKinFilter(Dataset &dataset) {
  auto const settingsNode = dataset.Info().Parameters()["me_kin_filter"];
  enabled_ = (settingsNode and not settingsNode.IsNull()
              and dataset.Info().IsSimulation());

  if (enabled_) {
    for (auto const &key: {"variable", "range"}) {
      auto const node = settingsNode[key];
      if (not node) {
        HZZException exception;
        exception << "Mandatory key \"" << key << "\" is missing in section "
            "\"me_kin_filter\" for dataset \"" << dataset.Info().Name()
            << "\".";
        throw exception;
      }
    }

    auto const variable = settingsNode["variable"].as<std::string>();
    if (variable != "Ht") {
      HZZException exception;
      exception << "Variable \"" << variable << "\" is not supported.";
      throw exception;
    }

    auto const rangeNode = settingsNode["range"];
    if (not rangeNode.IsSequence() or rangeNode.size() != 2) {
      HZZException exception;
      exception << "Range for variable \"" << variable << "\" is not a "
          "sequence or contains a wrong number of elements.";
      throw exception;
    }

    minValue_ = rangeNode[0].as<double>();
    maxValue_ = rangeNode[1].as<double>();

    if (minValue_ >= maxValue_) {
      HZZException exception;
      exception << "Wrong ordering in range (" << minValue_ << ", " << maxValue_
          << ").";
      throw exception;
    }

    LOG_DEBUG << "Will apply selection " << minValue_ << " < " << variable
        << " < " << maxValue_ << " in the ME final state.";


    auto &reader = dataset.Reader();
    srcNumPart_.reset(new TTreeReaderValue<UInt_t>(reader, "nGenPart"));
    srcPdgId_.reset(new TTreeReaderArray<Int_t>(reader, "GenPart_pdgId"));
    srcMotherIndex_.reset(new TTreeReaderArray<Int_t>(
        reader, "GenPart_genPartIdxMother"));
    srcPartPt_.reset(new TTreeReaderArray<Float_t>(reader, "GenPart_pt"));
  }
}


bool MeKinFilter::operator()() const {
  if (not enabled_)
    return true;

  double const value = ComputeHt();
  return (value >= minValue_ and value < maxValue_);
}


double MeKinFilter::ComputeHt() const {
  double Ht = 0.;

  for (unsigned i = 0; i < *srcNumPart_->Get(); ++i) {
    if (srcMotherIndex_->At(i) != 0)
      continue;

    int const absPdgId = std::abs(srcPdgId_->At(i));
    if (absPdgId > 5 and absPdgId != 21)
      continue;

    Ht += srcPartPt_->At(i);
  }

  return Ht;
}
