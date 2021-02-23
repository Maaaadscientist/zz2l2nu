#include <TriggerWeight.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <string>

#include <TFile.h>
#include <TList.h>
#include <yaml-cpp/yaml.h>

#include <HZZException.h>
#include <FileInPath.h>
#include <Logger.h>
#include <Utils.h>


TriggerWeight::TriggerWeight(Dataset &dataset, Options const &options,
    ElectronBuilder const *electronBuilder,
    MuonBuilder const *muonBuilder, int efficiencyType)
    : cache_{dataset.Reader()}, electronBuilder_{electronBuilder}, muonBuilder_{muonBuilder} {
  efficiencyType_ = efficiencyType;
  // The default weight index is chosen based on the requested systematic
  // variation
  auto const systLabel = options.GetAs<std::string>("syst");
  if (systLabel == "triggerEE_up")
    defaultWeightIndex_ = 1;
  else if (systLabel == "triggerEE_down")
    defaultWeightIndex_ = 2;
  else if (systLabel == "triggerMuMu_up")
    defaultWeightIndex_ = 3;
  else if (systLabel == "triggerMuMu_down")
    defaultWeightIndex_ = 4;
  else if (systLabel == "triggerEMu_up")
    defaultWeightIndex_ = 5;
  else if (systLabel == "triggerEMu_down")
    defaultWeightIndex_ = 6;
  else  
    defaultWeightIndex_ = 0;
  LOG_DEBUG << "Index of default trigger  weight: " << defaultWeightIndex_;
  auto const mumuPath = Options::NodeAs<std::string>(
      options.GetConfig(), {"trigger_efficiency", "mumu", "path"});
  auto const eePath = Options::NodeAs<std::string>(
      options.GetConfig(), {"trigger_efficiency", "ee", "path"});
  auto const emuPath = Options::NodeAs<std::string>(
      options.GetConfig(), {"trigger_efficiency", "emu", "path"});

  ReadHistogram(mumuPath, mumuScaleFactors_);
  ReadHistogram(eePath, eeScaleFactors_);
  ReadHistogram(emuPath, emuScaleFactors_);
}

TriggerWeight::~TriggerWeight() {}

std::string_view TriggerWeight::VariationName(int variation) const {
  switch (variation) {
    case 0:
      return "triggerEE_up";
    case 1:
      return "triggerEE_down";
    case 2:
      return "triggerMuMu_up";
    case 3:
      return "triggerMuMu_down";
    case 4:
      return "triggerEMu_up";
    case 5:
      return "triggerEMu_down";
    default:
      return "";
  }
}

void TriggerWeight::ReadHistogram(
    std::string const &pathWithName,
    std::map<std::string, std::unique_ptr<TH2D>> &th2Map) {

  std::filesystem::path path = FileInPath::Resolve(pathWithName);
  TFile *tmpFile = new TFile(path.c_str());
  auto const listOfKeys = tmpFile->GetListOfKeys();
  for (int i = 0; i < listOfKeys->GetSize(); i++){
    th2Map.emplace(listOfKeys->At(i)->GetName(), 
        utils::ReadHistogram<TH2D>(path, listOfKeys->At(i)->GetName()));
  }
}

void TriggerWeight::Update() const {
  auto const &tightElectrons = electronBuilder_->GetTight();
  auto const &tightMuons = muonBuilder_->GetTight();
  int leptonCat;
  Lepton const *l1, *l2;
  if (tightElectrons.size() == 2) {
    leptonCat = kEE;
    l1 = &tightElectrons[0];
    l2 = &tightElectrons[1];
  } else if (tightMuons.size() == 2) {
    leptonCat = kMuMu;
    l1 = &tightMuons[0];
    l2 = &tightMuons[1];
  } else if (tightElectrons.size() == 1 and tightMuons.size() == 1) {
    leptonCat = kEMu;
    l1 = &tightElectrons[0];
    l2 = &tightMuons[0];
  } else {
    for (int aSyst = 0 ; aSyst < NumVariations() + 1; aSyst++)
      weights_[aSyst] = 1.;
    l1 = l2 = nullptr;
    leptonCat = 3;
  }
  if (leptonCat == kEE or leptonCat == kMuMu or leptonCat == kEMu)
    for (int aSyst = 0 ; aSyst < NumVariations() + 1; aSyst++)
      weights_[aSyst] = GetEfficiency(l1, l2, leptonCat, aSyst);
}

double TriggerWeight::GetEfficiency(Lepton const *l1, Lepton const *l2, 
    int leptonCat, int syst) const {
  std::string cat = "";
  double eta1, eta2, pt1, pt2;
  eta1 = std::abs(l1->p4.Eta());
  eta2 = std::abs(l2->p4.Eta());
  pt1 = (l1->p4.Pt() > 200 ? 150. : l1->p4.Pt());
  pt2 = (l2->p4.Pt() > 200 ? 150. : l2->p4.Pt());
  if ((leptonCat == 0 || leptonCat == 1) && pt1 < pt2) {
    double tmp;
    tmp = pt2;
    pt2 = pt1;
    pt1 = tmp;
  }
  std::string systTypeEE = "nominal";
  std::string systTypeMuMu = "nominal";
  std::string systTypeEMu = "nominal";
  switch(syst){
    case 0:
      ;
      break;
    case 1:
      systTypeEE = "up";
      break;
    case 2:
      systTypeEE = "down";
      break;
    case 3:
      systTypeMuMu = "up";
      break;
    case 4:
      systTypeMuMu = "down";
      break;
    case 5:
      systTypeEMu = "up";
      break;
    case 6:
      systTypeEMu = "down";
      break;
  }
  std::string type;
  switch(efficiencyType_){
    case scaleFactor:
      type = "sf";
      break;
    case dataEfficiency:
      type = "effData";
      break;
    case mcEfficiency:
      type = "effMC";
      break;
    default:
      throw HZZException("Unknown type of trigger efficiencies.");
  }
  double eff;
  switch(leptonCat){
    case kEE: {
      if (eta1 < 1.479 and eta2 < 1.479) cat = "BB";
      else if (eta1 >= 1.479 and eta2 < 1.479) cat = "EB";
      else if (eta1 < 1.479 and eta2 >= 1.479) cat = "BE";
      else if (eta1 >= 1.479 and eta2 >= 1.479) cat = "EE";
      auto searchEE = eeScaleFactors_.find(type+"_"+cat+"_"+systTypeEE);
      if (searchEE != eeScaleFactors_.end()){
        auto const bin = searchEE->second->FindFixBin(pt1, pt2);
        eff = searchEE->second->GetBinContent(bin);
      }
      else
        throw HZZException("Unknown efficiency table\n");
      break;
  	}
    case kMuMu: {
      if (eta1 < 1.2 and eta2 < 1.2) cat = "BB";
      else if (eta1 >= 1.2 and eta2 < 1.2) cat = "EB";
      else if (eta1 < 1.2 and eta2 >= 1.2) cat = "BE";
      else if (eta1 >= 1.2 and eta2 >= 1.2) cat = "EE";
      auto searchMuMu = mumuScaleFactors_.find(type+"_"+cat+"_"+systTypeMuMu);
      if (searchMuMu != mumuScaleFactors_.end()) {
        auto const bin = searchMuMu->second->FindFixBin(pt1, pt2);
        eff =  searchMuMu->second->GetBinContent(bin);
      }
      else 
        throw HZZException("Unknown efficiency table\n");
      break;
    }
    case kEMu: {
      if (eta1 < 1.479 and eta2 < 1.2) cat = "BB";
      else if (eta1 >= 1.479 and eta2 < 1.2) cat = "EB";
      else if (eta1 < 1.479 and eta2 >= 1.2) cat = "BE";
      else if (eta1 >= 1.479 and eta2 >= 1.2) cat = "EE";
      auto searchEMu = emuScaleFactors_.find(type+"_"+cat+"_"+systTypeEMu);
      if (searchEMu != emuScaleFactors_.end()){
        auto const bin = searchEMu->second->FindFixBin(pt1, pt2);
        eff =  searchEMu->second->GetBinContent(bin);
      }
      else 
        throw HZZException("Unknown efficiency table\n");
      break;
    }
    default:
      throw HZZException("Unknown lepton category\n");
      return 0.;
  }
  if (eff == 0.)
    LOG_WARN <<"Zero trigger efficiency.\n" ;
  return eff;
}

double TriggerWeight::TransferFactor(
    Lepton const *l1, Lepton const *l2, int syst) const {
  if (efficiencyType_ == 0) 
    throw HZZException("Cannot use scale factors for transfer factor\n");
  double tf = std::sqrt(
      GetEfficiency(l1, l2, kEE, syst) * GetEfficiency(l1, l2, kMuMu, syst)) /
      GetEfficiency(l1, l2, kEMu, syst);
  if (std::abs(tf) < 100. && std::abs(tf) != 0.)
    return std::abs(tf);
  else {
    LOG_WARN << "Trigger transfer factor has a weight of:" << tf << ".\n";
    return 0.;
  }
}
