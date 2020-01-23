#include <LeptonWeight.h>

#include <sstream>
#include <stdexcept>

#include <TFile.h>

#include <FileInPath.h>
#include <Logger.h>

LeptonWeight::LeptonWeight(Dataset &, Options const &options,
                           ElectronBuilder const *electronBuilder,
                           MuonBuilder const *muonBuilder)
    : electronBuilder_{electronBuilder}, muonBuilder_{muonBuilder} {
  for (auto const &muonPath : 
    Options::NodeAs<std::vector<std::string>>(
      options.GetConfig(), {"lepton_efficiency", "muon"})) {
    muonTable_.emplace_back(ReadHistogram(muonPath));
  }
  for (auto const &electronPath : 
    Options::NodeAs<std::vector<std::string>>(
      options.GetConfig(), {"lepton_efficiency", "electron"})) {
    electronTable_.emplace_back(ReadHistogram(electronPath));
  }
  LOG_WARN << "Trigger scale factors are missing";
}

double LeptonWeight::NominalWeight() const {
  double eff = 1.;
  for (auto &electron : electronBuilder_->GetTight()) {
    eff *= ElectronSF(electron);
  }
  for (auto &muon : muonBuilder_->GetTight()) {
    eff *= MuonSF(muon);
  }
  return eff;
}

double LeptonWeight::ElectronSF(Electron const &electron) const {
  double eff = 1.;
  
  double pt = electron.p4.Pt();
  double eta = electron.etaSc;
    
  for (auto &electronTable : electronTable_) {

    int eleEtaBin = electronTable->GetXaxis()->FindFixBin(eta);
    int elePtBin = electronTable->GetYaxis()->FindFixBin(pt);
    if (pt > electronTable->GetYaxis()->GetXmax()) 
      elePtBin = electronTable->GetNbinsY();

    eff *= electronTable->GetBinContent(eleEtaBin, elePtBin);
  }

  return eff;
}

double LeptonWeight::MuonSF(Muon const &muon) const {
  double eff = 1.;
  
  double pt = muon.uncorrP4.Pt();
  double eta = muon.uncorrP4.Eta();

  for (auto &muonTable : muonTable_) {
    int muEtaBin = muonTable->GetXaxis()->FindFixBin(eta);
    int muPtBin = muonTable->GetYaxis()->FindFixBin(pt);
    if (pt > muonTable->GetYaxis()->GetXmax())
      muPtBin = muonTable->GetNbinsY();
    
    eff *= muonTable->GetBinContent(muEtaBin, muPtBin);
  }     
  
  return eff;
}

std::unique_ptr<TH2> LeptonWeight::ReadHistogram(
    std::string const &pathWithName) {
  auto const pos = pathWithName.find_last_of(':');
  
  if (pos == std::string::npos) {
    std::ostringstream message;
    message << "Histogram path does not contain ':'.";
    throw std::invalid_argument(message.str());
  }
  
  std::filesystem::path path = FileInPath::Resolve(pathWithName.substr(0, pos));
  std::string name = pathWithName.substr(pos + 1);
  
  TFile inputFile{path.c_str()};

  if (inputFile.IsZombie()) {
    std::ostringstream message;
    message << "Could not open file \"" << path << "\".";
    throw std::runtime_error(message.str());
  }

  std::unique_ptr<TH2> hist;
  hist.reset(dynamic_cast<TH2 *>(inputFile.Get(name.c_str())));
  
  if (not hist) {
    std::ostringstream message;
    message << "File " << path << 
      " does not contain required histogram \"" <<
      name << "\".";
    throw std::runtime_error(message.str());
  }

  hist->SetDirectory(nullptr);
  inputFile.Close();

  return hist;
}
