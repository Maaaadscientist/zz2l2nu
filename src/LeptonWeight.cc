#include <LeptonWeight.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <string>

#include <TH2.h>
#include <TFile.h>
#include <yaml-cpp/yaml.h>

#include <HZZException.h>
#include <FileInPath.h>
#include <Logger.h>
#include <Utils.h>


/**
 * \brief Auxiliary class that abstracts true representation of a 2D ROOT
 * histogram parameterized with pt and pseudorapidity
 *
 * The two dimensions of the underlying histogram can be ordered in any way.
 * It can also be parameterized with signed or absolute pseudorapidity. The
 * overflow bins in pt can be filled or be empty, which implies that values from
 * the last pt bins should be used. This class provides a consistent interface
 * to read values from the underlying histogram regardless of these aspects.
 */
class PtEtaHistogram {
 public:
  /**
   * \brief Construct from configuration
   *
   * The configuration must be a mapping with the folowing structure:
   * \code{.yaml}
   * path: {file_path}:{in_file_path}
   * schema: [{tag1}, {tag2}]
   * clip_pt: {bool}
   * \endcode
   * Here \c {file_path} is path to a ROOT file (which is resolved with
   * FileInPath) and \c {in_file_path} is the location of a 2D histogram within
   * that file. \c {tag1} and \c {tag2} take two of the following supported
   * values: "pt", "eta", "abs_eta"; they describe the dimensions of the
   * histogram. Optional key \c {clip_pt} specifies whether pt should be clipped
   * (which is the case if the overflow bins in pt are not set correctly).
   */
  PtEtaHistogram(YAML::Node const &config);

  /// Retrieves the value from the histogram for the given pt and eta
  double operator()(double pt, double eta) const;

 private:
  /**
   * \brief Reads a histogram with given path and name
   *
   * Checks for and reports errors. The returned histogram is owned by the
   * caller.
   */
  static std::unique_ptr<TH2> ReadHistogram(std::string const &pathsWithNames);

  /// Underlying histogram
  std::unique_ptr<TH2> histogram_;

  /// Ordering of dimensions in the underlying 2D histogram
  bool orderPtEta_;

  /// Indicates that absolute (as opposed to signed) pseudorapidity is used
  bool useAbsEta_;

  /// Indicates whether pt should be clipped to a maximal value
  bool clipPt_;

  /**
   * \brief Maximal value of pt for clipping
   *
   * Only used when clipPt_ is true.
   */
  double maxPt_;
};


PtEtaHistogram::PtEtaHistogram(YAML::Node const &config) {
  if (not config["path"])
    throw HZZException{
        "Configuration for a lepton scale factor component does not contain "
        "mandatory parameter \"path\"."};
  auto const path = config["path"].as<std::string>();
  histogram_ = ReadHistogram(path);

  auto const &schemaNode = config["schema"];
  if (not schemaNode) {
    HZZException exception;
    exception << "Configuration for lepton scale factor component with path \""
        << path << "\" does not contain mandatory paramter \"schema\".";
    throw exception;
  }
  if (not schemaNode.IsSequence() or schemaNode.size() != 2) {
    HZZException exception;
    exception << "Illegal schema in configuration for lepton scale factor "
        << "component with path \"" << path << "\". It must be a sequence "
        << "containing exactly two elements.";
    throw exception;
  }
  auto const schema = schemaNode.as<std::vector<std::string>>();
  int numPtTags = 0;
  for (auto const &tag : schema) {
    if (tag != "pt" and tag != "eta" and tag != "abs_eta") {
      HZZException exception;
      exception << "Schema in configuration for lepton scale factor component "
          << "with path \"" << path << "\" contains unknown tag \"" << tag
          << "\".";
      throw exception;
    }
    if (tag == "pt")
      ++numPtTags;
  }
  if (numPtTags != 1) {
    HZZException exception;
    exception << "Illegal schema in configuratino for lepton scale factor "
        << "component with path \"" << path << "\". It must contain exactly "
        << "one tag for pt and exactly one for pseudorapidity.";
    throw exception;
  }

  int etaTagIndex;
  if (schema[0] == "pt") {
    orderPtEta_ = true;
    etaTagIndex = 1;
  } else {
    orderPtEta_ = false;
    etaTagIndex = 0;
  }
  useAbsEta_ = (schema[etaTagIndex] == "abs_eta");

  auto const &clipNode = config["clip_pt"];
  clipPt_ = (clipNode and clipNode.as<bool>());
  if (clipPt_) {
    if (orderPtEta_)
      maxPt_ = histogram_->GetXaxis()->GetBinCenter(histogram_->GetNbinsX());
    else
      maxPt_ = histogram_->GetYaxis()->GetBinCenter(histogram_->GetNbinsY());
  }
}


double PtEtaHistogram::operator()(double pt, double eta) const {
  if (clipPt_ and pt > maxPt_)
    pt = maxPt_;
  if (useAbsEta_)
    eta = std::abs(eta);

  double x, y;
  if (orderPtEta_) {
    x = pt;
    y = eta;
  } else {
    x = eta;
    y = pt;
  }
  auto const bin = histogram_->FindFixBin(x, y);
  return histogram_->GetBinContent(bin);
}


std::unique_ptr<TH2> PtEtaHistogram::ReadHistogram(
    std::string const &pathWithName) {
  auto const pos = pathWithName.find_last_of(':');
  if (pos == std::string::npos)
    throw HZZException{"Histogram path does not contain ':'."};

  std::filesystem::path path = FileInPath::Resolve(pathWithName.substr(0, pos));
  std::string name = pathWithName.substr(pos + 1);
  return utils::ReadHistogram<TH2>(path, name);
}


LeptonWeight::LeptonWeight(Dataset &dataset, Options const &options,
                           ElectronBuilder const *electronBuilder,
                           MuonBuilder const *muonBuilder)
    : cache_{dataset.Reader()}, electronBuilder_{electronBuilder}, muonBuilder_{muonBuilder} {
  // The default weight index is chosen based on the requested systematic
  // variation
  auto const systLabel = options.GetAs<std::string>("syst");
  if (systLabel == "muonEff_syst_up")
    defaultWeightIndex_ = 1;
  else if (systLabel == "muonEff_syst_down")
    defaultWeightIndex_ = 2;
  else if (systLabel == "muonEff_stat_up")
    defaultWeightIndex_ = 3;
  else if (systLabel == "muonEff_stat_down")
    defaultWeightIndex_ = 4;
  else if (systLabel == "electronEff_syst_up")
    defaultWeightIndex_ = 5;
  else if (systLabel == "electronEff_syst_down")
    defaultWeightIndex_ = 6;
  else if (systLabel == "electronEff_stat_up")
    defaultWeightIndex_ = 7;
  else if (systLabel == "electronEff_stat_down")
    defaultWeightIndex_ = 8;
  else  
    defaultWeightIndex_ = 0;
  LOG_DEBUG << "Index of default leptonSF weight: " << defaultWeightIndex_;
  for (int i = 0; i < NumVariations() / 2 + 1; i++){		
    std::string systName;
    if (i == 0) systName = "nominal";
    else {
      systName = VariationName(i-1);
      auto const pos = systName.find_first_of("_"); //format: (electron or muon)Eff_(systType)_(up or down)
      systName = systName.substr(pos + 1);
    }
    auto const muonComponents = Options::NodeAs<YAML::Node>(
        options.GetConfig(), {"lepton_efficiency", "muon", systName});
    for (auto const &config : muonComponents)
      muonScaleFactors_[i].emplace_back(config);

    auto const electronComponents = Options::NodeAs<YAML::Node>(
        options.GetConfig(), {"lepton_efficiency", "electron", systName});
    for (auto const &config : electronComponents)
      electronScaleFactors_[i].emplace_back(config);
  }
  LOG_WARN << "Trigger scale factors are missing";
}

// Define the descructor at a point where PtEtaHistogram is a complete class
LeptonWeight::~LeptonWeight() {}

std::string_view LeptonWeight::VariationName(int variation) const {
  switch (variation) {
    case 0:
      return "muonEff_syst_up";
    case 1:
      return "muonEff_syst_down";
    case 2:
      return "muonEff_stat_up";
    case 3:
      return "muonEff_stat_down";
    case 4:
      return "electronEff_syst_up";
    case 5:
      return "electronEff_syst_down";
    case 6:
      return "electronEff_stat_up";
    case 7:
      return "electronEff_stat_down";
    default:
      return "";
  }
}

void LeptonWeight::Update() const {
  for(int i = 0;i < NumVariations() + 1; i++){
  	double sf = 1.;
    for (auto &electron : electronBuilder_->GetTight()) {
      sf *= ElectronSF(electron, i > NumVariations() ? i - NumVariations() / 2 : 0);
    }
    for (auto &muon : muonBuilder_->GetTight()) {
      sf *= MuonSF(muon, i <= NumVariations() ? i : 0);
    }
    weights_[i] = sf;
  }
}

double LeptonWeight::ElectronSF(Electron const &electron, int syst) const {
  double pt = electron.p4.Pt();
  double eta = electron.etaSc;
  double sf = 1.;
  for (auto const &component : electronScaleFactors_[syst])
    sf *= component(pt, eta);
  return sf;
}

double LeptonWeight::MuonSF(Muon const &muon, int syst) const {
  double pt = muon.uncorrP4.Pt();
  double eta = muon.uncorrP4.Eta();
  double sf = 1.;
  for (auto const &component : muonScaleFactors_[syst])
    sf *= component(pt, eta);
  return sf;
}
