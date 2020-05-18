#include <LeptonWeight.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <TH2.h>
#include <TFile.h>
#include <yaml-cpp/yaml.h>

#include <FileInPath.h>
#include <Logger.h>


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
    throw std::runtime_error(
        "Configuration for a lepton scale factor component does not contain "
        "mandatory parameter \"path\".");
  auto const path = config["path"].as<std::string>();
  histogram_ = ReadHistogram(path);

  auto const &schemaNode = config["schema"];
  if (not schemaNode) {
    std::ostringstream message;
    message << "Configuration for lepton scale factor component with path \""
        << path << "\" does not contain mandatory paramter \"schema\".";
    throw std::runtime_error(message.str());
  }
  if (not schemaNode.IsSequence() or schemaNode.size() != 2) {
    std::ostringstream message;
    message << "Illegal schema in configuration for lepton scale factor "
        << "component with path \"" << path << "\". It must be a sequence "
        << "containing exactly two elements.";
    throw std::runtime_error(message.str());
  }
  auto const schema = schemaNode.as<std::vector<std::string>>();
  int numPtTags = 0;
  for (auto const &tag : schema) {
    if (tag != "pt" and tag != "eta" and tag != "abs_eta") {
      std::ostringstream message;
      message << "Schema in configuration for lepton scale factor component "
          << "with path \"" << path << "\" contains unknown tag \"" << tag
          << "\".";
      throw std::runtime_error(message.str());
    }
    if (tag == "pt")
      ++numPtTags;
  }
  if (numPtTags != 1) {
    std::ostringstream message;
    message << "Illegal schema in configuratino for lepton scale factor "
        << "component with path \"" << path << "\". It must contain exactly "
        << "one tag for pt and exactly one for pseudorapidity.";
    throw std::runtime_error(message.str());
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

  std::unique_ptr<TH2> hist{dynamic_cast<TH2 *>(inputFile.Get(name.c_str()))};
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


LeptonWeight::LeptonWeight(Dataset &, Options const &options,
                           ElectronBuilder const *electronBuilder,
                           MuonBuilder const *muonBuilder)
    : electronBuilder_{electronBuilder}, muonBuilder_{muonBuilder} {

  auto const &muonComponents = Options::NodeAs<YAML::Node>(
      options.GetConfig(), {"lepton_efficiency", "muon"});
  for (auto const &config : muonComponents)
    muonScaleFactors_.emplace_back(config);

  auto const &electronComponents = Options::NodeAs<YAML::Node>(
      options.GetConfig(), {"lepton_efficiency", "electron"});
  for (auto const &config : electronComponents)
    electronScaleFactors_.emplace_back(config);

  LOG_WARN << "Trigger scale factors are missing";
}


// Define the descructor at a point where PtEtaHistogram is a complete class
LeptonWeight::~LeptonWeight() {}


double LeptonWeight::NominalWeight() const {
  double sf = 1.;
  for (auto &electron : electronBuilder_->GetTight()) {
    sf *= ElectronSF(electron);
  }
  for (auto &muon : muonBuilder_->GetTight()) {
    sf *= MuonSF(muon);
  }
  return sf;
}


double LeptonWeight::ElectronSF(Electron const &electron) const {
  double pt = electron.p4.Pt();
  double eta = electron.etaSc;
  double sf = 1.;
  for (auto const &component : electronScaleFactors_)
    sf *= component(pt, eta);
  return sf;
}


double LeptonWeight::MuonSF(Muon const &muon) const {
  double pt = muon.uncorrP4.Pt();
  double eta = muon.uncorrP4.Eta();
  double sf = 1.;
  for (auto const &component : muonScaleFactors_)
    sf *= component(pt, eta);
  return sf;
}

