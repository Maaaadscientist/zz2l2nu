#include <Dataset.h>

#include <algorithm>
#include <fstream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <boost/algorithm/string.hpp>

#include <FileInPath.h>
#include <Logger.h>


namespace fs = std::filesystem;


DatasetInfo::DatasetInfo(fs::path const &path, Options const &options)
    : definitionFile_{path},
      crossSection_{std::numeric_limits<double>::quiet_NaN()},
      numEventsTotal_{0},
      meanWeight_{std::numeric_limits<double>::quiet_NaN()} {

  if (options.GetConfig()["dataset_stems"])
    stemsFile_ = FileInPath::Resolve(
        options.GetConfig()["dataset_stems"].as<std::string>());

  if (not fs::exists(path) or not fs::is_regular_file(path)) {
    std::ostringstream message;
    message << "Dataset definition file " << path << " does not exist "
        "or is not a regular file.";
    throw std::runtime_error(message.str());
  }

  if (boost::algorithm::ends_with(path.string(), ".txt"))
    ParseText(path);
  else
    ReadYaml(path);
}


YAML::Node const DatasetInfo::FindStem(std::string_view name) const {
  if (stemsFile_.empty()) {
    std::ostringstream message;
    message << "Location of file with dataset stems has not been set.";
    throw std::runtime_error(message.str());
  }

  if (not fs::exists(stemsFile_) or not fs::is_regular_file(stemsFile_)) {
    std::ostringstream message;
    message << "File with dataset stems " << stemsFile_ << " does not exist "
       "or is not a regular file.";
    throw std::runtime_error(message.str());
  }

  auto stems = YAML::LoadFile(stemsFile_);

  if (not stems.IsSequence() or stems.size() == 0) {
    std::ostringstream message;
    message << "File with dataset stems " << stemsFile_ << " does not contain "
        "a sequence of stems or the sequence is empty.";
    throw std::runtime_error(message.str());
  }

  for (auto const &stem : stems) {
    if (not stem["name"]) {
      std::ostringstream message;
      message << "An entry in file with dataset stems " << stemsFile_
          << " does not contain mandatory parameter \"name\".";
      throw std::runtime_error(message.str());
    }

    if (stem["name"].as<std::string>() == name)
      return stem;
  }

  // If this point is reached, the stem has not been found
  std::ostringstream message;
  message << "No stem for name \"" << name << "\" was found in file "
      << stemsFile_ << ".";
  throw std::runtime_error(message.str());
}


YAML::Node const DatasetInfo::GetNode(YAML::Node const root,
                                      std::string const &key) const {
  auto const node = root[key];

  if (not node or not node.IsScalar()) {
    std::ostringstream message;
    message << "Dataset definition file " << definitionFile_
        << " does not contain mandatory scalar field \"" << key << "\".";
    throw std::runtime_error(message.str());
  }

  return node;
}


void DatasetInfo::ParseText(fs::path const &path) {
  LOG_WARN << "Reading an old-style catalogue file. There is only a limited "
      "support for this format.";

  // Regular expression that matches a blank line
  std::regex blankRegex("^\\s*$");

  // Regular expression that matches a line containing a file path. The first
  // group receives the path.
  std::regex filePathRegex("^\\s*(\\S+).*$");

  // Regular expression that matches a line with configuration. The first group
  // receives the name of the parameter, the second gets its value.
  std::regex configRegex("^\\s*[#\\*]\\s*(.+)\\s*:\\s*(.*)\\s*$");

  std::smatch match;


  std::ifstream datasetFile(path);
  std::string line;

  if (not datasetFile) {
    std::ostringstream message;
    message << "Cannot open dataset definition file " << path << ".";
    throw std::runtime_error(message.str());
  }

  std::map<std::string, std::string> parameters;
  bool readingHeader = true;

  while (std::getline(datasetFile, line)) {
    if (std::regex_match(line, blankRegex))
      continue;

    if (readingHeader) {
      if (std::regex_match(line, match, configRegex)) {
        auto const &name = match[1];
        auto const &value = match[2];
        auto const res = parameters.find(name);

        if (res != parameters.end())
          LOG_WARN << "Parameter \"" << name << "\" specified multiple times "
              "in dataset definition file " << path
              << ". Overwriting old value \"" << res->second << "\" with \""
              << value << "\".";

        parameters[name] = value;
      } else
        readingHeader = false;
    }

    if (not readingHeader) {
      if (std::regex_match(line, match, filePathRegex))
        files_.emplace_back(match[1]);
      else
        LOG_WARN << "In dataset definition file " << path
            << ", failed to parse line\n" << line << "\nSkipping it.";
    }
  }

  datasetFile.close();


  LOG_TRACE << "Parameters read from dataset definition file " << path << ":";

  for (auto const &p : parameters)
    LOG_TRACE << "  \"" << p.first << "\": \"" << p.second << "\"";

  LOG_TRACE << "Total number of paths of input files found: " << files_.size();

  if (files_.empty()) {
    std::ostringstream message;
    message << "No paths to input files read from dataset definition file "
        << path << ".";
    throw std::runtime_error(message.str());
  }


  // Save important parameters
  auto res = parameters.find("data type");

  if (res == parameters.end()) {
    LOG_WARN << "Dataset definition file " << path
        << " does not contain parameter \"" << res->first
        << "\". Going to assume this is simulation.";
    isSimulation_ = true;
  }

  std::string dataType{res->second};
  boost::to_lower(dataType);

  if (dataType == "mc")
    isSimulation_ = true;
  else if (dataType == "data")
    isSimulation_ = false;
  else {
    std::ostringstream message;
    message << "Illegal value \"" << res->second << "\" for parameter \""
        << res->first << "\" found in dataset definition file " << path << ".";
    throw std::runtime_error(message.str());
  }


  if (isSimulation_) {
    res = parameters.find("sample xsec");

    if (res == parameters.end()) {
      std::ostringstream message;
      message << "Dataset definition file " << path <<
          " does not contain mandatory parameter \"" << res->first << "\".";
      throw std::runtime_error(message.str());
    }

    crossSection_ = std::stod(res->second);


    res = parameters.find("primary events");

    if (res == parameters.end()) {
      std::ostringstream message;
      message << "Dataset definition file " << path <<
          " does not contain mandatory parameter \"" << res->first << "\".";
      throw std::runtime_error(message.str());
    }

    numEventsTotal_ = std::stoll(res->second);


    // Mean weight is never stored in catalogues. Assume 1.
    meanWeight_ = 1.;
  }


  // Save all parameters in the YAML representation. Translate them to the
  // format used in YAML dataset definition files if needed.
  for (auto const &[name, value] : parameters) {
    if (name == "data type")
      parameters_["is_sim"] = isSimulation_;
    else if (name == "sample xsec")
      parameters_["cross_section"] = crossSection_;
    else if (name == "primary events")
      parameters_["num_events"] = numEventsTotal_;
    else
      parameters_[name] = value;
  }
}


void DatasetInfo::ReadYaml(fs::path const &path) {
  YAML::Node info = YAML::LoadFile(path);

  if (info["stem"])
    // This is not a full dataset definition file. Extend it with the stem.
    SpliceYaml(info);

  YAML::Node filesNode = info["files"];

  if (not filesNode or not filesNode.IsSequence()) {
    std::ostringstream message;
    message << "Dataset definition file " << path
        << " does not contain mandatory sequence \"files\".";
    throw std::runtime_error(message.str());
  }

  for (auto const &element : filesNode)
    files_.emplace_back(element.as<std::string>());

  LOG_TRACE << "Total number of paths of input files found in dataset "
      "definition file " << path << ": " << files_.size();

  if (files_.empty()) {
    std::ostringstream message;
    message << "No paths to input files read from dataset definition file "
        << path << ".";
    throw std::runtime_error(message.str());
  }


  // Save important parameters
  isSimulation_ = GetNode(info, "is_sim").as<bool>();

  if (isSimulation_) {
    crossSection_ = GetNode(info, "cross_section").as<double>();
    numEventsTotal_ = GetNode(info, "num_events").as<int64_t>();

    YAML::Node const weightNode = info["mean_weight"];

    if (weightNode)
      meanWeight_ = weightNode.as<double>();
    else {
      LOG_DEBUG << "Mean weight is not found in dataset definition file. "
          "Setting it to 1.";
      meanWeight_ = 1.;
    }
  }


  // Save as parameters the full YAML configuration except for the list of
  // input files
  parameters_ = info;
  parameters_.remove("files");
}


void DatasetInfo::SpliceYaml(YAML::Node info) const {
  std::string const name{info["stem"].as<std::string>()};
  LOG_DEBUG << "Splicing dataset definition fragment " << definitionFile_
      << " with stem \"" << name << "\".";

  auto const stem = FindStem(name);
  info.remove("stem");

  for (auto it = stem.begin(); it != stem.end(); ++it)
    info[it->first] = it->second;
}


Dataset::Dataset(DatasetInfo info, int skipFiles, int maxFiles)
    : info_{std::move(info)}, chain_{"tupel/EventTree"} {

  if (skipFiles < 0) {
    std::ostringstream message;
    message << "Illegal value " << skipFiles
        << " given to argument \"skipFiles\".";
    throw std::runtime_error(message.str());
  }

  if (maxFiles < -1) {
    std::ostringstream message;
    message << "Illegal value " << maxFiles
        << " given to argument \"maxFiles\".";
    throw std::runtime_error(message.str());
  }

  int const numFilesTotal = info_.Files().size();
  int end;

  if (maxFiles == -1)
    end = numFilesTotal;
  else
    end = std::min(skipFiles + maxFiles, numFilesTotal);

  for (int i = skipFiles; i < end; ++i) {
    auto const &path = info_.Files()[i];
    selectedFiles_.emplace_back(path);
    chain_.AddFile(path.c_str());
    LOG_DEBUG << "File \"" << path << "\" added to the list of input files.";
  }

  if (selectedFiles_.empty())
    LOG_WARN << "An empty dataset was constructed.";

  reader_.SetTree(&chain_);
}

