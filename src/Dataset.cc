#include <Dataset.h>

#include <algorithm>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <boost/algorithm/string.hpp>

#include <Logger.h>


namespace fs = std::filesystem;


DatasetInfo::DatasetInfo(fs::path const &path)
    : definitionFile_{path} {

  if (not fs::exists(path) or not fs::is_regular_file(path)) {
    std::ostringstream message;
    message << "Dataset definition file " << path << " does not exist "
        "or is not a regular file.";
    throw std::runtime_error(message.str());
  }

  if (boost::algorithm::ends_with(path.string(), ".txt"))
    ParseText(path);
  else
    ParseYaml(path);
}


void DatasetInfo::ParseText(fs::path const &path) {

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


  // Save all parameters in the YAML representation. Translate them to the
  // format used in YAML dataset definition files if needed.
  for (auto const &[name, value] : parameters) {
    if (name == "data type")
      parameters_["is_sim"] = isSimulation_;
    else
      parameters_[name] = value;
  }
}


void DatasetInfo::ParseYaml(fs::path const &path) {
  YAML::Node info = YAML::LoadFile(path);
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
  YAML::Node node = info["is_sim"];

  if (not node or not node.IsScalar()) {
    std::ostringstream message;
    message << "Dataset definition file " << path
        << " does not contain mandatory scalar field \"is_sim\".";
    throw std::runtime_error(message.str());
  }

  isSimulation_ = node.as<bool>();


  // Save as parameters the full YAML configuration except for the list of
  // input files
  parameters_ = info;
  parameters_.remove("files");
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

