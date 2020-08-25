#include <Dataset.h>

#include <algorithm>
#include <limits>
#include <map>
#include <utility>

#include <FileInPath.h>
#include <HZZException.h>
#include <Logger.h>


namespace fs = std::filesystem;


DatasetInfo::DatasetInfo(fs::path const &path, Options const &options)
    : definitionFile_{path},
      crossSection_{std::numeric_limits<double>::quiet_NaN()},
      numEventsTotal_{0},
      meanWeight_{std::numeric_limits<double>::quiet_NaN()} {

  LOG_DEBUG << "Constructing dataset from definition file " << definitionFile_
      << ".";

  if (options.GetConfig()["dataset_stems"]) {
    auto const &node = options.GetConfig()["dataset_stems"];
    auto const &filenames = node.as<std::vector<std::string>>();

    for (auto const &filename : filenames)
      stemsFiles_.emplace_back(FileInPath::Resolve(filename));
  }

  if (not fs::exists(path) or not fs::is_regular_file(path)) {
    HZZException exception;
    exception << "Dataset definition file " << path << " does not exist "
        "or is not a regular file.";
    throw exception;
  }

  ReadYaml(path);
}


YAML::Node const DatasetInfo::FindStem(std::string_view name) const {
  if (stemsFiles_.empty()) {
    throw HZZException{
      "Locations of files with dataset stems have not been set."};
  }

  for (auto const &stemsFile : stemsFiles_) {
    if (not fs::exists(stemsFile) or not fs::is_regular_file(stemsFile)) {
      HZZException exception;
      exception << "File with dataset stems " << stemsFile << " does not exist "
         "or is not a regular file.";
      throw exception;
    }

    auto stems = YAML::LoadFile(stemsFile);

    if (not stems.IsSequence() or stems.size() == 0) {
      HZZException exception;
      exception << "File with dataset stems " << stemsFile << " does not contain "
          "a sequence of stems or the sequence is empty.";
      throw exception;
    }

    for (auto const &stem : stems) {
      if (not stem["name"]) {
        HZZException exception;
        exception << "An entry in file with dataset stems " << stemsFile
            << " does not contain mandatory parameter \"name\".";
        throw exception;
      }

      if (stem["name"].as<std::string>() == name)
        return stem;
    }
  }

  // If this point is reached, the stem has not been found
  HZZException exception;
  exception << "No stem for name \"" << name
      << "\" was found in the following files:\n";
  for (auto const &stemsFile : stemsFiles_)
    exception << stemsFile << '\n';
  throw exception;
}


YAML::Node const DatasetInfo::GetNode(YAML::Node const root,
                                      std::string const &key) const {
  auto const node = root[key];

  if (not node or not node.IsScalar()) {
    HZZException exception;
    exception << "Dataset definition file " << definitionFile_
        << " does not contain mandatory scalar field \"" << key << "\".";
    throw exception;
  }

  return node;
}


void DatasetInfo::ReadYaml(fs::path const &path) {
  YAML::Node info = YAML::LoadFile(path);

  if (info["stem"])
    // This is not a full dataset definition file. Extend it with the stem.
    SpliceYaml(info);

  YAML::Node filesNode = info["files"];

  if (not filesNode or not filesNode.IsSequence()) {
    HZZException exception;
    exception << "Dataset definition file " << path
        << " does not contain mandatory sequence \"files\".";
    throw exception;
  }

  for (auto const &element : filesNode)
    files_.emplace_back(element.as<std::string>());

  LOG_TRACE << "Total number of paths of input files found in dataset "
      "definition file " << path << ": " << files_.size();

  if (files_.empty()) {
    HZZException exception;
    exception << "No paths to input files read from dataset definition file "
        << path << ".";
    throw exception;
  }


  // Save important parameters
  name_ = GetNode(info, "name").as<std::string>();
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
    : info_{std::move(info)}, chain_{"Events"} {

  if (skipFiles < 0) {
    HZZException exception;
    exception << "Illegal value " << skipFiles
        << " given to argument \"skipFiles\".";
    throw exception;
  }

  if (maxFiles < -1) {
    HZZException exception;
    exception << "Illegal value " << maxFiles
        << " given to argument \"maxFiles\".";
    throw exception;
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
