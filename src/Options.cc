#include <Options.h>

#include <cstdlib>
#include <filesystem>
#include <iostream>

#include <FileInPath.h>
#include <Version.h>

namespace po = boost::program_options;


Options::Options(int argc, char **argv,
                 std::initializer_list<Group> const &optionGroups)
    : programName_{argv[0]} {
  
  po::options_description generalOptions{"General"};
  generalOptions.add_options()
    ("help,h", "Prints this help message")
    ("config", po::value<std::filesystem::path>(), "YAML configuration file")
    ("verbosity,v", po::value<int>()->default_value(1),
     "Verbosity level: warnings and errors (0), info (1), debug (2), "
     "trace (3)")
    ("version", "Prints the version and exits");
  allOptions_.add(generalOptions);

  for (auto const &group : optionGroups)
    allOptions_.add(group);

  try {
    po::store(po::parse_command_line(argc, argv, allOptions_), optionMap_);
    po:notify(optionMap_);
  } catch (po::error const &e) {
    LOG_ERROR << "Error while parsing command line arguments: " << e.what() <<
      ".";
    PrintUsage();
    std::exit(EXIT_FAILURE);
  }

  if (optionMap_.count("help") > 0) {
    PrintUsage();
    std::exit(EXIT_FAILURE);
  }

  if (optionMap_.count("version") > 0) {
    std::cerr << Version::Commit() << std::endl;
    std::exit(EXIT_FAILURE);
  }


  int const verbosityLevel = GetAsChecked<int>("verbosity",
    [](int level){return (level >= 0 and level <= 3);});
  Logger::SeverityLevel severityLevel;

  switch (verbosityLevel) {
    case 0:
      severityLevel = Logger::SeverityLevel::kWarning;
      break;

    case 1:
      severityLevel = Logger::SeverityLevel::kInfo;
      break;

    case 2:
      severityLevel = Logger::SeverityLevel::kDebug;
      break;

    case 3:
      severityLevel = Logger::SeverityLevel::kTrace;
      break;
  };

  Logger::SetLevel(severityLevel);


  if (optionMap_.count("config") > 0) {
    config_ = YAML::LoadFile(FileInPath::Resolve(
      GetAs<std::filesystem::path>("config")));
  }
}


YAML::Node const &Options::GetConfig() const {
  if (config_.IsNull())
    throw Error("Configuration is not available.");
  else
    return config_;
}


bool Options::Exists(std::string const &label) const {
  return (optionMap_.count(label) > 0);
}


void Options::PrintUsage() const {
  std::cerr << "Usage: " << programName_ << " [options]\n";
  std::cerr << allOptions_ << std::endl;
}

