#include <Options.h>

#include <cstdlib>
#include <iostream>

namespace po = boost::program_options;


Options::Options(int argc, char **argv,
                 std::initializer_list<Group> const &optionGroups)
    : programName_{argv[0]} {
  
  po::options_description generalOptions{"General"};
  generalOptions.add_options()("help,h", "Prints this help message");
  allOptions_.add(generalOptions);

  for (auto const &group : optionGroups)
    allOptions_.add(group);

  try {
    po::store(po::parse_command_line(argc, argv, allOptions_), optionMap_);
    po:notify(optionMap_);
  } catch (po::error const &e) {
    std::cerr << "Error while parsing command line arguments: " << e.what() <<
      "." << std::endl;
    PrintUsage();
    std::exit(EXIT_FAILURE);
  }

  if (optionMap_.count("help") > 0) {
    PrintUsage();
    std::exit(EXIT_SUCCESS);
  }
}


bool Options::Exists(std::string const &label) const {
  return (optionMap_.count(label) > 0);
}


void Options::PrintUsage() const {
  std::cerr << "Usage: " << programName_ << " [options]\n";
  std::cerr << allOptions_ << std::endl;
}

