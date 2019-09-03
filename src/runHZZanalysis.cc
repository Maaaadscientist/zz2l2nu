#include <cstdlib>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <DileptonTrees.h>
#include <InstrMetAnalysis.h>
#include <Logger.h>
#include <Looper.h>
#include <MainAnalysis.h>
#include <NrbAnalysis.h>
#include <Options.h>
#include <Version.h>

namespace po = boost::program_options;


enum class AnalysisType {
  Main,
  InstrMET,
  NRB,
  DileptonTrees
};


template<typename T>
void runAnalysis(int argc, char **argv,
                 po::options_description const &commonOptions) {
  Options options(
      argc, argv, {commonOptions, Looper<T>::OptionsDescription()});
  LOG_DEBUG << "Version: " << Version::Commit();
  Looper<T> looper{options};
  looper.Run();
}


int main(int argc, char **argv) {
  po::options_description analysisTypeOptions{"Analysis type"};
  analysisTypeOptions.add_options()
    ("analysis,a", po::value<std::string>()->default_value("Main"),
     "Analysis to run; allowed values are \"Main\", \"InstrMET\", \"NRB\", "
     "\"DileptonTrees\"");

  // Command line options are checked twice. At the first pass only check the
  // analysis type and update the list of expected options accordingly.
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(analysisTypeOptions)\
      .allow_unregistered().run(), vm);
  po::notify(vm);

  std::string analysisTypeArg{vm["analysis"].as<std::string>()};
  boost::to_lower(analysisTypeArg);
  AnalysisType analysisType;

  if (analysisTypeArg == "main")
    analysisType = AnalysisType::Main;
  else if (analysisTypeArg == "instrmet")
    analysisType = AnalysisType::InstrMET;
  else if (analysisTypeArg == "nrb")
    analysisType = AnalysisType::NRB;
  else if (analysisTypeArg == "dileptontrees")
    analysisType = AnalysisType::DileptonTrees;
  else {
    LOG_ERROR << "Unknown analysis type \"" <<
      vm["analysis"].as<std::string>() << "\"";
    std::exit(EXIT_FAILURE);
  }


  // Perform the second pass over command line options and run the requested
  // analysis
  switch (analysisType) {
    case AnalysisType::Main:
      runAnalysis<MainAnalysis>(argc, argv, analysisTypeOptions);
      break;

    case AnalysisType::InstrMET:
      runAnalysis<InstrMetAnalysis>(argc, argv, analysisTypeOptions);
      break;

    case AnalysisType::NRB:
      runAnalysis<NrbAnalysis>(argc, argv, analysisTypeOptions);
      break;

    case AnalysisType::DileptonTrees:
      runAnalysis<DileptonTrees>(argc, argv, analysisTypeOptions);
      break;
  }
}

