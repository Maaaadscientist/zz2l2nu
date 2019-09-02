#include <cstdlib>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <InstrMetAnalysis.h>
#include <Logger.h>
#include <Looper.h>
#include <MainAnalysis.h>
#include <NrbAnalysis.h>
#include <Options.h>
#include <Version.h>

using namespace std;
namespace po = boost::program_options;


enum class AnalysisType {
  Main,
  InstrMET,
  NRB
};


int main(int argc, char **argv) {

  po::options_description optionsDescription{"Analysis"};
  optionsDescription.add_options()
    ("catalog",
     po::value<std::string>()->default_value(
       "/user/npostiau/event_files/MC_ewk/Bonzais-catalog_test_ZZTo2L2Nu-ZZ2l2vPruner.txt"),
     "Path to catalog file")
    ("max-events", po::value<int64_t>()->default_value(-1),
     "Maximal number of events to read; -1 means all")
    ("skip-files", po::value<int>()->default_value(0),
     "Number of files to skip at the beginning of the catalog")
    ("max-files", po::value<int>()->default_value(1),
     "Maximal number of files to read")
    ("analysis,a", po::value<std::string>()->default_value("Main"),
     "Analysis to run; allowed values are \"Main\", \"InstrMET\", \"NRB\"")
    ("dd-photon", "Use data-driven photon+jets background")
    ("syst", po::value<string>()->default_value(""),
     "Requested systematic variation")
    ("all-control-plots", "Keep all control plots")
    ("output,o", po::value<std::string>()->default_value("outputFile.root"),
     "Name for output file with histograms")
    ("seed", po::value<unsigned>()->default_value(0),
     "Seed for random number generator; 0 means a unique seed")
    ("mela-weight", po::value<unsigned>(), "MELA reweighting index");

  Options options(argc, argv, {optionsDescription});

  std::string analysisTypeArg{options.GetAs<std::string>("analysis")};
  boost::to_lower(analysisTypeArg);
  AnalysisType analysisType;

  if (analysisTypeArg == "main")
    analysisType = AnalysisType::Main;
  else if (analysisTypeArg == "instrmet")
    analysisType = AnalysisType::InstrMET;
  else if (analysisTypeArg == "nrb")
    analysisType = AnalysisType::NRB;
  else {
    LOG_ERROR << "Unknown analysis type \"" <<
      options.GetAs<std::string>("analysis") << "\"";
    std::exit(EXIT_FAILURE);
  }


  LOG_DEBUG << "Version: " << Version::Commit();
  
  switch (analysisType) {
    case AnalysisType::Main: {
      Looper<MainAnalysis> looper{options};
      looper.Run();
      break;
    }

    case AnalysisType::InstrMET: {
      Looper<InstrMetAnalysis> looper{options};
      looper.Run();
      break;
    }

    case AnalysisType::NRB: {
      Looper<NrbAnalysis> looper{options};
      looper.Run();
      break;
    }
  }
}

