#include <cstdlib>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <Logger.h>
#include <LooperMain.h>
#include <Options.h>
#include <Version.h>

using namespace std;
namespace po = boost::program_options;


enum class AnalysisType {
  Main,
  InstrMET,
  NRB,
  TnP
};


int main(int argc, char **argv) {

  po::options_description optionsDescription{"Analysis"};
  optionsDescription.add_options()
    ("catalog,c",
     po::value<std::string>()->default_value(
       "/user/npostiau/event_files/MC_ewk/Bonzais-catalog_test_ZZTo2L2Nu-ZZ2l2vPruner.txt"),
     "Path to catalog file")
    ("max-events", po::value<long long>()->default_value(-1),
     "Maximal number of events to read; -1 means all")
    ("skip-files", po::value<int>()->default_value(0),
     "Number of files to skip at the beginning of the catalog")
    ("max-files", po::value<int>()->default_value(1),
     "Maximal number of files to read")
    ("analysis,a", po::value<std::string>()->default_value("Main"),
     "Analysis to run; allowed values are \"Main\", \"InstrMET\", \"NRB\", "
     "\"TnP\"")
    ("dd-photon", "Use data-driven photon+jets background")
    ("is-mc", po::value<bool>()->default_value(true), "Simulation or real data")
    ("xsec", po::value<float>()->default_value(-1.), "Sample cross section, pb")
    ("syst", po::value<string>()->default_value(""),
     "Requested systematic variation")
    ("all-control-plots", "Keep all control plots")
    ("output,o", po::value<std::string>()->default_value("outputFile.root"),
     "Name for output file with histograms")
    ("seed", po::value<unsigned>()->default_value(0),
     "Seed for random number generator; 0 means a unique seed");

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
  else if (analysisTypeArg == "tnp")
    analysisType = AnalysisType::TnP;
  else {
    LOG_ERROR << "Unknown analysis type \"" <<
      options.GetAs<std::string>("analysis") << "\"";
    std::exit(EXIT_FAILURE);
  }


  LOG_DEBUG << "Version: " << Version::Commit();
  
  LooperMain myHZZlooper(options);

  switch (analysisType) {
    case AnalysisType::Main:
      myHZZlooper.Loop();
      break;

    case AnalysisType::InstrMET:
      myHZZlooper.Loop_InstrMET();
      break;

    case AnalysisType::NRB:
      myHZZlooper.Loop_NRB();
      break;

    case AnalysisType::TnP:
      myHZZlooper.Loop_TnP();
  }
}

