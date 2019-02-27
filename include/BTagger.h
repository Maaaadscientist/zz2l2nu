#ifndef btagger_h
#define btagger_h

#include <TLorentzVectorWithIndex.h>
#include <BTagCalibrationStandalone.h>
#include <Tables.h>
#include <SmartSelectionMonitor.h>
#include <SmartSelectionMonitor_hzz.h>

#include <iostream>
#include <string>
#include <TMath.h>
#include <vector>

namespace btagger
{

  BTagCalibrationReader loadCalibrationReader();

  utils::tables loadEffTables();

  double apply_sf(std::vector<TLorentzVectorWithIndex> selCentralJets, std::vector<double> btags, std::vector<float> *JetAk04HadFlav, const utils::tables &tab, BTagCalibrationReader _btag_calibration_reader, TString systName);

}

#endif
