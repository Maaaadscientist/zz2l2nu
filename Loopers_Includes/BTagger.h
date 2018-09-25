#ifndef btagger_h
#define btagger_h

#include <iostream>
#include <string>
#include <TMath.h>
#include <vector>
#include "../Common/TLorentzVectorWithIndex.h"
#include "../Common/BTagCalibrationStandalone.h"
#include "../Common/Tables.h"
#include "SmartSelectionMonitor.h"
#include "SmartSelectionMonitor_hzz.h"

namespace btagger
{

  BTagCalibrationReader loadCalibrationReader();

  utils::tables loadEffTables();

  bool fill_eff(std::vector<TLorentzVectorWithIndex> selCentralJets, std::vector<double> btags, std::vector<float> *JetAk04HadFlav, double weight, SmartSelectionMonitor_hzz & mon);

  double apply_sf(std::vector<TLorentzVectorWithIndex> selCentralJets, std::vector<double> btags, std::vector<float> *JetAk04HadFlav, const utils::tables &tab, BTagCalibrationReader _btag_calibration_reader);

}

#endif
