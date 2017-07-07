#ifndef smartselectionmonitor_hzz_hh
#define smartselectionmonitor_hzz_hh

#include "SmartSelectionMonitor.h"

#include "HZZ2l2nuLooper.h"

//Subclass of SmartSelectionMonitor, for methods related to the HZZ analysis only.

class SmartSelectionMonitor_hzz : public SmartSelectionMonitor {

public:

  //declares all histograms and set up the labels, bin names etc
  bool declareHistos();

  //fills a histogram for each category of jets and final states
  bool fillHistoForAllCategories(TString name, double data, evt currentEvt, TString tag, double weight);

  //fills all analysis relevant histograms
  bool fillAnalysisHistos(evt currentEvt, TString tag, double weight);

private:

};

#endif
