#ifndef smartselectionmonitor_hzz_hh
#define smartselectionmonitor_hzz_hh

#include "SmartSelectionMonitor.h"

#include "HZZ2l2nuLooper.h"

//Subclass of SmartSelectionMonitor, for methods related to the HZZ analysis only.

class SmartSelectionMonitor_hzz : public SmartSelectionMonitor {

public:

  //fills all analysis relevant histograms
  bool fillAnalysisHistos(evt currentEvt, TString tag, double weight);

private:

};

#endif
