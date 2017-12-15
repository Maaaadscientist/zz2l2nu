#ifndef smartselectionmonitor_hzz_hh
#define smartselectionmonitor_hzz_hh

#include "SmartSelectionMonitor.h"

struct evt{
  TString s_jetCat;
  TString s_lepCat;
  double transverseMass;
  double MZ;
  double pTZ;
  double etaZ;
  double MET;
  int nJets;
};

//Subclass of SmartSelectionMonitor, for methods related to the HZZ analysis only.

class SmartSelectionMonitor_hzz : public SmartSelectionMonitor {

public:

  //declares all histograms and set up the labels, bin names etc
  bool declareHistos();
  bool declareHistos_InstrMET();

  //fills a histogram for each category of jets and final states
  bool fillHistoForAllCategories(TString name, double data, evt currentEvt, TString tag, double weight);

  //fills all analysis relevant histograms
  bool fillAnalysisHistos(evt currentEvt, TString tag, double weight);
  bool fillAnalysisHistos_InstrMET(evt currentEvt, TString tag, double weight);

private:

};

#endif
