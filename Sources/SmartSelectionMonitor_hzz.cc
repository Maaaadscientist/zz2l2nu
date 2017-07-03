#include "../Includes/SmartSelectionMonitor_hzz.h"

bool SmartSelectionMonitor_hzz::fillAnalysisHistos(evt currentEvt, TString tag, double weight){
  fillHisto("mT", tag+currentEvt.s_lepCat+currentEvt.s_jetCat, currentEvt.transverseMass, weight);
  fillHisto("M_Z", tag+currentEvt.s_lepCat+currentEvt.s_jetCat, currentEvt.MZ, weight);
  fillHisto("pT_Z", tag+currentEvt.s_lepCat+currentEvt.s_jetCat, currentEvt.pTZ, weight);
  fillHisto("MET", tag+currentEvt.s_lepCat+currentEvt.s_jetCat, currentEvt.MET, weight);
  fillHisto("nJets", tag+currentEvt.s_lepCat+currentEvt.s_jetCat, currentEvt.nJets, weight);
  fillHisto("mT", tag+currentEvt.s_jetCat, currentEvt.transverseMass, weight);
  fillHisto("M_Z", tag+currentEvt.s_jetCat, currentEvt.MZ, weight);
  fillHisto("pT_Z", tag+currentEvt.s_jetCat, currentEvt.pTZ, weight);
  fillHisto("MET", tag+currentEvt.s_jetCat, currentEvt.MET, weight);
  fillHisto("nJets", tag+currentEvt.s_jetCat, currentEvt.nJets, weight);
  fillHisto("mT", tag+currentEvt.s_lepCat, currentEvt.transverseMass, weight);
  fillHisto("M_Z", tag+currentEvt.s_lepCat, currentEvt.MZ, weight);
  fillHisto("pT_Z", tag+currentEvt.s_lepCat, currentEvt.pTZ, weight);
  fillHisto("MET", tag+currentEvt.s_lepCat, currentEvt.MET, weight);
  fillHisto("nJets", tag+currentEvt.s_lepCat, currentEvt.nJets, weight);
  fillHisto("mT", tag, currentEvt.transverseMass, weight);
  fillHisto("M_Z", tag, currentEvt.MZ, weight);
  fillHisto("pT_Z", tag, currentEvt.pTZ, weight);
  fillHisto("MET", tag, currentEvt.MET, weight);
  fillHisto("nJets", tag, currentEvt.nJets, weight);
  return true;
}
