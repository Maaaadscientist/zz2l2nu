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
  double METphi;
  int nJets;
  double lep1pT;
  double lep2pT;
  double lep1eta;
  double lep2eta;
  int runNumber;
  int nVtx;
};

struct photon_evt{
  TString s_jetCat;
  TString s_lepCat;
  double MT;
  double M;
  double pT;
  double eta;
  double MET;
  double METphi;
  double deltaPhi_MET_jet;
  double deltaPhi_MET_phot;
  double METoPT;
  double METpar;
  double METperp;
  double METsig;
  int nJets;
  int nb_selPhoton;
  int nVtx;
  double HoE;
  double sigmaIEtaIEta;
  double chIsoRhoCorr;
  double neuIsoRhoCorr;
  double phoIsoRhoCorr;
  double R9;
  double rho;
  double jet0_pT;
  double jet1_pT;
  double jet2_pT;
  double jet3_pT;
  double HT_selJets;
  double HT_PFTightJetsAbove5GeV;
};

//Subclass of SmartSelectionMonitor, for methods related to the HZZ analysis only.

class SmartSelectionMonitor_hzz : public SmartSelectionMonitor {

public:

  //declares all histograms and set up the labels, bin names etc
  bool declareHistos();
  bool declareHistos_InstrMET();

  //fills a histogram for each category of jets and final states
  bool fillHistoForAllCategories(TString name, double data, evt currentEvt, TString tag, double weight);
  bool fillHistoForAllCategories(TString name, double data, photon_evt currentEvt, TString tag, double weight);

  //fills all analysis relevant histograms
  bool fillAnalysisHistos(evt currentEvt, TString tag, double weight);
  bool fillAnalysisHistos_InstrMET(photon_evt currentEvt, TString tag, double weight);
  bool fillPhotonIDHistos_InstrMET(photon_evt currentEvt, TString tag, double weight);

private:

};

#endif
