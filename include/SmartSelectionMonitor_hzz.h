#ifndef smartselectionmonitor_hzz_hh
#define smartselectionmonitor_hzz_hh

#include <PhysicsObjects.h>
#include <SmartSelectionMonitor.h>
#include <Utils.h>

#include <type_traits>

struct base_evt{
  TString s_jetCat;
  TString s_lepCat;
  double MT;
  double M_Boson;
  double pT_Boson;
  double eta_Boson;
  double MET;
  double METphi;
  int nJets;
  int runNumber;
  int nVtx;
  double deltaPhi_MET_Jet;
  double deltaPhi_MET_Boson;
  double METoPT;
  double METpar;
  double METperp;
  //double MET_significance;
  double rho;
  double jet0_pT;
  double jet1_pT;
  double jet2_pT;
  double jet3_pT;
  double HT_selJets;

  void Fill_baseEvt(TString s_jetCat_, TString s_lepCat_, TLorentzVector boson_, TLorentzVector METVector_, std::vector<Jet> const &selJets_, double run_, double PV_npvsGood_, double fixedGridRhoFastjetAll_/*, double MET_significance_*/){ 
    s_jetCat = s_jetCat_;
    s_lepCat = s_lepCat_;
    MT = sqrt(pow(sqrt(pow(boson_.Pt(),2)+pow(boson_.M(),2))+sqrt(pow(METVector_.Pt(),2)+pow(91.1876,2)),2)-pow((boson_+METVector_).Pt(),2));;
    M_Boson = boson_.M();
    pT_Boson = boson_.Pt();
    eta_Boson = boson_.Eta();
    MET = METVector_.Pt();
    METphi = METVector_.Phi();
    nJets = selJets_.size();
    runNumber = run_;
    nVtx = PV_npvsGood_;
    double minDeltaPhiJetMET = 4.;
    for(int i = 0 ; i < selJets_.size() ; i++){
      if (fabs(utils::deltaPhi(selJets_[i].p4, METVector_)) < minDeltaPhiJetMET) minDeltaPhiJetMET = fabs(utils::deltaPhi(selJets_[i].p4, METVector_));
    }
    deltaPhi_MET_Jet = minDeltaPhiJetMET;
    deltaPhi_MET_Boson = fabs(utils::deltaPhi(boson_, METVector_));
    METoPT = METVector_.Pt()/(1.*boson_.Pt());
    double METorth_ =0, METpar_ = 0;
    TVector3 bosonDir = boson_.Vect();
    if (bosonDir.Perp()>0){
      TVector3 bosonDir2D(bosonDir.x()/bosonDir.Perp(), bosonDir.y()/bosonDir.Perp(), 0);
      TVector3 bosonPerp2D(-bosonDir2D.y(),bosonDir2D.x(),0);
      TVector3 METVector3D = METVector_.Vect();
      METpar_ = - METVector3D.Dot(bosonDir2D);
      METorth_ = METVector3D.Dot(bosonPerp2D);
    }
    METpar = METpar_;
    METperp = METorth_;
    //MET_significance = MET_significance_;
    rho = fixedGridRhoFastjetAll_;
    jet0_pT = ((selJets_.size() > 0) ? selJets_[0].p4.Pt() : 0);
    jet1_pT = ((selJets_.size() > 1) ? selJets_[1].p4.Pt() : 0);
    jet2_pT = ((selJets_.size() > 2) ? selJets_[2].p4.Pt() : 0);
    jet3_pT = ((selJets_.size() > 3) ? selJets_[3].p4.Pt() : 0);
    double allSelJets_HT =0.;
    for (unsigned int i =0; i < selJets_.size(); i++) allSelJets_HT += selJets_[i].p4.Pt();
    HT_selJets = allSelJets_HT;
  }

};

struct evt : base_evt{
  double lep1pT;
  double lep2pT;
  double lep1eta;
  double lep2eta;
  
  void Fill_evt(
      TString s_jetCat_, TString s_lepCat_, TLorentzVector boson_,
      TLorentzVector METVector_, std::vector<Jet> const &selJets_,
      double run_, double PV_npvsGood_, double fixedGridRhoFastjetAll_,
      /*double MET_significance_, */std::vector<Lepton> const &selLeptons_) {
    
    Fill_baseEvt(s_jetCat_, s_lepCat_, boson_, METVector_, selJets_, run_,
                 PV_npvsGood_, fixedGridRhoFastjetAll_/*, MET_significance_*/);

    if (selLeptons_.size() > 1) {
        lep1pT = selLeptons_[0].p4.Pt();
        lep2pT = selLeptons_[1].p4.Pt();
        lep1eta = selLeptons_[0].p4.Eta();
        lep2eta = selLeptons_[1].p4.Eta();
    }
  }
};

struct photon_evt : base_evt{
  double HoE;
  double sigmaIEtaIEta;
  double chIsoRhoCorr;
  double neuIsoRhoCorr;
  double phoIsoRhoCorr;
  double R9;

  void Fill_photonEvt(TString s_jetCat_, TString s_lepCat_, TLorentzVector boson_, TLorentzVector METVector_, std::vector<Jet> const &selJets_, double run_, double PV_npvsGood_, double fixedGridRhoFastjetAll_, /*double MET_significance_,*/ double HoE_ = -1., double sigmaIEtaIEta_ = -1., double chIsoRhoCorr_ = -1., double neuIsoRhoCorr_ = -1., double phoIsoRhoCorr_ = -1., double R9_ = -1.){
    Fill_baseEvt(s_jetCat_, s_lepCat_, boson_, METVector_, selJets_, run_, PV_npvsGood_, fixedGridRhoFastjetAll_/*, MET_significance_*/);
    HoE = HoE_;
    sigmaIEtaIEta = sigmaIEtaIEta_;
    chIsoRhoCorr = chIsoRhoCorr_;
    neuIsoRhoCorr = neuIsoRhoCorr_;
    phoIsoRhoCorr = phoIsoRhoCorr_;
    R9 = R9_;
  }


};

//Subclass of SmartSelectionMonitor, for methods related to the HZZ analysis only.

class SmartSelectionMonitor_hzz : public SmartSelectionMonitor {

public:

  //declares all histograms and set up the labels, bin names etc
  bool declareHistos();
  bool declareHistos_InstrMET();
  bool declareHistos_NRB();

  //fills a histogram for each category of jets and final states
  template<class T> bool fillHistoForAllCategories(TString name, double variable, T currentEvt, TString tag, double weight, bool divideByBinWidth = true);
  template<class T> bool fillProfileForAllCategories(TString name, double variableY, double variableX, T currentEvt, TString tag, double weight);

  //fills all analysis relevant histograms
  template<class T> bool fillAnalysisHistos_common(T currentEvt, TString tag, double weight, bool divideByBinWidth = true);
  bool fillAnalysisHistos(evt currentEvt, TString tag, double weight, bool divideByBinWidth = true);
  bool fillAnalysisHistos(photon_evt currentEvt, TString tag, double weight, bool divideByBinWidth = true);
  bool fillPhotonIDHistos_InstrMET(photon_evt currentEvt, TString tag, double weight, bool divideByBinWidth = true);
  bool fillInstrMETControlRegionHisto(base_evt currentEvt, TString tag, double weight, bool divideByBinWidth = true);

  //Method for writing histograms only when running on systematics. It appends each histogram name with the name of the systematic and, by default, only keeps the histograms with mT (but you may change it).
  void WriteForSysts(TString systName, bool keepEverything);

private:

};

#endif
