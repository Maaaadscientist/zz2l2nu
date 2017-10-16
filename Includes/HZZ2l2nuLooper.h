//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 23 18:51:00 2017 by ROOT version 6.06/01
// from TChain tupel/EventTree/
//////////////////////////////////////////////////////////

#ifndef HZZ2l2nuLooper_h
#define HZZ2l2nuLooper_h

#include <iostream>
#include <fstream>
#include <thread>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

using std::vector;
using namespace std;

struct evt{
  TString s_jetCat;
  TString s_lepCat;
  double transverseMass;
  double MZ;
  double pTZ;
  double MET;
  int nJets;
};

class HZZ2l2nuLooper {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   //global variables
   int maxEvents_;
   TString outputFile_;
   int isMC_;
   double sampleXsection_;
   double totalEventsInBaobab_;
   double sumWeightInBaobab_;
   double sumWeightInBonzai_;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           EvtIsRealData;
   UInt_t          EvtNum;
   UInt_t          EvtRunNum;
   Int_t           EvtLumiNum;
   Int_t           EvtBxNum;
   Int_t           EvtVtxCnt;
   Int_t           EvtPuCnt;
   Int_t           EvtPuCntTruth;
   vector<double>  *EvtWeights;
   Float_t         EvtFastJetRho;
   ULong64_t       TrigMET;
   UInt_t          TrigHlt;
   ULong64_t       TrigHltPhot;
   ULong64_t       TrigHltDiPhot;
   ULong64_t       TrigHltMu;
   ULong64_t       TrigHltDiMu;
   ULong64_t       TrigHltEl;
   ULong64_t       TrigHltDiEl;
   ULong64_t       TrigHltElMu;
   vector<float>   *METPt;
   vector<float>   *METPx;
   vector<float>   *METPy;
   vector<float>   *METPz;
   vector<float>   *METE;
   vector<float>   *METsigx2;
   vector<float>   *METsigxy;
   vector<float>   *METsigy2;
   vector<float>   *METsig;
   vector<float>   *GMETPt;
   vector<float>   *GMETPx;
   vector<float>   *GMETPy;
   vector<float>   *GMETPz;
   vector<float>   *GMETE;
   vector<float>   *GLepDr01Pt;
   vector<float>   *GLepDr01Eta;
   vector<float>   *GLepDr01Phi;
   vector<float>   *GLepDr01E;
   vector<int>     *GLepDr01Id;
   vector<int>     *GLepDr01St;
   vector<int>     *GLepDr01MomId;
   vector<bool>    *GLepDr01Prompt;
   vector<bool>    *GLepDr01TauProd;
   vector<float>   *GLepBarePt;
   vector<float>   *GLepBareEta;
   vector<float>   *GLepBarePhi;
   vector<float>   *GLepBareE;
   vector<int>     *GLepBareId;
   vector<int>     *GLepBareSt;
   vector<int>     *GLepBareMomId;
   vector<bool>    *GLepBarePrompt;
   vector<bool>    *GLepBareTauProd;
   vector<float>   *GLepSt3Pt;
   vector<float>   *GLepSt3Eta;
   vector<float>   *GLepSt3Phi;
   vector<float>   *GLepSt3E;
   vector<int>     *GLepSt3Id;
   vector<int>     *GLepSt3St;
   vector<int>     *GLepSt3Mother0Id;
   vector<int>     *GLepSt3MotherCnt;
   vector<float>   *GPhotPt;
   vector<float>   *GPhotEta;
   vector<float>   *GPhotPhi;
   vector<float>   *GPhotE;
   vector<int>     *GPhotMotherId;
   vector<int>     *GPhotSt;
   vector<float>   *GPhotIsoEDR03;
   vector<float>   *GPhotIsoEDR04;
   vector<float>   *GPhotIsoEDR05;
   vector<float>   *GPhotIsoSumPtDR03;
   vector<float>   *GPhotIsoSumPtDR04;
   vector<float>   *GPhotIsoSumPtDR05;
   vector<float>   *GLepClosePhotPt;
   vector<float>   *GLepClosePhotEta;
   vector<float>   *GLepClosePhotPhi;
   vector<float>   *GLepClosePhotE;
   vector<int>     *GLepClosePhotId;
   vector<int>     *GLepClosePhotMother0Id;
   vector<int>     *GLepClosePhotMotherCnt;
   vector<int>     *GLepClosePhotSt;
   vector<float>   *GJetAk04Pt;
   vector<float>   *GJetAk04Eta;
   vector<float>   *GJetAk04Phi;
   vector<float>   *GJetAk04E;
   vector<float>   *GJetAk04ChFrac;
   vector<int>     *GJetAk04ConstCnt;
   vector<int>     *GJetAk04ConstId;
   vector<float>   *GJetAk04ConstPt;
   vector<float>   *GJetAk04ConstEta;
   vector<float>   *GJetAk04ConstPhi;
   vector<float>   *GJetAk04ConstE;
   vector<float>   *GJetAk04MatchedPartonID;
   vector<float>   *GJetAk04MatchedPartonDR;
   vector<float>   *GJetAk08Pt;
   vector<float>   *GJetAk08Eta;
   vector<float>   *GJetAk08Phi;
   vector<float>   *GJetAk08E;
   vector<float>   *GJetAk08ChFrac;
   vector<int>     *GJetAk08ConstCnt;
   vector<int>     *GJetAk08ConstId;
   vector<float>   *GJetAk08ConstPt;
   vector<float>   *GJetAk08ConstEta;
   vector<float>   *GJetAk08ConstPhi;
   vector<float>   *GJetAk08ConstE;
   vector<float>   *GJetAk08MatchedPartonID;
   vector<float>   *GJetAk08MatchedPartonDR;
   vector<int>     *GPdfId1;
   vector<int>     *GPdfId2;
   vector<float>   *GPdfx1;
   vector<float>   *GPdfx2;
   vector<float>   *GPdfScale;
   Float_t         GBinningValue;
   Int_t           GNup;
   vector<float>   *MuPt;
   vector<float>   *MuEta;
   vector<float>   *MuPhi;
   vector<float>   *MuE;
   vector<unsigned int> *MuId;
   vector<unsigned int> *MuIdTight;
   vector<float>   *MuCh;
   vector<float>   *MuVtxZ;
   vector<float>   *MuDxy;
   vector<float>   *MuIsoRho;
   vector<float>   *MuPfIso;
   vector<float>   *MuType;
   vector<float>   *MuIsoTkIsoAbs;
   vector<float>   *MuIsoTkIsoRel;
   vector<float>   *MuIsoCalAbs;
   vector<float>   *MuIsoCombRel;
   vector<float>   *MuTkNormChi2;
   vector<int>     *MuTkHitCnt;
   vector<int>     *MuMatchedStationCnt;
   vector<float>   *MuDz;
   vector<int>     *MuPixelHitCnt;
   vector<int>     *MuTkLayerCnt;
   vector<float>   *MuPfIsoChHad;
   vector<float>   *MuPfIsoNeutralHad;
   vector<float>   *MuPfIsoRawRel;
   vector<unsigned int> *MuHltMatch;
   vector<float>   *ElPt;
   vector<float>   *ElEta;
   vector<float>   *ElEtaSc;
   vector<float>   *ElPhi;
   vector<float>   *ElE;
   vector<unsigned int> *ElId;
   vector<float>   *ElCh;
   vector<float>   *ElScRawE;
   vector<float>   *ElCorrE;
   vector<float>   *ElEcalIso;
   vector<float>   *ElEcalPfIso;
   vector<float>   *ElMvaTrig;
   vector<float>   *ElMvaNonTrig;
   vector<float>   *ElMvaPresel;
   vector<float>   *ElDEtaTkScAtVtx;
   vector<float>   *ElDPhiTkScAtVtx;
   vector<float>   *ElHoE;
   vector<float>   *ElSigmaIetaIeta;
   vector<float>   *ElSigmaIetaIetaFull5x5;
   vector<float>   *ElEinvMinusPinv;
   vector<float>   *ElD0;
   vector<float>   *ElDz;
   vector<int>     *ElExpectedMissingInnerHitCnt;
   vector<int>     *ElPassConvVeto;
   vector<unsigned int> *ElHltMatch;
   vector<float>   *ElPfIsoChHad;
   vector<float>   *ElPfIsoNeutralHad;
   vector<float>   *ElPfIsoIso;
   vector<float>   *ElPfIsoPuChHad;
   vector<float>   *ElPfIsoRaw;
   vector<float>   *ElPfIsoDbeta;
   vector<float>   *ElPfIsoRho;
   vector<float>   *ElAEff;
   vector<float>   *ElDr03TkSumPt;
   vector<float>   *ElDr03EcalRecHitSumEt;
   vector<float>   *ElDr03HcalTowerSumEt;
   vector<float>   *TauPt;
   vector<float>   *TauEta;
   vector<float>   *TauPhi;
   vector<float>   *TauE;
   vector<float>   *TauCh;
   vector<unsigned int> *TauDecayModeFinding;
   vector<float>   *TauCombinedIsolationDeltaBetaCorrRaw3Hits;
   vector<unsigned int> *TauDiscMuonLoose;
   vector<unsigned int> *TauDiscMuonTight;
   vector<unsigned int> *TauDiscElVLoose;
   vector<unsigned int> *TauDiscElLoose;
   vector<unsigned int> *TauDiscElVTight;
   vector<unsigned int> *TauDiscElTight;
   vector<float>   *PhotPt;
   vector<float>   *PhotEta;
   vector<float>   *PhotPhi;
   vector<float>   *PhotScRawE;
   vector<float>   *PhotScEta;
   vector<float>   *PhotScPhi;
   vector<float>   *PhotIsoEcal;
   vector<float>   *PhotIsoHcal;
   vector<float>   *PhotIsoTk;
   vector<float>   *PhotPfIsoChHad;
   vector<float>   *PhotPfIsoNeutralHad;
   vector<float>   *PhotPfIsoPhot;
   vector<float>   *PhotPfIsoPuChHad;
   vector<float>   *PhotPfIsoEcalClus;
   vector<float>   *PhotPfIsoHcalClus;
   vector<float>   *PhotE3x3;
   vector<float>   *PhotE1x5;
   vector<float>   *PhotE2x5;
   vector<float>   *PhotE5x5;
   vector<float>   *PhotSigmaIetaIeta;
   vector<float>   *PhotHoE;
   vector<float>   *PhotEtaWidth;
   vector<float>   *PhotPhiWidth;
   vector<float>   *PhotR9;
   vector<float>   *PhotE1x3;
   vector<float>   *PhotE2x2;
   vector<float>   *PhotS4;
   vector<float>   *PhotE1x5Full5x5;
   vector<float>   *PhotE2x5Full5x5;
   vector<float>   *PhotE3x3Full5x5;
   vector<float>   *PhotE5x5Full5x5;
   vector<float>   *PhotSigmaIetaIetaFull5x5;
   vector<float>   *PhotR9Full5x5;
   vector<unsigned int> *PhotId;
   vector<bool>    *PhotHasPixelSeed;
   vector<int>     *PhotPassElVeto;
   vector<float>   *JetAk04Pt;
   vector<float>   *JetAk04Eta;
   vector<float>   *JetAk04Phi;
   vector<float>   *JetAk04E;
   vector<float>   *JetAk04Id;
   vector<bool>    *JetAk04PuId;
   vector<float>   *JetAk04PuMva;
   vector<float>   *JetAk04RawPt;
   vector<float>   *JetAk04RawE;
   vector<float>   *JetAk04HfHadE;
   vector<float>   *JetAk04HfEmE;
   vector<float>   *JetAk04ChHadFrac;
   vector<float>   *JetAk04NeutralHadAndHfFrac;
   vector<float>   *JetAk04ChEmFrac;
   vector<float>   *JetAk04NeutralEmFrac;
   vector<float>   *JetAk04ChMult;
   vector<float>   *JetAk04ConstCnt;
   vector<float>   *JetAk04Beta;
   vector<float>   *JetAk04BetaClassic;
   vector<float>   *JetAk04BetaStar;
   vector<float>   *JetAk04BetaStarClassic;
   vector<float>   *JetAk04Rms;
   vector<float>   *JetAk04BTagCsv;
   vector<float>   *JetAk04BTagCsvV1;
   vector<float>   *JetAk04BTagCsvSLV1;
   vector<float>   *JetAk04BDiscCisvV2;
   vector<float>   *JetAk04BDiscJp;
   vector<float>   *JetAk04BDiscBjp;
   vector<float>   *JetAk04BDiscTche;
   vector<float>   *JetAk04BDiscTchp;
   vector<float>   *JetAk04BDiscSsvhe;
   vector<float>   *JetAk04BDiscSsvhp;
   vector<float>   *JetAk04PartFlav;
   vector<float>   *JetAk04HadFlav;
   vector<float>   *JetAk04JecUncUp;
   vector<float>   *JetAk04JecUncDwn;
   vector<int>     *JetAk04ConstId;
   vector<float>   *JetAk04ConstPt;
   vector<float>   *JetAk04ConstEta;
   vector<float>   *JetAk04ConstPhi;
   vector<float>   *JetAk04ConstE;
   vector<int>     *JetAk04GenJet;
   vector<float>   *JetAk08Pt;
   vector<float>   *JetAk08Eta;
   vector<float>   *JetAk08Phi;
   vector<float>   *JetAk08E;
   vector<float>   *JetAk08Id;
   vector<float>   *JetAk08RawPt;
   vector<float>   *JetAk08RawE;
   vector<float>   *JetAk08HfHadE;
   vector<float>   *JetAk08HfEmE;
   vector<float>   *JetAk08ChHadFrac;
   vector<float>   *JetAk08NeutralHadAndHfFrac;
   vector<float>   *JetAk08ChEmFrac;
   vector<float>   *JetAk08NeutralEmFrac;
   vector<float>   *JetAk08ChMult;
   vector<float>   *JetAk08ConstCnt;
   vector<float>   *JetAk08BTagCsv;
   vector<float>   *JetAk08BTagCsvV1;
   vector<float>   *JetAk08BTagCsvSLV1;
   vector<float>   *JetAk08BDiscCisvV2;
   vector<float>   *JetAk08BDiscJp;
   vector<float>   *JetAk08BDiscBjp;
   vector<float>   *JetAk08BDiscTche;
   vector<float>   *JetAk08BDiscTchp;
   vector<float>   *JetAk08BDiscSsvhe;
   vector<float>   *JetAk08BDiscSsvhp;
   vector<float>   *JetAk08PartFlav;
   vector<float>   *JetAk08HadFlav;
   vector<float>   *JetAk08JecUncUp;
   vector<float>   *JetAk08JecUncDwn;
   vector<int>     *JetAk08ConstId;
   vector<float>   *JetAk08ConstPt;
   vector<float>   *JetAk08ConstEta;
   vector<float>   *JetAk08ConstPhi;
   vector<float>   *JetAk08ConstE;
   vector<int>     *JetAk08GenJet;
   vector<float>   *JetAk08PrunedMass;
   vector<float>   *JetAk08FilteredMass;
   vector<float>   *JetAk08SoftDropMass;
   vector<float>   *JetAk08TrimmedMass;
   vector<float>   *JetAk08Tau1;
   vector<float>   *JetAk08Tau2;
   vector<float>   *JetAk08Tau3;

   // List of branches
   TBranch        *b_EvtIsRealData;   //!
   TBranch        *b_EvtNum;   //!
   TBranch        *b_EvtRunNum;   //!
   TBranch        *b_EvtLumiNum;   //!
   TBranch        *b_EvtBxNum;   //!
   TBranch        *b_EvtVtxCnt;   //!
   TBranch        *b_EvtPuCnt;   //!
   TBranch        *b_EvtPuCntTruth;   //!
   TBranch        *b_EvtWeights;   //!
   TBranch        *b_EvtFastJetRho;   //!
   TBranch        *b_TrigMET;   //!
   TBranch        *b_TrigHlt;   //!
   TBranch        *b_TrigHltPhot;   //!
   TBranch        *b_TrigHltDiPhot;   //!
   TBranch        *b_TrigHltMu;   //!
   TBranch        *b_TrigHltDiMu;   //!
   TBranch        *b_TrigHltEl;   //!
   TBranch        *b_TrigHltDiEl;   //!
   TBranch        *b_TrigHltElMu;   //!
   TBranch        *b_METPt;   //!
   TBranch        *b_METPx;   //!
   TBranch        *b_METPy;   //!
   TBranch        *b_METPz;   //!
   TBranch        *b_METE;   //!
   TBranch        *b_METsigx2;   //!
   TBranch        *b_METsigxy;   //!
   TBranch        *b_METsigy2;   //!
   TBranch        *b_METsig;   //!
   TBranch        *b_GMETPt;   //!
   TBranch        *b_GMETPx;   //!
   TBranch        *b_GMETPy;   //!
   TBranch        *b_GMETPz;   //!
   TBranch        *b_GMETE;   //!
   TBranch        *b_GLepDr01Pt;   //!
   TBranch        *b_GLepDr01Eta;   //!
   TBranch        *b_GLepDr01Phi;   //!
   TBranch        *b_GLepDr01E;   //!
   TBranch        *b_GLepDr01Id;   //!
   TBranch        *b_GLepDr01St;   //!
   TBranch        *b_GLepDr01MomId;   //!
   TBranch        *b_GLepDr01Prompt;   //!
   TBranch        *b_GLepDr01TauProd;   //!
   TBranch        *b_GLepBarePt;   //!
   TBranch        *b_GLepBareEta;   //!
   TBranch        *b_GLepBarePhi;   //!
   TBranch        *b_GLepBareE;   //!
   TBranch        *b_GLepBareId;   //!
   TBranch        *b_GLepBareSt;   //!
   TBranch        *b_GLepBareMomId;   //!
   TBranch        *b_GLepBarePrompt;   //!
   TBranch        *b_GLepBareTauProd;   //!
   TBranch        *b_GLepSt3Pt;   //!
   TBranch        *b_GLepSt3Eta;   //!
   TBranch        *b_GLepSt3Phi;   //!
   TBranch        *b_GLepSt3E;   //!
   TBranch        *b_GLepSt3Id;   //!
   TBranch        *b_GLepSt3St;   //!
   TBranch        *b_GLepSt3Mother0Id;   //!
   TBranch        *b_GLepSt3MotherCnt;   //!
   TBranch        *b_GPhotPt;   //!
   TBranch        *b_GPhotEta;   //!
   TBranch        *b_GPhotPhi;   //!
   TBranch        *b_GPhotE;   //!
   TBranch        *b_GPhotMotherId;   //!
   TBranch        *b_GPhotSt;   //!
   TBranch        *b_GPhotIsoEDR03;   //!
   TBranch        *b_GPhotIsoEDR04;   //!
   TBranch        *b_GPhotIsoEDR05;   //!
   TBranch        *b_GPhotIsoSumPtDR03;   //!
   TBranch        *b_GPhotIsoSumPtDR04;   //!
   TBranch        *b_GPhotIsoSumPtDR05;   //!
   TBranch        *b_GLepClosePhotPt;   //!
   TBranch        *b_GLepClosePhotEta;   //!
   TBranch        *b_GLepClosePhotPhi;   //!
   TBranch        *b_GLepClosePhotE;   //!
   TBranch        *b_GLepClosePhotId;   //!
   TBranch        *b_GLepClosePhotMother0Id;   //!
   TBranch        *b_GLepClosePhotMotherCnt;   //!
   TBranch        *b_GLepClosePhotSt;   //!
   TBranch        *b_GJetAk04Pt;   //!
   TBranch        *b_GJetAk04Eta;   //!
   TBranch        *b_GJetAk04Phi;   //!
   TBranch        *b_GJetAk04E;   //!
   TBranch        *b_GJetAk04ChFrac;   //!
   TBranch        *b_GJetAk04ConstCnt;   //!
   TBranch        *b_GJetAk04ConstId;   //!
   TBranch        *b_GJetAk04ConstPt;   //!
   TBranch        *b_GJetAk04ConstEta;   //!
   TBranch        *b_GJetAk04ConstPhi;   //!
   TBranch        *b_GJetAk04ConstE;   //!
   TBranch        *b_GJetAk04MatchedPartonID;   //!
   TBranch        *b_GJetAk04MatchedPartonDR;   //!
   TBranch        *b_GJetAk08Pt;   //!
   TBranch        *b_GJetAk08Eta;   //!
   TBranch        *b_GJetAk08Phi;   //!
   TBranch        *b_GJetAk08E;   //!
   TBranch        *b_GJetAk08ChFrac;   //!
   TBranch        *b_GJetAk08ConstCnt;   //!
   TBranch        *b_GJetAk08ConstId;   //!
   TBranch        *b_GJetAk08ConstPt;   //!
   TBranch        *b_GJetAk08ConstEta;   //!
   TBranch        *b_GJetAk08ConstPhi;   //!
   TBranch        *b_GJetAk08ConstE;   //!
   TBranch        *b_GJetAk08MatchedPartonID;   //!
   TBranch        *b_GJetAk08MatchedPartonDR;   //!
   TBranch        *b_GPdfId1;   //!
   TBranch        *b_GPdfId2;   //!
   TBranch        *b_GPdfx1;   //!
   TBranch        *b_GPdfx2;   //!
   TBranch        *b_GPdfScale;   //!
   TBranch        *b_GBinningValue;   //!
   TBranch        *b_GNup;   //!
   TBranch        *b_MuPt;   //!
   TBranch        *b_MuEta;   //!
   TBranch        *b_MuPhi;   //!
   TBranch        *b_MuE;   //!
   TBranch        *b_MuId;   //!
   TBranch        *b_MuIdTight;   //!
   TBranch        *b_MuCh;   //!
   TBranch        *b_MuVtxZ;   //!
   TBranch        *b_MuDxy;   //!
   TBranch        *b_MuIsoRho;   //!
   TBranch        *b_MuPfIso;   //!
   TBranch        *b_MuType;   //!
   TBranch        *b_MuIsoTkIsoAbs;   //!
   TBranch        *b_MuIsoTkIsoRel;   //!
   TBranch        *b_MuIsoCalAbs;   //!
   TBranch        *b_MuIsoCombRel;   //!
   TBranch        *b_MuTkNormChi2;   //!
   TBranch        *b_MuTkHitCnt;   //!
   TBranch        *b_MuMatchedStationCnt;   //!
   TBranch        *b_MuDz;   //!
   TBranch        *b_MuPixelHitCnt;   //!
   TBranch        *b_MuTkLayerCnt;   //!
   TBranch        *b_MuPfIsoChHad;   //!
   TBranch        *b_MuPfIsoNeutralHad;   //!
   TBranch        *b_MuPfIsoRawRel;   //!
   TBranch        *b_MuHltMatch;   //!
   TBranch        *b_ElPt;   //!
   TBranch        *b_ElEta;   //!
   TBranch        *b_ElEtaSc;   //!
   TBranch        *b_ElPhi;   //!
   TBranch        *b_ElE;   //!
   TBranch        *b_ElId;   //!
   TBranch        *b_ElCh;   //!
   TBranch        *b_ElScRawE;   //!
   TBranch        *b_ElCorrE;   //!
   TBranch        *b_ElEcalIso;   //!
   TBranch        *b_ElEcalPfIso;   //!
   TBranch        *b_ElMvaTrig;   //!
   TBranch        *b_ElMvaNonTrig;   //!
   TBranch        *b_ElMvaPresel;   //!
   TBranch        *b_ElDEtaTkScAtVtx;   //!
   TBranch        *b_ElDPhiTkScAtVtx;   //!
   TBranch        *b_ElHoE;   //!
   TBranch        *b_ElSigmaIetaIeta;   //!
   TBranch        *b_ElSigmaIetaIetaFull5x5;   //!
   TBranch        *b_ElEinvMinusPinv;   //!
   TBranch        *b_ElD0;   //!
   TBranch        *b_ElDz;   //!
   TBranch        *b_ElExpectedMissingInnerHitCnt;   //!
   TBranch        *b_ElPassConvVeto;   //!
   TBranch        *b_ElHltMatch;   //!
   TBranch        *b_ElPfIsoChHad;   //!
   TBranch        *b_ElPfIsoNeutralHad;   //!
   TBranch        *b_ElPfIsoIso;   //!
   TBranch        *b_ElPfIsoPuChHad;   //!
   TBranch        *b_ElPfIsoRaw;   //!
   TBranch        *b_ElPfIsoDbeta;   //!
   TBranch        *b_ElPfIsoRho;   //!
   TBranch        *b_ElAEff;   //!
   TBranch        *b_ElDr03TkSumPt;   //!
   TBranch        *b_ElDr03EcalRecHitSumEt;   //!
   TBranch        *b_ElDr03HcalTowerSumEt;   //!
   TBranch        *b_TauPt;   //!
   TBranch        *b_TauEta;   //!
   TBranch        *b_TauPhi;   //!
   TBranch        *b_TauE;   //!
   TBranch        *b_TauCh;   //!
   TBranch        *b_TauDecayModeFinding;   //!
   TBranch        *b_TauCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_TauDiscMuonLoose;   //!
   TBranch        *b_TauDiscMuonTight;   //!
   TBranch        *b_TauDiscElVLoose;   //!
   TBranch        *b_TauDiscElLoose;   //!
   TBranch        *b_TauDiscElVTight;   //!
   TBranch        *b_TauDiscElTight;   //!
   TBranch        *b_PhotPt;   //!
   TBranch        *b_PhotEta;   //!
   TBranch        *b_PhotPhi;   //!
   TBranch        *b_PhotScRawE;   //!
   TBranch        *b_PhotScEta;   //!
   TBranch        *b_PhotScPhi;   //!
   TBranch        *b_PhotIsoEcal;   //!
   TBranch        *b_PhotIsoHcal;   //!
   TBranch        *b_PhotIsoTk;   //!
   TBranch        *b_PhotPfIsoChHad;   //!
   TBranch        *b_PhotPfIsoNeutralHad;   //!
   TBranch        *b_PhotPfIsoPhot;   //!
   TBranch        *b_PhotPfIsoPuChHad;   //!
   TBranch        *b_PhotPfIsoEcalClus;   //!
   TBranch        *b_PhotPfIsoHcalClus;   //!
   TBranch        *b_PhotE3x3;   //!
   TBranch        *b_PhotE1x5;   //!
   TBranch        *b_PhotE2x5;   //!
   TBranch        *b_PhotE5x5;   //!
   TBranch        *b_PhotSigmaIetaIeta;   //!
   TBranch        *b_PhotHoE;   //!
   TBranch        *b_PhotEtaWidth;   //!
   TBranch        *b_PhotPhiWidth;   //!
   TBranch        *b_PhotR9;   //!
   TBranch        *b_PhotE1x3;   //!
   TBranch        *b_PhotE2x2;   //!
   TBranch        *b_PhotS4;   //!
   TBranch        *b_PhotE1x5Full5x5;   //!
   TBranch        *b_PhotE2x5Full5x5;   //!
   TBranch        *b_PhotE3x3Full5x5;   //!
   TBranch        *b_PhotE5x5Full5x5;   //!
   TBranch        *b_PhotSigmaIetaIetaFull5x5;   //!
   TBranch        *b_PhotR9Full5x5;   //!
   TBranch        *b_PhotId;   //!
   TBranch        *b_PhotHasPixelSeed;   //!
   TBranch        *b_PhotPassElVeto;   //!
   TBranch        *b_JetAk04Pt;   //!
   TBranch        *b_JetAk04Eta;   //!
   TBranch        *b_JetAk04Phi;   //!
   TBranch        *b_JetAk04E;   //!
   TBranch        *b_JetAk04Id;   //!
   TBranch        *b_JetAk04PuId;   //!
   TBranch        *b_JetAk04PuMva;   //!
   TBranch        *b_JetAk04RawPt;   //!
   TBranch        *b_JetAk04RawE;   //!
   TBranch        *b_JetAk04HfHadE;   //!
   TBranch        *b_JetAk04HfEmE;   //!
   TBranch        *b_JetAk04ChHadFrac;   //!
   TBranch        *b_JetAk04NeutralHadAndHfFrac;   //!
   TBranch        *b_JetAk04ChEmFrac;   //!
   TBranch        *b_JetAk04NeutralEmFrac;   //!
   TBranch        *b_JetAk04ChMult;   //!
   TBranch        *b_JetAk04ConstCnt;   //!
   TBranch        *b_JetAk04Beta;   //!
   TBranch        *b_JetAk04BetaClassic;   //!
   TBranch        *b_JetAk04BetaStar;   //!
   TBranch        *b_JetAk04BetaStarClassic;   //!
   TBranch        *b_JetAk04Rms;   //!
   TBranch        *b_JetAk04BTagCsv;   //!
   TBranch        *b_JetAk04BTagCsvV1;   //!
   TBranch        *b_JetAk04BTagCsvSLV1;   //!
   TBranch        *b_JetAk04BDiscCisvV2;   //!
   TBranch        *b_JetAk04BDiscJp;   //!
   TBranch        *b_JetAk04BDiscBjp;   //!
   TBranch        *b_JetAk04BDiscTche;   //!
   TBranch        *b_JetAk04BDiscTchp;   //!
   TBranch        *b_JetAk04BDiscSsvhe;   //!
   TBranch        *b_JetAk04BDiscSsvhp;   //!
   TBranch        *b_JetAk04PartFlav;   //!
   TBranch        *b_JetAk04HadFlav;   //!
   TBranch        *b_JetAk04JecUncUp;   //!
   TBranch        *b_JetAk04JecUncDwn;   //!
   TBranch        *b_JetAk04ConstId;   //!
   TBranch        *b_JetAk04ConstPt;   //!
   TBranch        *b_JetAk04ConstEta;   //!
   TBranch        *b_JetAk04ConstPhi;   //!
   TBranch        *b_JetAk04ConstE;   //!
   TBranch        *b_JetAk04GenJet;   //!
   TBranch        *b_JetAk08Pt;   //!
   TBranch        *b_JetAk08Eta;   //!
   TBranch        *b_JetAk08Phi;   //!
   TBranch        *b_JetAk08E;   //!
   TBranch        *b_JetAk08Id;   //!
   TBranch        *b_JetAk08RawPt;   //!
   TBranch        *b_JetAk08RawE;   //!
   TBranch        *b_JetAk08HfHadE;   //!
   TBranch        *b_JetAk08HfEmE;   //!
   TBranch        *b_JetAk08ChHadFrac;   //!
   TBranch        *b_JetAk08NeutralHadAndHfFrac;   //!
   TBranch        *b_JetAk08ChEmFrac;   //!
   TBranch        *b_JetAk08NeutralEmFrac;   //!
   TBranch        *b_JetAk08ChMult;   //!
   TBranch        *b_JetAk08ConstCnt;   //!
   TBranch        *b_JetAk08BTagCsv;   //!
   TBranch        *b_JetAk08BTagCsvV1;   //!
   TBranch        *b_JetAk08BTagCsvSLV1;   //!
   TBranch        *b_JetAk08BDiscCisvV2;   //!
   TBranch        *b_JetAk08BDiscJp;   //!
   TBranch        *b_JetAk08BDiscBjp;   //!
   TBranch        *b_JetAk08BDiscTche;   //!
   TBranch        *b_JetAk08BDiscTchp;   //!
   TBranch        *b_JetAk08BDiscSsvhe;   //!
   TBranch        *b_JetAk08BDiscSsvhp;   //!
   TBranch        *b_JetAk08PartFlav;   //!
   TBranch        *b_JetAk08HadFlav;   //!
   TBranch        *b_JetAk08JecUncUp;   //!
   TBranch        *b_JetAk08JecUncDwn;   //!
   TBranch        *b_JetAk08ConstId;   //!
   TBranch        *b_JetAk08ConstPt;   //!
   TBranch        *b_JetAk08ConstEta;   //!
   TBranch        *b_JetAk08ConstPhi;   //!
   TBranch        *b_JetAk08ConstE;   //!
   TBranch        *b_JetAk08GenJet;   //!
   TBranch        *b_JetAk08PrunedMass;   //!
   TBranch        *b_JetAk08FilteredMass;   //!
   TBranch        *b_JetAk08SoftDropMass;   //!
   TBranch        *b_JetAk08TrimmedMass;   //!
   TBranch        *b_JetAk08Tau1;   //!
   TBranch        *b_JetAk08Tau2;   //!
   TBranch        *b_JetAk08Tau3;   //!

   HZZ2l2nuLooper(TString, int, int, TString, int, int, float);
   virtual ~HZZ2l2nuLooper();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     FillNbEntries(TChain *);
   virtual void     FillTheTChain(TChain *, TString, int, int);
};

#ifdef HZZ2l2nuLooper_cxx
HZZ2l2nuLooper::HZZ2l2nuLooper(TString fileName, int skipFile, int maxFiles, TString outputFile, int maxEvents,int isMC, float crossSection) : fChain(0)
{

  outputFile_ = outputFile;
  maxEvents_ = maxEvents;
  isMC_ = isMC;
  sampleXsection_  = crossSection;

  totalEventsInBaobab_=-1;
  sumWeightInBaobab_=-1;
  sumWeightInBonzai_=-1;

  //First  get the tot number of events from the BonzaiHeader
  TChain * chainHeader = new TChain("tupel/BonzaiHeader","");
  FillTheTChain(chainHeader, fileName, skipFile, maxFiles);
  FillNbEntries(chainHeader);
  delete chainHeader;

  TChain * chain = new TChain("tupel/EventTree","");
  FillTheTChain(chain, fileName, skipFile, maxFiles);
  TTree *tree = chain;
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("tupel/EventTree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("tupel/EventTree","");
      chain->Add("mySkim.root/tupel/EventTree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}


HZZ2l2nuLooper::~HZZ2l2nuLooper()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void HZZ2l2nuLooper::FillTheTChain(TChain *theChain, TString theInputCatalog, int skipFiles, int maxFiles){
  cout << "catalog name=" << theInputCatalog << endl;

  std::ifstream f(theInputCatalog);
  if(!f.good()){
    std::cerr << "Failed to open file "<< theInputCatalog << "!\n";
    return;
  }


  int iline = 0;
  int nfiles = 0;
  std::string firstFile_ = "";
  while(f.good()){
    ++iline;
    std::string l;
    std::string::size_type p;

    std::getline(f, l);

    //trim white spaces:
    p = l.find_first_not_of(" \t");
    if(p!=std::string::npos) l.erase(0, p);
    p = l.find_last_not_of(" \t\n\r");
    if(p!=std::string::npos) l.erase(p + 1);
    else l.clear();

    //skip empty lines and comment lines:
    if (!l.size() || l[0] == '#') continue;

    if (!l.size() || l[0] == '*') continue;

    //extract first column (file name):
    p = l.find_first_of(" \t");
    if(p!=std::string::npos) l.erase(p);

    //sanity check:
    const char ext[6] = ".root";

    if(l.size() < sizeof(ext) || l.substr(l.size() - sizeof(ext) + 1) != ext){
      std::cerr << "Line " << iline << " of catalog file " << theInputCatalog << " was skipped.\n";
      continue;
    }



    if(skipFiles <= 0){
      ++nfiles;
      if((maxFiles > 0) &&  (nfiles > maxFiles)) break;
      std::cout << "Add file " << l.c_str() << " to the list of input files.\n";
      theChain->Add(l.c_str());
      if(firstFile_.size()==0) firstFile_ = l;
    } else{
    --skipFiles;
  }
}

return ;

}

Int_t HZZ2l2nuLooper::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HZZ2l2nuLooper::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HZZ2l2nuLooper::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   EvtWeights = 0;
   METPt = 0;
   METPx = 0;
   METPy = 0;
   METPz = 0;
   METE = 0;
   METsigx2 = 0;
   METsigxy = 0;
   METsigy2 = 0;
   METsig = 0;
   GMETPt = 0;
   GMETPx = 0;
   GMETPy = 0;
   GMETPz = 0;
   GMETE = 0;
   GLepDr01Pt = 0;
   GLepDr01Eta = 0;
   GLepDr01Phi = 0;
   GLepDr01E = 0;
   GLepDr01Id = 0;
   GLepDr01St = 0;
   GLepDr01MomId = 0;
   GLepDr01Prompt = 0;
   GLepDr01TauProd = 0;
   GLepBarePt = 0;
   GLepBareEta = 0;
   GLepBarePhi = 0;
   GLepBareE = 0;
   GLepBareId = 0;
   GLepBareSt = 0;
   GLepBareMomId = 0;
   GLepBarePrompt = 0;
   GLepBareTauProd = 0;
   GLepSt3Pt = 0;
   GLepSt3Eta = 0;
   GLepSt3Phi = 0;
   GLepSt3E = 0;
   GLepSt3Id = 0;
   GLepSt3St = 0;
   GLepSt3Mother0Id = 0;
   GLepSt3MotherCnt = 0;
   GPhotPt = 0;
   GPhotEta = 0;
   GPhotPhi = 0;
   GPhotE = 0;
   GPhotMotherId = 0;
   GPhotSt = 0;
   GPhotIsoEDR03 = 0;
   GPhotIsoEDR04 = 0;
   GPhotIsoEDR05 = 0;
   GPhotIsoSumPtDR03 = 0;
   GPhotIsoSumPtDR04 = 0;
   GPhotIsoSumPtDR05 = 0;
   GLepClosePhotPt = 0;
   GLepClosePhotEta = 0;
   GLepClosePhotPhi = 0;
   GLepClosePhotE = 0;
   GLepClosePhotId = 0;
   GLepClosePhotMother0Id = 0;
   GLepClosePhotMotherCnt = 0;
   GLepClosePhotSt = 0;
   GJetAk04Pt = 0;
   GJetAk04Eta = 0;
   GJetAk04Phi = 0;
   GJetAk04E = 0;
   GJetAk04ChFrac = 0;
   GJetAk04ConstCnt = 0;
   GJetAk04ConstId = 0;
   GJetAk04ConstPt = 0;
   GJetAk04ConstEta = 0;
   GJetAk04ConstPhi = 0;
   GJetAk04ConstE = 0;
   GJetAk04MatchedPartonID = 0;
   GJetAk04MatchedPartonDR = 0;
   GJetAk08Pt = 0;
   GJetAk08Eta = 0;
   GJetAk08Phi = 0;
   GJetAk08E = 0;
   GJetAk08ChFrac = 0;
   GJetAk08ConstCnt = 0;
   GJetAk08ConstId = 0;
   GJetAk08ConstPt = 0;
   GJetAk08ConstEta = 0;
   GJetAk08ConstPhi = 0;
   GJetAk08ConstE = 0;
   GJetAk08MatchedPartonID = 0;
   GJetAk08MatchedPartonDR = 0;
   GPdfId1 = 0;
   GPdfId2 = 0;
   GPdfx1 = 0;
   GPdfx2 = 0;
   GPdfScale = 0;
   MuPt = 0;
   MuEta = 0;
   MuPhi = 0;
   MuE = 0;
   MuId = 0;
   MuIdTight = 0;
   MuCh = 0;
   MuVtxZ = 0;
   MuDxy = 0;
   MuIsoRho = 0;
   MuPfIso = 0;
   MuType = 0;
   MuIsoTkIsoAbs = 0;
   MuIsoTkIsoRel = 0;
   MuIsoCalAbs = 0;
   MuIsoCombRel = 0;
   MuTkNormChi2 = 0;
   MuTkHitCnt = 0;
   MuMatchedStationCnt = 0;
   MuDz = 0;
   MuPixelHitCnt = 0;
   MuTkLayerCnt = 0;
   MuPfIsoChHad = 0;
   MuPfIsoNeutralHad = 0;
   MuPfIsoRawRel = 0;
   MuHltMatch = 0;
   ElPt = 0;
   ElEta = 0;
   ElEtaSc = 0;
   ElPhi = 0;
   ElE = 0;
   ElId = 0;
   ElCh = 0;
   ElScRawE = 0;
   ElCorrE = 0;
   ElEcalIso = 0;
   ElEcalPfIso = 0;
   ElMvaTrig = 0;
   ElMvaNonTrig = 0;
   ElMvaPresel = 0;
   ElDEtaTkScAtVtx = 0;
   ElDPhiTkScAtVtx = 0;
   ElHoE = 0;
   ElSigmaIetaIeta = 0;
   ElSigmaIetaIetaFull5x5 = 0;
   ElEinvMinusPinv = 0;
   ElD0 = 0;
   ElDz = 0;
   ElExpectedMissingInnerHitCnt = 0;
   ElPassConvVeto = 0;
   ElHltMatch = 0;
   ElPfIsoChHad = 0;
   ElPfIsoNeutralHad = 0;
   ElPfIsoIso = 0;
   ElPfIsoPuChHad = 0;
   ElPfIsoRaw = 0;
   ElPfIsoDbeta = 0;
   ElPfIsoRho = 0;
   ElAEff = 0;
   ElDr03TkSumPt = 0;
   ElDr03EcalRecHitSumEt = 0;
   ElDr03HcalTowerSumEt = 0;
   TauPt = 0;
   TauEta = 0;
   TauPhi = 0;
   TauE = 0;
   TauCh = 0;
   TauDecayModeFinding = 0;
   TauCombinedIsolationDeltaBetaCorrRaw3Hits = 0;
   TauDiscMuonLoose = 0;
   TauDiscMuonTight = 0;
   TauDiscElVLoose = 0;
   TauDiscElLoose = 0;
   TauDiscElVTight = 0;
   TauDiscElTight = 0;
   PhotPt = 0;
   PhotEta = 0;
   PhotPhi = 0;
   PhotScRawE = 0;
   PhotScEta = 0;
   PhotScPhi = 0;
   PhotIsoEcal = 0;
   PhotIsoHcal = 0;
   PhotIsoTk = 0;
   PhotPfIsoChHad = 0;
   PhotPfIsoNeutralHad = 0;
   PhotPfIsoPhot = 0;
   PhotPfIsoPuChHad = 0;
   PhotPfIsoEcalClus = 0;
   PhotPfIsoHcalClus = 0;
   PhotE3x3 = 0;
   PhotE1x5 = 0;
   PhotE2x5 = 0;
   PhotE5x5 = 0;
   PhotSigmaIetaIeta = 0;
   PhotHoE = 0;
   PhotEtaWidth = 0;
   PhotPhiWidth = 0;
   PhotR9 = 0;
   PhotE1x3 = 0;
   PhotE2x2 = 0;
   PhotS4 = 0;
   PhotE1x5Full5x5 = 0;
   PhotE2x5Full5x5 = 0;
   PhotE3x3Full5x5 = 0;
   PhotE5x5Full5x5 = 0;
   PhotSigmaIetaIetaFull5x5 = 0;
   PhotR9Full5x5 = 0;
   PhotId = 0;
   PhotHasPixelSeed = 0;
   PhotPassElVeto = 0;
   JetAk04Pt = 0;
   JetAk04Eta = 0;
   JetAk04Phi = 0;
   JetAk04E = 0;
   JetAk04Id = 0;
   JetAk04PuId = 0;
   JetAk04PuMva = 0;
   JetAk04RawPt = 0;
   JetAk04RawE = 0;
   JetAk04HfHadE = 0;
   JetAk04HfEmE = 0;
   JetAk04ChHadFrac = 0;
   JetAk04NeutralHadAndHfFrac = 0;
   JetAk04ChEmFrac = 0;
   JetAk04NeutralEmFrac = 0;
   JetAk04ChMult = 0;
   JetAk04ConstCnt = 0;
   JetAk04Beta = 0;
   JetAk04BetaClassic = 0;
   JetAk04BetaStar = 0;
   JetAk04BetaStarClassic = 0;
   JetAk04Rms = 0;
   JetAk04BTagCsv = 0;
   JetAk04BTagCsvV1 = 0;
   JetAk04BTagCsvSLV1 = 0;
   JetAk04BDiscCisvV2 = 0;
   JetAk04BDiscJp = 0;
   JetAk04BDiscBjp = 0;
   JetAk04BDiscTche = 0;
   JetAk04BDiscTchp = 0;
   JetAk04BDiscSsvhe = 0;
   JetAk04BDiscSsvhp = 0;
   JetAk04PartFlav = 0;
   JetAk04HadFlav = 0;
   JetAk04JecUncUp = 0;
   JetAk04JecUncDwn = 0;
   JetAk04ConstId = 0;
   JetAk04ConstPt = 0;
   JetAk04ConstEta = 0;
   JetAk04ConstPhi = 0;
   JetAk04ConstE = 0;
   JetAk04GenJet = 0;
   JetAk08Pt = 0;
   JetAk08Eta = 0;
   JetAk08Phi = 0;
   JetAk08E = 0;
   JetAk08Id = 0;
   JetAk08RawPt = 0;
   JetAk08RawE = 0;
   JetAk08HfHadE = 0;
   JetAk08HfEmE = 0;
   JetAk08ChHadFrac = 0;
   JetAk08NeutralHadAndHfFrac = 0;
   JetAk08ChEmFrac = 0;
   JetAk08NeutralEmFrac = 0;
   JetAk08ChMult = 0;
   JetAk08ConstCnt = 0;
   JetAk08BTagCsv = 0;
   JetAk08BTagCsvV1 = 0;
   JetAk08BTagCsvSLV1 = 0;
   JetAk08BDiscCisvV2 = 0;
   JetAk08BDiscJp = 0;
   JetAk08BDiscBjp = 0;
   JetAk08BDiscTche = 0;
   JetAk08BDiscTchp = 0;
   JetAk08BDiscSsvhe = 0;
   JetAk08BDiscSsvhp = 0;
   JetAk08PartFlav = 0;
   JetAk08HadFlav = 0;
   JetAk08JecUncUp = 0;
   JetAk08JecUncDwn = 0;
   JetAk08ConstId = 0;
   JetAk08ConstPt = 0;
   JetAk08ConstEta = 0;
   JetAk08ConstPhi = 0;
   JetAk08ConstE = 0;
   JetAk08GenJet = 0;
   JetAk08PrunedMass = 0;
   JetAk08FilteredMass = 0;
   JetAk08SoftDropMass = 0;
   JetAk08TrimmedMass = 0;
   JetAk08Tau1 = 0;
   JetAk08Tau2 = 0;
   JetAk08Tau3 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EvtIsRealData", &EvtIsRealData, &b_EvtIsRealData);
   fChain->SetBranchAddress("EvtNum", &EvtNum, &b_EvtNum);
   fChain->SetBranchAddress("EvtRunNum", &EvtRunNum, &b_EvtRunNum);
   fChain->SetBranchAddress("EvtLumiNum", &EvtLumiNum, &b_EvtLumiNum);
   fChain->SetBranchAddress("EvtBxNum", &EvtBxNum, &b_EvtBxNum);
   fChain->SetBranchAddress("EvtVtxCnt", &EvtVtxCnt, &b_EvtVtxCnt);
   fChain->SetBranchAddress("EvtPuCnt", &EvtPuCnt, &b_EvtPuCnt);
   fChain->SetBranchAddress("EvtPuCntTruth", &EvtPuCntTruth, &b_EvtPuCntTruth);
   fChain->SetBranchAddress("EvtWeights", &EvtWeights, &b_EvtWeights);
   fChain->SetBranchAddress("EvtFastJetRho", &EvtFastJetRho, &b_EvtFastJetRho);
   fChain->SetBranchAddress("TrigMET", &TrigMET, &b_TrigMET);
   fChain->SetBranchAddress("TrigHlt", &TrigHlt, &b_TrigHlt);
   fChain->SetBranchAddress("TrigHltPhot", &TrigHltPhot, &b_TrigHltPhot);
   fChain->SetBranchAddress("TrigHltDiPhot", &TrigHltDiPhot, &b_TrigHltDiPhot);
   fChain->SetBranchAddress("TrigHltMu", &TrigHltMu, &b_TrigHltMu);
   fChain->SetBranchAddress("TrigHltDiMu", &TrigHltDiMu, &b_TrigHltDiMu);
   fChain->SetBranchAddress("TrigHltEl", &TrigHltEl, &b_TrigHltEl);
   fChain->SetBranchAddress("TrigHltDiEl", &TrigHltDiEl, &b_TrigHltDiEl);
   fChain->SetBranchAddress("TrigHltElMu", &TrigHltElMu, &b_TrigHltElMu);
   fChain->SetBranchAddress("METPt", &METPt, &b_METPt);
   fChain->SetBranchAddress("METPx", &METPx, &b_METPx);
   fChain->SetBranchAddress("METPy", &METPy, &b_METPy);
   fChain->SetBranchAddress("METPz", &METPz, &b_METPz);
   fChain->SetBranchAddress("METE", &METE, &b_METE);
   fChain->SetBranchAddress("METsigx2", &METsigx2, &b_METsigx2);
   fChain->SetBranchAddress("METsigxy", &METsigxy, &b_METsigxy);
   fChain->SetBranchAddress("METsigy2", &METsigy2, &b_METsigy2);
   fChain->SetBranchAddress("METsig", &METsig, &b_METsig);
   fChain->SetBranchAddress("GMETPt", &GMETPt, &b_GMETPt);
   fChain->SetBranchAddress("GMETPx", &GMETPx, &b_GMETPx);
   fChain->SetBranchAddress("GMETPy", &GMETPy, &b_GMETPy);
   fChain->SetBranchAddress("GMETPz", &GMETPz, &b_GMETPz);
   fChain->SetBranchAddress("GMETE", &GMETE, &b_GMETE);
   fChain->SetBranchAddress("GLepDr01Pt", &GLepDr01Pt, &b_GLepDr01Pt);
   fChain->SetBranchAddress("GLepDr01Eta", &GLepDr01Eta, &b_GLepDr01Eta);
   fChain->SetBranchAddress("GLepDr01Phi", &GLepDr01Phi, &b_GLepDr01Phi);
   fChain->SetBranchAddress("GLepDr01E", &GLepDr01E, &b_GLepDr01E);
   fChain->SetBranchAddress("GLepDr01Id", &GLepDr01Id, &b_GLepDr01Id);
   fChain->SetBranchAddress("GLepDr01St", &GLepDr01St, &b_GLepDr01St);
   fChain->SetBranchAddress("GLepDr01MomId", &GLepDr01MomId, &b_GLepDr01MomId);
   fChain->SetBranchAddress("GLepDr01Prompt", &GLepDr01Prompt, &b_GLepDr01Prompt);
   fChain->SetBranchAddress("GLepDr01TauProd", &GLepDr01TauProd, &b_GLepDr01TauProd);
   fChain->SetBranchAddress("GLepBarePt", &GLepBarePt, &b_GLepBarePt);
   fChain->SetBranchAddress("GLepBareEta", &GLepBareEta, &b_GLepBareEta);
   fChain->SetBranchAddress("GLepBarePhi", &GLepBarePhi, &b_GLepBarePhi);
   fChain->SetBranchAddress("GLepBareE", &GLepBareE, &b_GLepBareE);
   fChain->SetBranchAddress("GLepBareId", &GLepBareId, &b_GLepBareId);
   fChain->SetBranchAddress("GLepBareSt", &GLepBareSt, &b_GLepBareSt);
   fChain->SetBranchAddress("GLepBareMomId", &GLepBareMomId, &b_GLepBareMomId);
   fChain->SetBranchAddress("GLepBarePrompt", &GLepBarePrompt, &b_GLepBarePrompt);
   fChain->SetBranchAddress("GLepBareTauProd", &GLepBareTauProd, &b_GLepBareTauProd);
   fChain->SetBranchAddress("GLepSt3Pt", &GLepSt3Pt, &b_GLepSt3Pt);
   fChain->SetBranchAddress("GLepSt3Eta", &GLepSt3Eta, &b_GLepSt3Eta);
   fChain->SetBranchAddress("GLepSt3Phi", &GLepSt3Phi, &b_GLepSt3Phi);
   fChain->SetBranchAddress("GLepSt3E", &GLepSt3E, &b_GLepSt3E);
   fChain->SetBranchAddress("GLepSt3Id", &GLepSt3Id, &b_GLepSt3Id);
   fChain->SetBranchAddress("GLepSt3St", &GLepSt3St, &b_GLepSt3St);
   fChain->SetBranchAddress("GLepSt3Mother0Id", &GLepSt3Mother0Id, &b_GLepSt3Mother0Id);
   fChain->SetBranchAddress("GLepSt3MotherCnt", &GLepSt3MotherCnt, &b_GLepSt3MotherCnt);
   fChain->SetBranchAddress("GPhotPt", &GPhotPt, &b_GPhotPt);
   fChain->SetBranchAddress("GPhotEta", &GPhotEta, &b_GPhotEta);
   fChain->SetBranchAddress("GPhotPhi", &GPhotPhi, &b_GPhotPhi);
   fChain->SetBranchAddress("GPhotE", &GPhotE, &b_GPhotE);
   fChain->SetBranchAddress("GPhotMotherId", &GPhotMotherId, &b_GPhotMotherId);
   fChain->SetBranchAddress("GPhotSt", &GPhotSt, &b_GPhotSt);
   fChain->SetBranchAddress("GPhotIsoEDR03", &GPhotIsoEDR03, &b_GPhotIsoEDR03);
   fChain->SetBranchAddress("GPhotIsoEDR04", &GPhotIsoEDR04, &b_GPhotIsoEDR04);
   fChain->SetBranchAddress("GPhotIsoEDR05", &GPhotIsoEDR05, &b_GPhotIsoEDR05);
   fChain->SetBranchAddress("GPhotIsoSumPtDR03", &GPhotIsoSumPtDR03, &b_GPhotIsoSumPtDR03);
   fChain->SetBranchAddress("GPhotIsoSumPtDR04", &GPhotIsoSumPtDR04, &b_GPhotIsoSumPtDR04);
   fChain->SetBranchAddress("GPhotIsoSumPtDR05", &GPhotIsoSumPtDR05, &b_GPhotIsoSumPtDR05);
   fChain->SetBranchAddress("GLepClosePhotPt", &GLepClosePhotPt, &b_GLepClosePhotPt);
   fChain->SetBranchAddress("GLepClosePhotEta", &GLepClosePhotEta, &b_GLepClosePhotEta);
   fChain->SetBranchAddress("GLepClosePhotPhi", &GLepClosePhotPhi, &b_GLepClosePhotPhi);
   fChain->SetBranchAddress("GLepClosePhotE", &GLepClosePhotE, &b_GLepClosePhotE);
   fChain->SetBranchAddress("GLepClosePhotId", &GLepClosePhotId, &b_GLepClosePhotId);
   fChain->SetBranchAddress("GLepClosePhotMother0Id", &GLepClosePhotMother0Id, &b_GLepClosePhotMother0Id);
   fChain->SetBranchAddress("GLepClosePhotMotherCnt", &GLepClosePhotMotherCnt, &b_GLepClosePhotMotherCnt);
   fChain->SetBranchAddress("GLepClosePhotSt", &GLepClosePhotSt, &b_GLepClosePhotSt);
   fChain->SetBranchAddress("GJetAk04Pt", &GJetAk04Pt, &b_GJetAk04Pt);
   fChain->SetBranchAddress("GJetAk04Eta", &GJetAk04Eta, &b_GJetAk04Eta);
   fChain->SetBranchAddress("GJetAk04Phi", &GJetAk04Phi, &b_GJetAk04Phi);
   fChain->SetBranchAddress("GJetAk04E", &GJetAk04E, &b_GJetAk04E);
   fChain->SetBranchAddress("GJetAk04ChFrac", &GJetAk04ChFrac, &b_GJetAk04ChFrac);
   fChain->SetBranchAddress("GJetAk04ConstCnt", &GJetAk04ConstCnt, &b_GJetAk04ConstCnt);
   fChain->SetBranchAddress("GJetAk04ConstId", &GJetAk04ConstId, &b_GJetAk04ConstId);
   fChain->SetBranchAddress("GJetAk04ConstPt", &GJetAk04ConstPt, &b_GJetAk04ConstPt);
   fChain->SetBranchAddress("GJetAk04ConstEta", &GJetAk04ConstEta, &b_GJetAk04ConstEta);
   fChain->SetBranchAddress("GJetAk04ConstPhi", &GJetAk04ConstPhi, &b_GJetAk04ConstPhi);
   fChain->SetBranchAddress("GJetAk04ConstE", &GJetAk04ConstE, &b_GJetAk04ConstE);
   fChain->SetBranchAddress("GJetAk04MatchedPartonID", &GJetAk04MatchedPartonID, &b_GJetAk04MatchedPartonID);
   fChain->SetBranchAddress("GJetAk04MatchedPartonDR", &GJetAk04MatchedPartonDR, &b_GJetAk04MatchedPartonDR);
   fChain->SetBranchAddress("GJetAk08Pt", &GJetAk08Pt, &b_GJetAk08Pt);
   fChain->SetBranchAddress("GJetAk08Eta", &GJetAk08Eta, &b_GJetAk08Eta);
   fChain->SetBranchAddress("GJetAk08Phi", &GJetAk08Phi, &b_GJetAk08Phi);
   fChain->SetBranchAddress("GJetAk08E", &GJetAk08E, &b_GJetAk08E);
   fChain->SetBranchAddress("GJetAk08ChFrac", &GJetAk08ChFrac, &b_GJetAk08ChFrac);
   fChain->SetBranchAddress("GJetAk08ConstCnt", &GJetAk08ConstCnt, &b_GJetAk08ConstCnt);
   fChain->SetBranchAddress("GJetAk08ConstId", &GJetAk08ConstId, &b_GJetAk08ConstId);
   fChain->SetBranchAddress("GJetAk08ConstPt", &GJetAk08ConstPt, &b_GJetAk08ConstPt);
   fChain->SetBranchAddress("GJetAk08ConstEta", &GJetAk08ConstEta, &b_GJetAk08ConstEta);
   fChain->SetBranchAddress("GJetAk08ConstPhi", &GJetAk08ConstPhi, &b_GJetAk08ConstPhi);
   fChain->SetBranchAddress("GJetAk08ConstE", &GJetAk08ConstE, &b_GJetAk08ConstE);
   fChain->SetBranchAddress("GJetAk08MatchedPartonID", &GJetAk08MatchedPartonID, &b_GJetAk08MatchedPartonID);
   fChain->SetBranchAddress("GJetAk08MatchedPartonDR", &GJetAk08MatchedPartonDR, &b_GJetAk08MatchedPartonDR);
   fChain->SetBranchAddress("GPdfId1", &GPdfId1, &b_GPdfId1);
   fChain->SetBranchAddress("GPdfId2", &GPdfId2, &b_GPdfId2);
   fChain->SetBranchAddress("GPdfx1", &GPdfx1, &b_GPdfx1);
   fChain->SetBranchAddress("GPdfx2", &GPdfx2, &b_GPdfx2);
   fChain->SetBranchAddress("GPdfScale", &GPdfScale, &b_GPdfScale);
   fChain->SetBranchAddress("GBinningValue", &GBinningValue, &b_GBinningValue);
   fChain->SetBranchAddress("GNup", &GNup, &b_GNup);
   fChain->SetBranchAddress("MuPt", &MuPt, &b_MuPt);
   fChain->SetBranchAddress("MuEta", &MuEta, &b_MuEta);
   fChain->SetBranchAddress("MuPhi", &MuPhi, &b_MuPhi);
   fChain->SetBranchAddress("MuE", &MuE, &b_MuE);
   fChain->SetBranchAddress("MuId", &MuId, &b_MuId);
   fChain->SetBranchAddress("MuIdTight", &MuIdTight, &b_MuIdTight);
   fChain->SetBranchAddress("MuCh", &MuCh, &b_MuCh);
   fChain->SetBranchAddress("MuVtxZ", &MuVtxZ, &b_MuVtxZ);
   fChain->SetBranchAddress("MuDxy", &MuDxy, &b_MuDxy);
   fChain->SetBranchAddress("MuIsoRho", &MuIsoRho, &b_MuIsoRho);
   fChain->SetBranchAddress("MuPfIso", &MuPfIso, &b_MuPfIso);
   fChain->SetBranchAddress("MuType", &MuType, &b_MuType);
   fChain->SetBranchAddress("MuIsoTkIsoAbs", &MuIsoTkIsoAbs, &b_MuIsoTkIsoAbs);
   fChain->SetBranchAddress("MuIsoTkIsoRel", &MuIsoTkIsoRel, &b_MuIsoTkIsoRel);
   fChain->SetBranchAddress("MuIsoCalAbs", &MuIsoCalAbs, &b_MuIsoCalAbs);
   fChain->SetBranchAddress("MuIsoCombRel", &MuIsoCombRel, &b_MuIsoCombRel);
   fChain->SetBranchAddress("MuTkNormChi2", &MuTkNormChi2, &b_MuTkNormChi2);
   fChain->SetBranchAddress("MuTkHitCnt", &MuTkHitCnt, &b_MuTkHitCnt);
   fChain->SetBranchAddress("MuMatchedStationCnt", &MuMatchedStationCnt, &b_MuMatchedStationCnt);
   fChain->SetBranchAddress("MuDz", &MuDz, &b_MuDz);
   fChain->SetBranchAddress("MuPixelHitCnt", &MuPixelHitCnt, &b_MuPixelHitCnt);
   fChain->SetBranchAddress("MuTkLayerCnt", &MuTkLayerCnt, &b_MuTkLayerCnt);
   fChain->SetBranchAddress("MuPfIsoChHad", &MuPfIsoChHad, &b_MuPfIsoChHad);
   fChain->SetBranchAddress("MuPfIsoNeutralHad", &MuPfIsoNeutralHad, &b_MuPfIsoNeutralHad);
   fChain->SetBranchAddress("MuPfIsoRawRel", &MuPfIsoRawRel, &b_MuPfIsoRawRel);
   fChain->SetBranchAddress("MuHltMatch", &MuHltMatch, &b_MuHltMatch);
   fChain->SetBranchAddress("ElPt", &ElPt, &b_ElPt);
   fChain->SetBranchAddress("ElEta", &ElEta, &b_ElEta);
   fChain->SetBranchAddress("ElEtaSc", &ElEtaSc, &b_ElEtaSc);
   fChain->SetBranchAddress("ElPhi", &ElPhi, &b_ElPhi);
   fChain->SetBranchAddress("ElE", &ElE, &b_ElE);
   fChain->SetBranchAddress("ElId", &ElId, &b_ElId);
   fChain->SetBranchAddress("ElCh", &ElCh, &b_ElCh);
   fChain->SetBranchAddress("ElScRawE", &ElScRawE, &b_ElScRawE);
   fChain->SetBranchAddress("ElCorrE", &ElCorrE, &b_ElCorrE);
   fChain->SetBranchAddress("ElEcalIso", &ElEcalIso, &b_ElEcalIso);
   fChain->SetBranchAddress("ElEcalPfIso", &ElEcalPfIso, &b_ElEcalPfIso);
   fChain->SetBranchAddress("ElMvaTrig", &ElMvaTrig, &b_ElMvaTrig);
   fChain->SetBranchAddress("ElMvaNonTrig", &ElMvaNonTrig, &b_ElMvaNonTrig);
   fChain->SetBranchAddress("ElMvaPresel", &ElMvaPresel, &b_ElMvaPresel);
   fChain->SetBranchAddress("ElDEtaTkScAtVtx", &ElDEtaTkScAtVtx, &b_ElDEtaTkScAtVtx);
   fChain->SetBranchAddress("ElDPhiTkScAtVtx", &ElDPhiTkScAtVtx, &b_ElDPhiTkScAtVtx);
   fChain->SetBranchAddress("ElHoE", &ElHoE, &b_ElHoE);
   fChain->SetBranchAddress("ElSigmaIetaIeta", &ElSigmaIetaIeta, &b_ElSigmaIetaIeta);
   fChain->SetBranchAddress("ElSigmaIetaIetaFull5x5", &ElSigmaIetaIetaFull5x5, &b_ElSigmaIetaIetaFull5x5);
   fChain->SetBranchAddress("ElEinvMinusPinv", &ElEinvMinusPinv, &b_ElEinvMinusPinv);
   fChain->SetBranchAddress("ElD0", &ElD0, &b_ElD0);
   fChain->SetBranchAddress("ElDz", &ElDz, &b_ElDz);
   fChain->SetBranchAddress("ElExpectedMissingInnerHitCnt", &ElExpectedMissingInnerHitCnt, &b_ElExpectedMissingInnerHitCnt);
   fChain->SetBranchAddress("ElPassConvVeto", &ElPassConvVeto, &b_ElPassConvVeto);
   fChain->SetBranchAddress("ElHltMatch", &ElHltMatch, &b_ElHltMatch);
   fChain->SetBranchAddress("ElPfIsoChHad", &ElPfIsoChHad, &b_ElPfIsoChHad);
   fChain->SetBranchAddress("ElPfIsoNeutralHad", &ElPfIsoNeutralHad, &b_ElPfIsoNeutralHad);
   fChain->SetBranchAddress("ElPfIsoIso", &ElPfIsoIso, &b_ElPfIsoIso);
   fChain->SetBranchAddress("ElPfIsoPuChHad", &ElPfIsoPuChHad, &b_ElPfIsoPuChHad);
   fChain->SetBranchAddress("ElPfIsoRaw", &ElPfIsoRaw, &b_ElPfIsoRaw);
   fChain->SetBranchAddress("ElPfIsoDbeta", &ElPfIsoDbeta, &b_ElPfIsoDbeta);
   fChain->SetBranchAddress("ElPfIsoRho", &ElPfIsoRho, &b_ElPfIsoRho);
   fChain->SetBranchAddress("ElAEff", &ElAEff, &b_ElAEff);
   fChain->SetBranchAddress("ElDr03TkSumPt", &ElDr03TkSumPt, &b_ElDr03TkSumPt);
   fChain->SetBranchAddress("ElDr03EcalRecHitSumEt", &ElDr03EcalRecHitSumEt, &b_ElDr03EcalRecHitSumEt);
   fChain->SetBranchAddress("ElDr03HcalTowerSumEt", &ElDr03HcalTowerSumEt, &b_ElDr03HcalTowerSumEt);
   fChain->SetBranchAddress("TauPt", &TauPt, &b_TauPt);
   fChain->SetBranchAddress("TauEta", &TauEta, &b_TauEta);
   fChain->SetBranchAddress("TauPhi", &TauPhi, &b_TauPhi);
   fChain->SetBranchAddress("TauE", &TauE, &b_TauE);
   fChain->SetBranchAddress("TauCh", &TauCh, &b_TauCh);
   fChain->SetBranchAddress("TauDecayModeFinding", &TauDecayModeFinding, &b_TauDecayModeFinding);
   fChain->SetBranchAddress("TauCombinedIsolationDeltaBetaCorrRaw3Hits", &TauCombinedIsolationDeltaBetaCorrRaw3Hits, &b_TauCombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("TauDiscMuonLoose", &TauDiscMuonLoose, &b_TauDiscMuonLoose);
   fChain->SetBranchAddress("TauDiscMuonTight", &TauDiscMuonTight, &b_TauDiscMuonTight);
   fChain->SetBranchAddress("TauDiscElVLoose", &TauDiscElVLoose, &b_TauDiscElVLoose);
   fChain->SetBranchAddress("TauDiscElLoose", &TauDiscElLoose, &b_TauDiscElLoose);
   fChain->SetBranchAddress("TauDiscElVTight", &TauDiscElVTight, &b_TauDiscElVTight);
   fChain->SetBranchAddress("TauDiscElTight", &TauDiscElTight, &b_TauDiscElTight);
   fChain->SetBranchAddress("PhotPt", &PhotPt, &b_PhotPt);
   fChain->SetBranchAddress("PhotEta", &PhotEta, &b_PhotEta);
   fChain->SetBranchAddress("PhotPhi", &PhotPhi, &b_PhotPhi);
   fChain->SetBranchAddress("PhotScRawE", &PhotScRawE, &b_PhotScRawE);
   fChain->SetBranchAddress("PhotScEta", &PhotScEta, &b_PhotScEta);
   fChain->SetBranchAddress("PhotScPhi", &PhotScPhi, &b_PhotScPhi);
   fChain->SetBranchAddress("PhotIsoEcal", &PhotIsoEcal, &b_PhotIsoEcal);
   fChain->SetBranchAddress("PhotIsoHcal", &PhotIsoHcal, &b_PhotIsoHcal);
   fChain->SetBranchAddress("PhotIsoTk", &PhotIsoTk, &b_PhotIsoTk);
   fChain->SetBranchAddress("PhotPfIsoChHad", &PhotPfIsoChHad, &b_PhotPfIsoChHad);
   fChain->SetBranchAddress("PhotPfIsoNeutralHad", &PhotPfIsoNeutralHad, &b_PhotPfIsoNeutralHad);
   fChain->SetBranchAddress("PhotPfIsoPhot", &PhotPfIsoPhot, &b_PhotPfIsoPhot);
   fChain->SetBranchAddress("PhotPfIsoPuChHad", &PhotPfIsoPuChHad, &b_PhotPfIsoPuChHad);
   fChain->SetBranchAddress("PhotPfIsoEcalClus", &PhotPfIsoEcalClus, &b_PhotPfIsoEcalClus);
   fChain->SetBranchAddress("PhotPfIsoHcalClus", &PhotPfIsoHcalClus, &b_PhotPfIsoHcalClus);
   fChain->SetBranchAddress("PhotE3x3", &PhotE3x3, &b_PhotE3x3);
   fChain->SetBranchAddress("PhotE1x5", &PhotE1x5, &b_PhotE1x5);
   fChain->SetBranchAddress("PhotE2x5", &PhotE2x5, &b_PhotE2x5);
   fChain->SetBranchAddress("PhotE5x5", &PhotE5x5, &b_PhotE5x5);
   fChain->SetBranchAddress("PhotSigmaIetaIeta", &PhotSigmaIetaIeta, &b_PhotSigmaIetaIeta);
   fChain->SetBranchAddress("PhotHoE", &PhotHoE, &b_PhotHoE);
   fChain->SetBranchAddress("PhotEtaWidth", &PhotEtaWidth, &b_PhotEtaWidth);
   fChain->SetBranchAddress("PhotPhiWidth", &PhotPhiWidth, &b_PhotPhiWidth);
   fChain->SetBranchAddress("PhotR9", &PhotR9, &b_PhotR9);
   fChain->SetBranchAddress("PhotE1x3", &PhotE1x3, &b_PhotE1x3);
   fChain->SetBranchAddress("PhotE2x2", &PhotE2x2, &b_PhotE2x2);
   fChain->SetBranchAddress("PhotS4", &PhotS4, &b_PhotS4);
   fChain->SetBranchAddress("PhotE1x5Full5x5", &PhotE1x5Full5x5, &b_PhotE1x5Full5x5);
   fChain->SetBranchAddress("PhotE2x5Full5x5", &PhotE2x5Full5x5, &b_PhotE2x5Full5x5);
   fChain->SetBranchAddress("PhotE3x3Full5x5", &PhotE3x3Full5x5, &b_PhotE3x3Full5x5);
   fChain->SetBranchAddress("PhotE5x5Full5x5", &PhotE5x5Full5x5, &b_PhotE5x5Full5x5);
   fChain->SetBranchAddress("PhotSigmaIetaIetaFull5x5", &PhotSigmaIetaIetaFull5x5, &b_PhotSigmaIetaIetaFull5x5);
   fChain->SetBranchAddress("PhotR9Full5x5", &PhotR9Full5x5, &b_PhotR9Full5x5);
   fChain->SetBranchAddress("PhotId", &PhotId, &b_PhotId);
   fChain->SetBranchAddress("PhotHasPixelSeed", &PhotHasPixelSeed, &b_PhotHasPixelSeed);
   fChain->SetBranchAddress("PhotPassElVeto", &PhotPassElVeto, &b_PhotPassElVeto);
   fChain->SetBranchAddress("JetAk04Pt", &JetAk04Pt, &b_JetAk04Pt);
   fChain->SetBranchAddress("JetAk04Eta", &JetAk04Eta, &b_JetAk04Eta);
   fChain->SetBranchAddress("JetAk04Phi", &JetAk04Phi, &b_JetAk04Phi);
   fChain->SetBranchAddress("JetAk04E", &JetAk04E, &b_JetAk04E);
   fChain->SetBranchAddress("JetAk04Id", &JetAk04Id, &b_JetAk04Id);
   fChain->SetBranchAddress("JetAk04PuId", &JetAk04PuId, &b_JetAk04PuId);
   fChain->SetBranchAddress("JetAk04PuMva", &JetAk04PuMva, &b_JetAk04PuMva);
   fChain->SetBranchAddress("JetAk04RawPt", &JetAk04RawPt, &b_JetAk04RawPt);
   fChain->SetBranchAddress("JetAk04RawE", &JetAk04RawE, &b_JetAk04RawE);
   fChain->SetBranchAddress("JetAk04HfHadE", &JetAk04HfHadE, &b_JetAk04HfHadE);
   fChain->SetBranchAddress("JetAk04HfEmE", &JetAk04HfEmE, &b_JetAk04HfEmE);
   fChain->SetBranchAddress("JetAk04ChHadFrac", &JetAk04ChHadFrac, &b_JetAk04ChHadFrac);
   fChain->SetBranchAddress("JetAk04NeutralHadAndHfFrac", &JetAk04NeutralHadAndHfFrac, &b_JetAk04NeutralHadAndHfFrac);
   fChain->SetBranchAddress("JetAk04ChEmFrac", &JetAk04ChEmFrac, &b_JetAk04ChEmFrac);
   fChain->SetBranchAddress("JetAk04NeutralEmFrac", &JetAk04NeutralEmFrac, &b_JetAk04NeutralEmFrac);
   fChain->SetBranchAddress("JetAk04ChMult", &JetAk04ChMult, &b_JetAk04ChMult);
   fChain->SetBranchAddress("JetAk04ConstCnt", &JetAk04ConstCnt, &b_JetAk04ConstCnt);
   fChain->SetBranchAddress("JetAk04Beta", &JetAk04Beta, &b_JetAk04Beta);
   fChain->SetBranchAddress("JetAk04BetaClassic", &JetAk04BetaClassic, &b_JetAk04BetaClassic);
   fChain->SetBranchAddress("JetAk04BetaStar", &JetAk04BetaStar, &b_JetAk04BetaStar);
   fChain->SetBranchAddress("JetAk04BetaStarClassic", &JetAk04BetaStarClassic, &b_JetAk04BetaStarClassic);
   fChain->SetBranchAddress("JetAk04Rms", &JetAk04Rms, &b_JetAk04Rms);
   fChain->SetBranchAddress("JetAk04BTagCsv", &JetAk04BTagCsv, &b_JetAk04BTagCsv);
   fChain->SetBranchAddress("JetAk04BTagCsvV1", &JetAk04BTagCsvV1, &b_JetAk04BTagCsvV1);
   fChain->SetBranchAddress("JetAk04BTagCsvSLV1", &JetAk04BTagCsvSLV1, &b_JetAk04BTagCsvSLV1);
   fChain->SetBranchAddress("JetAk04BDiscCisvV2", &JetAk04BDiscCisvV2, &b_JetAk04BDiscCisvV2);
   fChain->SetBranchAddress("JetAk04BDiscJp", &JetAk04BDiscJp, &b_JetAk04BDiscJp);
   fChain->SetBranchAddress("JetAk04BDiscBjp", &JetAk04BDiscBjp, &b_JetAk04BDiscBjp);
   fChain->SetBranchAddress("JetAk04BDiscTche", &JetAk04BDiscTche, &b_JetAk04BDiscTche);
   fChain->SetBranchAddress("JetAk04BDiscTchp", &JetAk04BDiscTchp, &b_JetAk04BDiscTchp);
   fChain->SetBranchAddress("JetAk04BDiscSsvhe", &JetAk04BDiscSsvhe, &b_JetAk04BDiscSsvhe);
   fChain->SetBranchAddress("JetAk04BDiscSsvhp", &JetAk04BDiscSsvhp, &b_JetAk04BDiscSsvhp);
   fChain->SetBranchAddress("JetAk04PartFlav", &JetAk04PartFlav, &b_JetAk04PartFlav);
   fChain->SetBranchAddress("JetAk04HadFlav", &JetAk04HadFlav, &b_JetAk04HadFlav);
   fChain->SetBranchAddress("JetAk04JecUncUp", &JetAk04JecUncUp, &b_JetAk04JecUncUp);
   fChain->SetBranchAddress("JetAk04JecUncDwn", &JetAk04JecUncDwn, &b_JetAk04JecUncDwn);
   fChain->SetBranchAddress("JetAk04ConstId", &JetAk04ConstId, &b_JetAk04ConstId);
   fChain->SetBranchAddress("JetAk04ConstPt", &JetAk04ConstPt, &b_JetAk04ConstPt);
   fChain->SetBranchAddress("JetAk04ConstEta", &JetAk04ConstEta, &b_JetAk04ConstEta);
   fChain->SetBranchAddress("JetAk04ConstPhi", &JetAk04ConstPhi, &b_JetAk04ConstPhi);
   fChain->SetBranchAddress("JetAk04ConstE", &JetAk04ConstE, &b_JetAk04ConstE);
   fChain->SetBranchAddress("JetAk04GenJet", &JetAk04GenJet, &b_JetAk04GenJet);
   fChain->SetBranchAddress("JetAk08Pt", &JetAk08Pt, &b_JetAk08Pt);
   fChain->SetBranchAddress("JetAk08Eta", &JetAk08Eta, &b_JetAk08Eta);
   fChain->SetBranchAddress("JetAk08Phi", &JetAk08Phi, &b_JetAk08Phi);
   fChain->SetBranchAddress("JetAk08E", &JetAk08E, &b_JetAk08E);
   fChain->SetBranchAddress("JetAk08Id", &JetAk08Id, &b_JetAk08Id);
   fChain->SetBranchAddress("JetAk08RawPt", &JetAk08RawPt, &b_JetAk08RawPt);
   fChain->SetBranchAddress("JetAk08RawE", &JetAk08RawE, &b_JetAk08RawE);
   fChain->SetBranchAddress("JetAk08HfHadE", &JetAk08HfHadE, &b_JetAk08HfHadE);
   fChain->SetBranchAddress("JetAk08HfEmE", &JetAk08HfEmE, &b_JetAk08HfEmE);
   fChain->SetBranchAddress("JetAk08ChHadFrac", &JetAk08ChHadFrac, &b_JetAk08ChHadFrac);
   fChain->SetBranchAddress("JetAk08NeutralHadAndHfFrac", &JetAk08NeutralHadAndHfFrac, &b_JetAk08NeutralHadAndHfFrac);
   fChain->SetBranchAddress("JetAk08ChEmFrac", &JetAk08ChEmFrac, &b_JetAk08ChEmFrac);
   fChain->SetBranchAddress("JetAk08NeutralEmFrac", &JetAk08NeutralEmFrac, &b_JetAk08NeutralEmFrac);
   fChain->SetBranchAddress("JetAk08ChMult", &JetAk08ChMult, &b_JetAk08ChMult);
   fChain->SetBranchAddress("JetAk08ConstCnt", &JetAk08ConstCnt, &b_JetAk08ConstCnt);
   fChain->SetBranchAddress("JetAk08BTagCsv", &JetAk08BTagCsv, &b_JetAk08BTagCsv);
   fChain->SetBranchAddress("JetAk08BTagCsvV1", &JetAk08BTagCsvV1, &b_JetAk08BTagCsvV1);
   fChain->SetBranchAddress("JetAk08BTagCsvSLV1", &JetAk08BTagCsvSLV1, &b_JetAk08BTagCsvSLV1);
   fChain->SetBranchAddress("JetAk08BDiscCisvV2", &JetAk08BDiscCisvV2, &b_JetAk08BDiscCisvV2);
   fChain->SetBranchAddress("JetAk08BDiscJp", &JetAk08BDiscJp, &b_JetAk08BDiscJp);
   fChain->SetBranchAddress("JetAk08BDiscBjp", &JetAk08BDiscBjp, &b_JetAk08BDiscBjp);
   fChain->SetBranchAddress("JetAk08BDiscTche", &JetAk08BDiscTche, &b_JetAk08BDiscTche);
   fChain->SetBranchAddress("JetAk08BDiscTchp", &JetAk08BDiscTchp, &b_JetAk08BDiscTchp);
   fChain->SetBranchAddress("JetAk08BDiscSsvhe", &JetAk08BDiscSsvhe, &b_JetAk08BDiscSsvhe);
   fChain->SetBranchAddress("JetAk08BDiscSsvhp", &JetAk08BDiscSsvhp, &b_JetAk08BDiscSsvhp);
   fChain->SetBranchAddress("JetAk08PartFlav", &JetAk08PartFlav, &b_JetAk08PartFlav);
   fChain->SetBranchAddress("JetAk08HadFlav", &JetAk08HadFlav, &b_JetAk08HadFlav);
   fChain->SetBranchAddress("JetAk08JecUncUp", &JetAk08JecUncUp, &b_JetAk08JecUncUp);
   fChain->SetBranchAddress("JetAk08JecUncDwn", &JetAk08JecUncDwn, &b_JetAk08JecUncDwn);
   fChain->SetBranchAddress("JetAk08ConstId", &JetAk08ConstId, &b_JetAk08ConstId);
   fChain->SetBranchAddress("JetAk08ConstPt", &JetAk08ConstPt, &b_JetAk08ConstPt);
   fChain->SetBranchAddress("JetAk08ConstEta", &JetAk08ConstEta, &b_JetAk08ConstEta);
   fChain->SetBranchAddress("JetAk08ConstPhi", &JetAk08ConstPhi, &b_JetAk08ConstPhi);
   fChain->SetBranchAddress("JetAk08ConstE", &JetAk08ConstE, &b_JetAk08ConstE);
   fChain->SetBranchAddress("JetAk08GenJet", &JetAk08GenJet, &b_JetAk08GenJet);
   fChain->SetBranchAddress("JetAk08PrunedMass", &JetAk08PrunedMass, &b_JetAk08PrunedMass);
   fChain->SetBranchAddress("JetAk08FilteredMass", &JetAk08FilteredMass, &b_JetAk08FilteredMass);
   fChain->SetBranchAddress("JetAk08SoftDropMass", &JetAk08SoftDropMass, &b_JetAk08SoftDropMass);
   fChain->SetBranchAddress("JetAk08TrimmedMass", &JetAk08TrimmedMass, &b_JetAk08TrimmedMass);
   fChain->SetBranchAddress("JetAk08Tau1", &JetAk08Tau1, &b_JetAk08Tau1);
   fChain->SetBranchAddress("JetAk08Tau2", &JetAk08Tau2, &b_JetAk08Tau2);
   fChain->SetBranchAddress("JetAk08Tau3", &JetAk08Tau3, &b_JetAk08Tau3);
   Notify();
}

Bool_t HZZ2l2nuLooper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HZZ2l2nuLooper::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HZZ2l2nuLooper::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void HZZ2l2nuLooper::FillNbEntries(TChain  *inputChain)
{

  TTree *treeBonzaiHeader = inputChain;
  int   InEvtCount; InEvtCount=0;
  vector<double>   *InEvtWeightSums; InEvtWeightSums=0;
  vector<double>   *EvtWeightSums; EvtWeightSums=0;

  TBranch *b_InEvtCount;
  TBranch *b_InEvtWeightSums;
  TBranch *b_EvtWeightSums;

  treeBonzaiHeader->SetBranchAddress("InEvtCount", &InEvtCount, &b_InEvtCount);
  treeBonzaiHeader->SetBranchAddress("InEvtWeightSums", &InEvtWeightSums, &b_InEvtWeightSums);
  treeBonzaiHeader->SetBranchAddress("EvtWeightSums", &EvtWeightSums, &b_EvtWeightSums);


  int nbEntriesInHeader = treeBonzaiHeader->GetEntries();
  if (nbEntriesInHeader<1) {
    cout << "ALERT: Nb of entries in bonzai header different smaller from 1 ! " << endl;
    return;
  }
  else{
    totalEventsInBaobab_ = 0;
    sumWeightInBaobab_ = 0;
    sumWeightInBonzai_ = 0;
    for (int i=0 ; i<nbEntriesInHeader ; i++){
      treeBonzaiHeader->GetEntry(i);
      totalEventsInBaobab_ += InEvtCount;
      sumWeightInBaobab_ += (InEvtWeightSums->size()>0 ? InEvtWeightSums->at(0) : -99999999);
      sumWeightInBonzai_ += (EvtWeightSums->size()>0 ? EvtWeightSums->at(0) : -99999999);
    }
    cout << "total events in baobab = " << totalEventsInBaobab_ << endl;
    cout << "sum weight in baobab = " << sumWeightInBaobab_ << endl;
    cout << "sum weight in bonzais = " << sumWeightInBonzai_ << endl;

  }
  return;
}
#endif // #ifdef HZZ2l2nuLooper_cxx

#endif
