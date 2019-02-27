#ifndef LooperMain_h
#define LooperMain_h

#include <RoccoR.h>
#include <Utils.h>

#include <cstdlib>
#include <iostream>
#include <fstream>

#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TRandom3.h> 

#include <Options.h>
#include <SmartSelectionMonitor_hzz.h>
#include <TLorentzVectorWithIndex.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

using std::vector;
using namespace std;

class LooperMain {
public :
   TTreeReader fReader;  //!the tree reader
   TChain *fChain = nullptr;   //!pointer to the analyzed TChain
   

   //global variables
   Options const &options_;
   int maxEvents_;
   TString outputFile_;
   int isMC_;
   int isPhotonDatadriven_;
   double sampleXsection_;
   double totalEventsInBaobab_;
   double sumWeightInBaobab_;
   double sumWeightInBonzai_;
   TString syst_;
   bool keepAllControlPlots_;
   bool runOnBaobabs_;
   
   
   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> EvtIsRealData = {fReader, "EvtIsRealData"};
   TTreeReaderValue<UInt_t> EvtNum = {fReader, "EvtNum"};
   TTreeReaderValue<UInt_t> EvtRunNum = {fReader, "EvtRunNum"};
   TTreeReaderValue<Int_t> EvtLumiNum = {fReader, "EvtLumiNum"};
   TTreeReaderValue<Int_t> EvtBxNum = {fReader, "EvtBxNum"};
   TTreeReaderValue<Int_t> EvtVtxCnt = {fReader, "EvtVtxCnt"};
   TTreeReaderValue<Int_t> EvtPuCnt = {fReader, "EvtPuCnt"};
   TTreeReaderValue<Int_t> EvtPuCntTruth = {fReader, "EvtPuCntTruth"};
   TTreeReaderArray<double> EvtWeights = {fReader, "EvtWeights"};
   TTreeReaderValue<Float_t> EvtFastJetRho = {fReader, "EvtFastJetRho"};
   TTreeReaderValue<ULong64_t> TrigMET = {fReader, "TrigMET"};
   TTreeReaderValue<UInt_t> TrigHlt = {fReader, "TrigHlt"};
   TTreeReaderValue<ULong64_t> TrigHltPhot = {fReader, "TrigHltPhot"};
   TTreeReaderValue<ULong64_t> TrigHltDiPhot = {fReader, "TrigHltDiPhot"};
   TTreeReaderValue<ULong64_t> TrigHltMu = {fReader, "TrigHltMu"};
   TTreeReaderValue<ULong64_t> TrigHltDiMu = {fReader, "TrigHltDiMu"};
   TTreeReaderValue<ULong64_t> TrigHltEl = {fReader, "TrigHltEl"};
   TTreeReaderValue<ULong64_t> TrigHltDiEl = {fReader, "TrigHltDiEl"};
   TTreeReaderValue<ULong64_t> TrigHltElMu = {fReader, "TrigHltElMu"};
   TTreeReaderArray<unsigned int> TrigHltPhot_prescale = {fReader, "TrigHltPhot_prescale"};
   TTreeReaderArray<unsigned int> TrigHltDiPhot_prescale = {fReader, "TrigHltDiPhot_prescale"};
   TTreeReaderArray<unsigned int> TrigHltMu_prescale = {fReader, "TrigHltMu_prescale"};
   TTreeReaderArray<unsigned int> TrigHltDiMu_prescale = {fReader, "TrigHltDiMu_prescale"};
   TTreeReaderArray<unsigned int> TrigHltEl_prescale = {fReader, "TrigHltEl_prescale"};
   TTreeReaderArray<unsigned int> TrigHltDiEl_prescale = {fReader, "TrigHltDiEl_prescale"};
   TTreeReaderArray<unsigned int> TrigHltElMu_prescale = {fReader, "TrigHltElMu_prescale"};
   TTreeReaderArray<float> METPt = {fReader, "METPt"};
   TTreeReaderArray<float> METPhi = {fReader, "METPhi"};
   TTreeReaderArray<float> METPtType1 = {fReader, "METPtType1"};
   TTreeReaderArray<float> METPhiType1 = {fReader, "METPhiType1"};
   TTreeReaderArray<float> METPtType1XY = {fReader, "METPtType1XY"};
   TTreeReaderArray<float> METPhiType1XY = {fReader, "METPhiType1XY"};
   TTreeReaderArray<float> METPtRaw = {fReader, "METPtRaw"};
   TTreeReaderArray<float> METPhiRaw = {fReader, "METPhiRaw"};
   TTreeReaderArray<float> METsigx2 = {fReader, "METsigx2"};
   TTreeReaderArray<float> METsigxy = {fReader, "METsigxy"};
   TTreeReaderArray<float> METsigy2 = {fReader, "METsigy2"};
   TTreeReaderArray<float> METsig = {fReader, "METsig"};
   TTreeReaderArray<float> GMETPt = {fReader, "GMETPt"};
   TTreeReaderArray<float> GMETPhi = {fReader, "GMETPhi"};
   TTreeReaderArray<float> GLepDr01Pt = {fReader, "GLepDr01Pt"};
   TTreeReaderArray<float> GLepDr01Eta = {fReader, "GLepDr01Eta"};
   TTreeReaderArray<float> GLepDr01Phi = {fReader, "GLepDr01Phi"};
   TTreeReaderArray<float> GLepDr01E = {fReader, "GLepDr01E"};
   TTreeReaderArray<int> GLepDr01Id = {fReader, "GLepDr01Id"};
   TTreeReaderArray<int> GLepDr01St = {fReader, "GLepDr01St"};
   TTreeReaderArray<int> GLepDr01MomId = {fReader, "GLepDr01MomId"};
   TTreeReaderValue<vector<bool>> GLepDr01Prompt = {fReader, "GLepDr01Prompt"};
   TTreeReaderValue<vector<bool>> GLepDr01TauProd = {fReader, "GLepDr01TauProd"};
   TTreeReaderArray<float> GLepBarePt = {fReader, "GLepBarePt"};
   TTreeReaderArray<float> GLepBareEta = {fReader, "GLepBareEta"};
   TTreeReaderArray<float> GLepBarePhi = {fReader, "GLepBarePhi"};
   TTreeReaderArray<float> GLepBareE = {fReader, "GLepBareE"};
   TTreeReaderArray<int> GLepBareId = {fReader, "GLepBareId"};
   TTreeReaderArray<int> GLepBareSt = {fReader, "GLepBareSt"};
   TTreeReaderArray<int> GLepBareMomId = {fReader, "GLepBareMomId"};
   TTreeReaderValue<vector<bool>> GLepBarePrompt = {fReader, "GLepBarePrompt"};
   TTreeReaderValue<vector<bool>> GLepBareTauProd = {fReader, "GLepBareTauProd"};
   TTreeReaderArray<float> GLepSt3Pt = {fReader, "GLepSt3Pt"};
   TTreeReaderArray<float> GLepSt3Eta = {fReader, "GLepSt3Eta"};
   TTreeReaderArray<float> GLepSt3Phi = {fReader, "GLepSt3Phi"};
   TTreeReaderArray<float> GLepSt3E = {fReader, "GLepSt3E"};
   TTreeReaderArray<int> GLepSt3Id = {fReader, "GLepSt3Id"};
   TTreeReaderArray<int> GLepSt3St = {fReader, "GLepSt3St"};
   TTreeReaderArray<int> GLepSt3Mother0Id = {fReader, "GLepSt3Mother0Id"};
   TTreeReaderArray<int> GLepSt3MotherCnt = {fReader, "GLepSt3MotherCnt"};
   TTreeReaderArray<float> GPhotPt = {fReader, "GPhotPt"};
   TTreeReaderArray<float> GPhotEta = {fReader, "GPhotEta"};
   TTreeReaderArray<float> GPhotPhi = {fReader, "GPhotPhi"};
   TTreeReaderArray<float> GPhotE = {fReader, "GPhotE"};
   TTreeReaderArray<int> GPhotMotherId = {fReader, "GPhotMotherId"};
   TTreeReaderArray<float> GPhotPrompt = {fReader, "GPhotPrompt"};
   TTreeReaderArray<int> GPhotSt = {fReader, "GPhotSt"};
   TTreeReaderArray<float> GPhotIsoEDR03 = {fReader, "GPhotIsoEDR03"};
   TTreeReaderArray<float> GPhotIsoEDR04 = {fReader, "GPhotIsoEDR04"};
   TTreeReaderArray<float> GPhotIsoEDR05 = {fReader, "GPhotIsoEDR05"};
   TTreeReaderArray<float> GPhotIsoSumPtDR03 = {fReader, "GPhotIsoSumPtDR03"};
   TTreeReaderArray<float> GPhotIsoSumPtDR04 = {fReader, "GPhotIsoSumPtDR04"};
   TTreeReaderArray<float> GPhotIsoSumPtDR05 = {fReader, "GPhotIsoSumPtDR05"};
   TTreeReaderArray<float> GLepClosePhotPt = {fReader, "GLepClosePhotPt"};
   TTreeReaderArray<float> GLepClosePhotEta = {fReader, "GLepClosePhotEta"};
   TTreeReaderArray<float> GLepClosePhotPhi = {fReader, "GLepClosePhotPhi"};
   TTreeReaderArray<float> GLepClosePhotE = {fReader, "GLepClosePhotE"};
   TTreeReaderArray<int> GLepClosePhotId = {fReader, "GLepClosePhotId"};
   TTreeReaderArray<int> GLepClosePhotMother0Id = {fReader, "GLepClosePhotMother0Id"};
   TTreeReaderArray<int> GLepClosePhotMotherCnt = {fReader, "GLepClosePhotMotherCnt"};
   TTreeReaderArray<int> GLepClosePhotSt = {fReader, "GLepClosePhotSt"};
   TTreeReaderArray<float> GJetAk04Pt = {fReader, "GJetAk04Pt"};
   TTreeReaderArray<float> GJetAk04Eta = {fReader, "GJetAk04Eta"};
   TTreeReaderArray<float> GJetAk04Phi = {fReader, "GJetAk04Phi"};
   TTreeReaderArray<float> GJetAk04E = {fReader, "GJetAk04E"};
   TTreeReaderArray<float> GJetAk04ChFrac = {fReader, "GJetAk04ChFrac"};
   TTreeReaderArray<int> GJetAk04ConstCnt = {fReader, "GJetAk04ConstCnt"};
   TTreeReaderArray<int> GJetAk04ConstId = {fReader, "GJetAk04ConstId"};
   TTreeReaderArray<float> GJetAk04ConstPt = {fReader, "GJetAk04ConstPt"};
   TTreeReaderArray<float> GJetAk04ConstEta = {fReader, "GJetAk04ConstEta"};
   TTreeReaderArray<float> GJetAk04ConstPhi = {fReader, "GJetAk04ConstPhi"};
   TTreeReaderArray<float> GJetAk04ConstE = {fReader, "GJetAk04ConstE"};
   TTreeReaderArray<float> GJetAk04MatchedPartonID = {fReader, "GJetAk04MatchedPartonID"};
   TTreeReaderArray<float> GJetAk04MatchedPartonDR = {fReader, "GJetAk04MatchedPartonDR"};
   TTreeReaderArray<float> GJetAk08Pt = {fReader, "GJetAk08Pt"};
   TTreeReaderArray<float> GJetAk08Eta = {fReader, "GJetAk08Eta"};
   TTreeReaderArray<float> GJetAk08Phi = {fReader, "GJetAk08Phi"};
   TTreeReaderArray<float> GJetAk08E = {fReader, "GJetAk08E"};
   TTreeReaderArray<float> GJetAk08ChFrac = {fReader, "GJetAk08ChFrac"};
   TTreeReaderArray<int> GJetAk08ConstCnt = {fReader, "GJetAk08ConstCnt"};
   TTreeReaderArray<int> GJetAk08ConstId = {fReader, "GJetAk08ConstId"};
   TTreeReaderArray<float> GJetAk08ConstPt = {fReader, "GJetAk08ConstPt"};
   TTreeReaderArray<float> GJetAk08ConstEta = {fReader, "GJetAk08ConstEta"};
   TTreeReaderArray<float> GJetAk08ConstPhi = {fReader, "GJetAk08ConstPhi"};
   TTreeReaderArray<float> GJetAk08ConstE = {fReader, "GJetAk08ConstE"};
   TTreeReaderArray<float> GJetAk08MatchedPartonID = {fReader, "GJetAk08MatchedPartonID"};
   TTreeReaderArray<float> GJetAk08MatchedPartonDR = {fReader, "GJetAk08MatchedPartonDR"};
   TTreeReaderArray<int> GPdfId1 = {fReader, "GPdfId1"};
   TTreeReaderArray<int> GPdfId2 = {fReader, "GPdfId2"};
   TTreeReaderArray<float> GPdfx1 = {fReader, "GPdfx1"};
   TTreeReaderArray<float> GPdfx2 = {fReader, "GPdfx2"};
   TTreeReaderArray<float> GPdfScale = {fReader, "GPdfScale"};
   TTreeReaderValue<Float_t> GBinningValue = {fReader, "GBinningValue"};
   TTreeReaderValue<Int_t> GNup = {fReader, "GNup"};
   TTreeReaderArray<float> MuPt = {fReader, "MuPt"};
   TTreeReaderArray<float> MuEta = {fReader, "MuEta"};
   TTreeReaderArray<float> MuPhi = {fReader, "MuPhi"};
   TTreeReaderArray<float> MuE = {fReader, "MuE"};
   TTreeReaderArray<unsigned int> MuId = {fReader, "MuId"};
   TTreeReaderArray<unsigned int> MuIdTight = {fReader, "MuIdTight"};
   TTreeReaderArray<unsigned int> MuIdSoft = {fReader, "MuIdSoft"};
   TTreeReaderArray<unsigned int> MuIdHighPt = {fReader, "MuIdHighPt"};
   TTreeReaderArray<unsigned int> MuIdTkHighPt = {fReader, "MuIdTkHighPt"};
   TTreeReaderArray<float> MuCh = {fReader, "MuCh"};
   TTreeReaderArray<float> MuVtxZ = {fReader, "MuVtxZ"};
   TTreeReaderArray<float> MuDxy = {fReader, "MuDxy"};
   TTreeReaderArray<float> MuIsoRho = {fReader, "MuIsoRho"};
   TTreeReaderArray<float> MuPfIso = {fReader, "MuPfIso"};
   TTreeReaderArray<float> MuType = {fReader, "MuType"};
   TTreeReaderArray<float> MuIsoTkIsoAbs = {fReader, "MuIsoTkIsoAbs"};
   TTreeReaderArray<float> MuIsoTkIsoRel = {fReader, "MuIsoTkIsoRel"};
   TTreeReaderArray<float> MuIsoCalAbs = {fReader, "MuIsoCalAbs"};
   TTreeReaderArray<float> MuIsoCombRel = {fReader, "MuIsoCombRel"};
   TTreeReaderArray<float> MuTkNormChi2 = {fReader, "MuTkNormChi2"};
   TTreeReaderArray<int> MuTkHitCnt = {fReader, "MuTkHitCnt"};
   TTreeReaderArray<int> MuMatchedStationCnt = {fReader, "MuMatchedStationCnt"};
   TTreeReaderArray<float> MuDz = {fReader, "MuDz"};
   TTreeReaderArray<int> MuPixelHitCnt = {fReader, "MuPixelHitCnt"};
   TTreeReaderArray<int> MuTkLayerCnt = {fReader, "MuTkLayerCnt"};
   TTreeReaderArray<float> MuPfIsoChHad = {fReader, "MuPfIsoChHad"};
   TTreeReaderArray<float> MuPfIsoNeutralHad = {fReader, "MuPfIsoNeutralHad"};
   TTreeReaderArray<float> MuPfIsoRawRel = {fReader, "MuPfIsoRawRel"};
   TTreeReaderArray<float> MuRelIsoSafeCloseBy = {fReader, "MuRelIsoSafeCloseBy"};
   TTreeReaderArray<float> MuTrkrelIsoCloseByCleanned = {fReader, "MuTrkrelIsoCloseByCleanned"};
   TTreeReaderArray<unsigned int> MuHltMatch = {fReader, "MuHltMatch"};
   TTreeReaderArray<float> ElPt = {fReader, "ElPt"};
   TTreeReaderArray<float> ElEta = {fReader, "ElEta"};
   TTreeReaderArray<float> ElEtaSc = {fReader, "ElEtaSc"};
   TTreeReaderArray<float> ElPhi = {fReader, "ElPhi"};
   TTreeReaderArray<float> ElE = {fReader, "ElE"};
   TTreeReaderArray<unsigned int> ElId = {fReader, "ElId"};
   TTreeReaderArray<float> ElCh = {fReader, "ElCh"};
   TTreeReaderArray<float> ElScRawE = {fReader, "ElScRawE"};
   TTreeReaderArray<float> ElCorrE = {fReader, "ElCorrE"};
   TTreeReaderArray<float> ElEcalIso = {fReader, "ElEcalIso"};
   TTreeReaderArray<float> ElEcalPfIso = {fReader, "ElEcalPfIso"};
   TTreeReaderArray<float> ElMvaTrig = {fReader, "ElMvaTrig"};
   TTreeReaderArray<float> ElMvaNonTrig = {fReader, "ElMvaNonTrig"};
   TTreeReaderArray<float> ElMvaPresel = {fReader, "ElMvaPresel"};
   TTreeReaderArray<float> ElDEtaTkScAtVtx = {fReader, "ElDEtaTkScAtVtx"};
   TTreeReaderArray<float> ElDPhiTkScAtVtx = {fReader, "ElDPhiTkScAtVtx"};
   TTreeReaderArray<float> ElHoE = {fReader, "ElHoE"};
   TTreeReaderArray<float> ElSigmaIetaIeta = {fReader, "ElSigmaIetaIeta"};
   TTreeReaderArray<float> ElSigmaIetaIetaFull5x5 = {fReader, "ElSigmaIetaIetaFull5x5"};
   TTreeReaderArray<float> ElEinvMinusPinv = {fReader, "ElEinvMinusPinv"};
   TTreeReaderArray<float> ElD0 = {fReader, "ElD0"};
   TTreeReaderArray<float> ElDz = {fReader, "ElDz"};
   TTreeReaderArray<int> ElExpectedMissingInnerHitCnt = {fReader, "ElExpectedMissingInnerHitCnt"};
   TTreeReaderArray<int> ElPassConvVeto = {fReader, "ElPassConvVeto"};
   TTreeReaderArray<unsigned int> ElHltMatch = {fReader, "ElHltMatch"};
   TTreeReaderArray<float> ElPfIsoChHad = {fReader, "ElPfIsoChHad"};
   TTreeReaderArray<float> ElPfIsoNeutralHad = {fReader, "ElPfIsoNeutralHad"};
   TTreeReaderArray<float> ElPfIsoIso = {fReader, "ElPfIsoIso"};
   TTreeReaderArray<float> ElPfIsoPuChHad = {fReader, "ElPfIsoPuChHad"};
   TTreeReaderArray<float> ElPfIsoRaw = {fReader, "ElPfIsoRaw"};
   TTreeReaderArray<float> ElPfIsoDbeta = {fReader, "ElPfIsoDbeta"};
   TTreeReaderArray<float> ElPfIsoRho = {fReader, "ElPfIsoRho"};
   TTreeReaderArray<float> ElAEff = {fReader, "ElAEff"};
   TTreeReaderArray<float> ElDr03TkSumPt = {fReader, "ElDr03TkSumPt"};
   TTreeReaderArray<float> ElDr03EcalRecHitSumEt = {fReader, "ElDr03EcalRecHitSumEt"};
   TTreeReaderArray<float> ElDr03HcalTowerSumEt = {fReader, "ElDr03HcalTowerSumEt"};
   TTreeReaderArray<float> TauPt = {fReader, "TauPt"};
   TTreeReaderArray<float> TauEta = {fReader, "TauEta"};
   TTreeReaderArray<float> TauPhi = {fReader, "TauPhi"};
   TTreeReaderArray<float> TauE = {fReader, "TauE"};
   TTreeReaderArray<float> TauCh = {fReader, "TauCh"};
   TTreeReaderArray<unsigned int> TauDecayModeFinding = {fReader, "TauDecayModeFinding"};
   TTreeReaderArray<float> TauCombinedIsolationDeltaBetaCorrRaw3Hits = {fReader, "TauCombinedIsolationDeltaBetaCorrRaw3Hits"};
   TTreeReaderArray<unsigned int> TauDiscMuonLoose = {fReader, "TauDiscMuonLoose"};
   TTreeReaderArray<unsigned int> TauDiscMuonTight = {fReader, "TauDiscMuonTight"};
   TTreeReaderArray<unsigned int> TauDiscElVLoose = {fReader, "TauDiscElVLoose"};
   TTreeReaderArray<unsigned int> TauDiscElLoose = {fReader, "TauDiscElLoose"};
   TTreeReaderArray<unsigned int> TauDiscElVTight = {fReader, "TauDiscElVTight"};
   TTreeReaderArray<unsigned int> TauDiscElTight = {fReader, "TauDiscElTight"};
   TTreeReaderArray<float> PhotPt = {fReader, "PhotPt"};
   TTreeReaderArray<float> PhotEta = {fReader, "PhotEta"};
   TTreeReaderArray<float> PhotPhi = {fReader, "PhotPhi"};
   TTreeReaderArray<float> PhotScRawE = {fReader, "PhotScRawE"};
   TTreeReaderArray<float> PhotScEta = {fReader, "PhotScEta"};
   TTreeReaderArray<float> PhotScPhi = {fReader, "PhotScPhi"};
   TTreeReaderArray<float> PhotIsoEcal = {fReader, "PhotIsoEcal"};
   TTreeReaderArray<float> PhotIsoHcal = {fReader, "PhotIsoHcal"};
   TTreeReaderArray<float> PhotIsoTk = {fReader, "PhotIsoTk"};
   TTreeReaderArray<float> PhotPfIsoChHad = {fReader, "PhotPfIsoChHad"};
   TTreeReaderArray<float> PhotPfIsoNeutralHad = {fReader, "PhotPfIsoNeutralHad"};
   TTreeReaderArray<float> PhotPfIsoPhot = {fReader, "PhotPfIsoPhot"};
   TTreeReaderArray<float> PhotPfIsoPuChHad = {fReader, "PhotPfIsoPuChHad"};
   TTreeReaderArray<float> PhotPfIsoEcalClus = {fReader, "PhotPfIsoEcalClus"};
   TTreeReaderArray<float> PhotPfIsoHcalClus = {fReader, "PhotPfIsoHcalClus"};
   TTreeReaderArray<float> PhotE3x3 = {fReader, "PhotE3x3"};
   TTreeReaderArray<float> PhotE1x5 = {fReader, "PhotE1x5"};
   TTreeReaderArray<float> PhotE2x5 = {fReader, "PhotE2x5"};
   TTreeReaderArray<float> PhotE5x5 = {fReader, "PhotE5x5"};
   TTreeReaderArray<float> PhotSigmaIetaIeta = {fReader, "PhotSigmaIetaIeta"};
   TTreeReaderArray<float> PhotSigmaIetaIphi = {fReader, "PhotSigmaIetaIphi"};
   TTreeReaderArray<float> PhotSigmaIphiIphi = {fReader, "PhotSigmaIphiIphi"};
   TTreeReaderArray<float> PhotHoE = {fReader, "PhotHoE"};
   TTreeReaderArray<float> PhotHadTowOverEm = {fReader, "PhotHadTowOverEm"};
   TTreeReaderArray<float> PhotEtaWidth = {fReader, "PhotEtaWidth"};
   TTreeReaderArray<float> PhotPhiWidth = {fReader, "PhotPhiWidth"};
   TTreeReaderArray<float> PhotR9 = {fReader, "PhotR9"};
   TTreeReaderArray<float> PhotE1x3 = {fReader, "PhotE1x3"};
   TTreeReaderArray<float> PhotE2x2 = {fReader, "PhotE2x2"};
   TTreeReaderArray<float> PhotS4 = {fReader, "PhotS4"};
   TTreeReaderArray<float> PhotE1x5Full5x5 = {fReader, "PhotE1x5Full5x5"};
   TTreeReaderArray<float> PhotE2x5Full5x5 = {fReader, "PhotE2x5Full5x5"};
   TTreeReaderArray<float> PhotE3x3Full5x5 = {fReader, "PhotE3x3Full5x5"};
   TTreeReaderArray<float> PhotE5x5Full5x5 = {fReader, "PhotE5x5Full5x5"};
   TTreeReaderArray<float> PhotSigmaIetaIetaFull5x5 = {fReader, "PhotSigmaIetaIetaFull5x5"};
   TTreeReaderArray<float> PhotR9Full5x5 = {fReader, "PhotR9Full5x5"};
   TTreeReaderArray<unsigned int> PhotId = {fReader, "PhotId"};
   TTreeReaderValue<vector<bool>> PhotHasPixelSeed = {fReader, "PhotHasPixelSeed"};
   TTreeReaderArray<int> PhotPassElVeto = {fReader, "PhotPassElVeto"};
   TTreeReaderArray<float> JetAk04Pt = {fReader, "JetAk04Pt"};
   TTreeReaderArray<float> JetAk04Eta = {fReader, "JetAk04Eta"};
   TTreeReaderArray<float> JetAk04Phi = {fReader, "JetAk04Phi"};
   TTreeReaderArray<float> JetAk04E = {fReader, "JetAk04E"};
   TTreeReaderArray<float> JetAk04Id = {fReader, "JetAk04Id"};
   TTreeReaderValue<vector<bool>> JetAk04PuId = {fReader, "JetAk04PuId"};
   TTreeReaderArray<float> JetAk04PuMva = {fReader, "JetAk04PuMva"};
   TTreeReaderArray<float> JetAk04RawPt = {fReader, "JetAk04RawPt"};
   TTreeReaderArray<float> JetAk04RawE = {fReader, "JetAk04RawE"};
   TTreeReaderArray<float> JetAk04HfHadE = {fReader, "JetAk04HfHadE"};
   TTreeReaderArray<float> JetAk04HfEmE = {fReader, "JetAk04HfEmE"};
   TTreeReaderArray<float> JetAk04ChHadFrac = {fReader, "JetAk04ChHadFrac"};
   TTreeReaderArray<float> JetAk04NeutralHadAndHfFrac = {fReader, "JetAk04NeutralHadAndHfFrac"};
   TTreeReaderArray<float> JetAk04ChEmFrac = {fReader, "JetAk04ChEmFrac"};
   TTreeReaderArray<float> JetAk04NeutralEmFrac = {fReader, "JetAk04NeutralEmFrac"};
   TTreeReaderArray<float> JetAk04ChMult = {fReader, "JetAk04ChMult"};
   TTreeReaderArray<float> JetAk04NeutMult = {fReader, "JetAk04NeutMult"};
   TTreeReaderArray<float> JetAk04ConstCnt = {fReader, "JetAk04ConstCnt"};
   TTreeReaderArray<float> JetAk04Beta = {fReader, "JetAk04Beta"};
   TTreeReaderArray<float> JetAk04BetaClassic = {fReader, "JetAk04BetaClassic"};
   TTreeReaderArray<float> JetAk04BetaStar = {fReader, "JetAk04BetaStar"};
   TTreeReaderArray<float> JetAk04BetaStarClassic = {fReader, "JetAk04BetaStarClassic"};
   TTreeReaderArray<float> JetAk04Rms = {fReader, "JetAk04Rms"};
   TTreeReaderArray<float> JetAk04BTagCsv = {fReader, "JetAk04BTagCsv"};
   TTreeReaderArray<float> JetAk04BTagCsvV1 = {fReader, "JetAk04BTagCsvV1"};
   TTreeReaderArray<float> JetAk04BTagCsvSLV1 = {fReader, "JetAk04BTagCsvSLV1"};
   TTreeReaderArray<float> JetAk04BDiscCisvV2 = {fReader, "JetAk04BDiscCisvV2"};
   TTreeReaderArray<float> JetAk04BDiscJp = {fReader, "JetAk04BDiscJp"};
   TTreeReaderArray<float> JetAk04BDiscBjp = {fReader, "JetAk04BDiscBjp"};
   TTreeReaderArray<float> JetAk04BDiscTche = {fReader, "JetAk04BDiscTche"};
   TTreeReaderArray<float> JetAk04BDiscTchp = {fReader, "JetAk04BDiscTchp"};
   TTreeReaderArray<float> JetAk04BDiscSsvhe = {fReader, "JetAk04BDiscSsvhe"};
   TTreeReaderArray<float> JetAk04BDiscSsvhp = {fReader, "JetAk04BDiscSsvhp"};
   TTreeReaderArray<float> JetAk04PartFlav = {fReader, "JetAk04PartFlav"};
   TTreeReaderArray<float> JetAk04HadFlav = {fReader, "JetAk04HadFlav"};
   TTreeReaderArray<float> JetAk04JecUncUp = {fReader, "JetAk04JecUncUp"};
   TTreeReaderArray<float> JetAk04JecUncDwn = {fReader, "JetAk04JecUncDwn"};
   TTreeReaderArray<int> JetAk04ConstId = {fReader, "JetAk04ConstId"};
   TTreeReaderArray<float> JetAk04ConstPt = {fReader, "JetAk04ConstPt"};
   TTreeReaderArray<float> JetAk04ConstEta = {fReader, "JetAk04ConstEta"};
   TTreeReaderArray<float> JetAk04ConstPhi = {fReader, "JetAk04ConstPhi"};
   TTreeReaderArray<float> JetAk04ConstE = {fReader, "JetAk04ConstE"};
   TTreeReaderArray<int> JetAk04GenJet = {fReader, "JetAk04GenJet"};
   TTreeReaderArray<float> JetAk08Pt = {fReader, "JetAk08Pt"};
   TTreeReaderArray<float> JetAk08Eta = {fReader, "JetAk08Eta"};
   TTreeReaderArray<float> JetAk08Phi = {fReader, "JetAk08Phi"};
   TTreeReaderArray<float> JetAk08E = {fReader, "JetAk08E"};
   TTreeReaderArray<float> JetAk08Id = {fReader, "JetAk08Id"};
   TTreeReaderArray<float> JetAk08RawPt = {fReader, "JetAk08RawPt"};
   TTreeReaderArray<float> JetAk08RawE = {fReader, "JetAk08RawE"};
   TTreeReaderArray<float> JetAk08HfHadE = {fReader, "JetAk08HfHadE"};
   TTreeReaderArray<float> JetAk08HfEmE = {fReader, "JetAk08HfEmE"};
   TTreeReaderArray<float> JetAk08ChHadFrac = {fReader, "JetAk08ChHadFrac"};
   TTreeReaderArray<float> JetAk08NeutralHadAndHfFrac = {fReader, "JetAk08NeutralHadAndHfFrac"};
   TTreeReaderArray<float> JetAk08ChEmFrac = {fReader, "JetAk08ChEmFrac"};
   TTreeReaderArray<float> JetAk08NeutralEmFrac = {fReader, "JetAk08NeutralEmFrac"};
   TTreeReaderArray<float> JetAk08ChMult = {fReader, "JetAk08ChMult"};
   TTreeReaderArray<float> JetAk08ConstCnt = {fReader, "JetAk08ConstCnt"};
   TTreeReaderArray<float> JetAk08BTagCsv = {fReader, "JetAk08BTagCsv"};
   TTreeReaderArray<float> JetAk08BTagCsvV1 = {fReader, "JetAk08BTagCsvV1"};
   TTreeReaderArray<float> JetAk08BTagCsvSLV1 = {fReader, "JetAk08BTagCsvSLV1"};
   TTreeReaderArray<float> JetAk08BDiscCisvV2 = {fReader, "JetAk08BDiscCisvV2"};
   TTreeReaderArray<float> JetAk08BDiscJp = {fReader, "JetAk08BDiscJp"};
   TTreeReaderArray<float> JetAk08BDiscBjp = {fReader, "JetAk08BDiscBjp"};
   TTreeReaderArray<float> JetAk08BDiscTche = {fReader, "JetAk08BDiscTche"};
   TTreeReaderArray<float> JetAk08BDiscTchp = {fReader, "JetAk08BDiscTchp"};
   TTreeReaderArray<float> JetAk08BDiscSsvhe = {fReader, "JetAk08BDiscSsvhe"};
   TTreeReaderArray<float> JetAk08BDiscSsvhp = {fReader, "JetAk08BDiscSsvhp"};
   TTreeReaderArray<float> JetAk08PartFlav = {fReader, "JetAk08PartFlav"};
   TTreeReaderArray<float> JetAk08HadFlav = {fReader, "JetAk08HadFlav"};
   TTreeReaderArray<float> JetAk08JecUncUp = {fReader, "JetAk08JecUncUp"};
   TTreeReaderArray<float> JetAk08JecUncDwn = {fReader, "JetAk08JecUncDwn"};
   TTreeReaderArray<int> JetAk08ConstId = {fReader, "JetAk08ConstId"};
   TTreeReaderArray<float> JetAk08ConstPt = {fReader, "JetAk08ConstPt"};
   TTreeReaderArray<float> JetAk08ConstEta = {fReader, "JetAk08ConstEta"};
   TTreeReaderArray<float> JetAk08ConstPhi = {fReader, "JetAk08ConstPhi"};
   TTreeReaderArray<float> JetAk08ConstE = {fReader, "JetAk08ConstE"};
   TTreeReaderArray<int> JetAk08GenJet = {fReader, "JetAk08GenJet"};
   TTreeReaderArray<float> JetAk08PrunedMass = {fReader, "JetAk08PrunedMass"};
   TTreeReaderArray<float> JetAk08FilteredMass = {fReader, "JetAk08FilteredMass"};
   TTreeReaderArray<float> JetAk08SoftDropMass = {fReader, "JetAk08SoftDropMass"};
   TTreeReaderArray<float> JetAk08TrimmedMass = {fReader, "JetAk08TrimmedMass"};
   TTreeReaderArray<float> JetAk08Tau1 = {fReader, "JetAk08Tau1"};
   TTreeReaderArray<float> JetAk08Tau2 = {fReader, "JetAk08Tau2"};
   TTreeReaderArray<float> JetAk08Tau3 = {fReader, "JetAk08Tau3"};

   LooperMain(Options const &options);
   virtual ~LooperMain();
   virtual void     Loop();
   virtual void     Loop_InstrMET();
   virtual void     Loop_TnP();
   virtual void     Loop_NRB();
   virtual bool     passTrigger(int triggerType);
   virtual void     FillNbEntries(TChain *);
   virtual void     FillTheTChain(TChain *, TString, int, int);
   virtual std::vector<float> *computeCorrectedMuPt(bool);
   virtual int findTheMatchingGenParticle(int indexOfRecoParticle, float maxDeltaR);

  /**
   * \brief Fills histograms with jets passing b-tagging selection
   */
  void FillBTagEfficiency(std::vector<TLorentzVectorWithIndex> selCentralJets,
    std::vector<double> btags, TTreeReaderArray<float> const &JetAk04HadFlav,
    double weight, SmartSelectionMonitor_hzz &mon) const;

private :
   RoccoR *rocCorrect;
   TRandom3 randomGenerator;
};

#if defined(HZZ2l2nuLooper_cxx) || defined(InstrMETLooper_cxx) || defined(TnPLooper_cxx)
LooperMain::LooperMain(Options const &options)
    : options_(options), fChain(0),
      randomGenerator(options.GetAs<unsigned>("seed")) {

  outputFile_ = options_.GetAs<std::string>("output");
  maxEvents_ = options_.GetAs<long long>("max-events");
  keepAllControlPlots_ = options_.Exists("all-control-plots");
  isMC_ = options_.GetAs<bool>("is-mc");
  sampleXsection_  = options_.GetAs<float>("xsec");
  isPhotonDatadriven_ = options_.Exists("dd-photon");
  syst_ = options_.GetAs<std::string>("syst");

  TString const fileName{options.GetAs<std::string>("catalog")};
  int const skipFiles{options_.GetAs<int>("skip-files")};
  int const maxFiles{options_.GetAs<int>("max-files")};

  if (fileName.BeginsWith("Baobab"))
    runOnBaobabs_ = true;
  else
    runOnBaobabs_ = false;

  std::cout << "The Input Catalog is " << fileName << std::endl;
  std::cout << "The output file is " << outputFile_ << std::endl;
  std::cout << "Will run on a max of " << maxEvents_ << " events" << std::endl;

  if (syst_ == "")
    std::cout << "Will not use systematic uncertainties" << std::endl;
  else
    std::cout << "Will use the systematic " << syst_ << std::endl;

  if (isMC_)
    std::cout << "This file is MC with a cross section of " <<
      sampleXsection_ << std::endl;

  totalEventsInBaobab_=-1;
  sumWeightInBaobab_=-1;
  sumWeightInBonzai_=-1;

  //initialize the Roc correction
  std::string const installPath(std::getenv("HZZ2L2NU_BASE"));
  rocCorrect = new RoccoR(installPath + "/data/rcdata.2016.v3/");

  //First  get the tot number of events from the BonzaiHeader
  TChain * chainHeader = new TChain("tupel/BonzaiHeader","");
  FillTheTChain(chainHeader, fileName, skipFiles, maxFiles);
  FillNbEntries(chainHeader);
  delete chainHeader;

  fChain = new TChain("tupel/EventTree","");
  FillTheTChain(fChain, fileName, skipFiles, maxFiles);
  fReader.SetTree(fChain);
}


LooperMain::~LooperMain()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void LooperMain::FillTheTChain(TChain *theChain, TString theInputCatalog, int skipFiles, int maxFiles){
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
    if (!l.size() || l[0] == '#' || l[0] == '*') continue;

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


void LooperMain::FillNbEntries(TChain  *inputChain)
{

  TTree *treeBonzaiHeader = inputChain;
  int   InEvtCount=0;
  vector<double>   *InEvtWeightSums = new std::vector<double>;
  vector<double>   *EvtWeightSums = new std::vector<double>;

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
  delete InEvtWeightSums;
  delete EvtWeightSums;
  return;
}
std::vector<float> *LooperMain::computeCorrectedMuPt(bool isMC){
  std::vector<float> *correctedPt = new std::vector<float>;
  for (unsigned int i=0 ; i < MuPt.GetSize() ; i++){
    if (!(std::abs(MuEta[i])<2.4 && MuPt[i]<200)){ //apply the correction only in its domain of validity 
      correctedPt->push_back(MuPt[i]);
      continue;
    }
    float  momentumScaleCorr = 1;
    if (isMC){
      int genMatch = findTheMatchingGenParticle(i, 0.01); //for muons a deltaR of 0.01 is actually conservative 
      if (genMatch > -1)
        momentumScaleCorr = rocCorrect->kScaleFromGenMC(
          MuCh[i], MuPt[i], MuEta[i], MuPhi[i],
          MuTkLayerCnt[i], GLepBarePt[genMatch],
          randomGenerator.Uniform(), 0, 0);
      else
        momentumScaleCorr = rocCorrect->kScaleAndSmearMC(
          MuCh[i], MuPt[i], MuEta[i], MuPhi[i],
          MuTkLayerCnt[i], randomGenerator.Uniform(),
          randomGenerator.Uniform(), 0, 0);
      correctedPt->push_back(momentumScaleCorr*MuPt[i]);
    }
    else {
      momentumScaleCorr = rocCorrect->kScaleDT(MuCh[i], MuPt[i], MuEta[i], MuPhi[i], 0, 0);
      correctedPt->push_back(momentumScaleCorr*MuPt[i]);
    }
  }
  return correctedPt;
}
int LooperMain::findTheMatchingGenParticle(int indexOfRecoParticle, float maxDeltaR){
  float partPhi = MuPhi[indexOfRecoParticle];
  float partEta = MuEta[indexOfRecoParticle];
  float minDeltaR = 100;
  float indexMatchedGen = -1;
  for (unsigned int i=0 ; i < GLepBareEta.GetSize() ; i++){
    float deltaR = utils::deltaR(partEta, partPhi, GLepBareEta[i], GLepBarePhi[i]);
    if (deltaR < minDeltaR){
      minDeltaR = deltaR;
      indexMatchedGen = i;
    } 
  }
  if (minDeltaR<maxDeltaR) return indexMatchedGen;
  else return -1;
} 
#endif // #if defined(HZZ2l2nuLooper_cxx) || defined(InstrMETLooper_cxx) || defined(TnPLooper_cxx)

#endif
