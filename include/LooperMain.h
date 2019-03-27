#ifndef LooperMain_h
#define LooperMain_h

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
   
   
   // Readers to access the data
   TTreeReaderValue<UInt_t> EvtNum = {fReader, "EvtNum"};
   TTreeReaderValue<UInt_t> EvtRunNum = {fReader, "EvtRunNum"};
   TTreeReaderValue<Int_t> EvtLumiNum = {fReader, "EvtLumiNum"};
   TTreeReaderValue<Int_t> EvtVtxCnt = {fReader, "EvtVtxCnt"};
   TTreeReaderValue<Int_t> EvtPuCnt = {fReader, "EvtPuCnt"};
   TTreeReaderValue<Int_t> EvtPuCntTruth = {fReader, "EvtPuCntTruth"};
   TTreeReaderArray<double> EvtWeights = {fReader, "EvtWeights"};
   TTreeReaderValue<Float_t> EvtFastJetRho = {fReader, "EvtFastJetRho"};
   TTreeReaderValue<ULong64_t> TrigMET = {fReader, "TrigMET"};
   TTreeReaderValue<ULong64_t> TrigHltPhot = {fReader, "TrigHltPhot"};
   TTreeReaderValue<ULong64_t> TrigHltMu = {fReader, "TrigHltMu"};
   TTreeReaderValue<ULong64_t> TrigHltDiMu = {fReader, "TrigHltDiMu"};
   TTreeReaderValue<ULong64_t> TrigHltEl = {fReader, "TrigHltEl"};
   TTreeReaderValue<ULong64_t> TrigHltDiEl = {fReader, "TrigHltDiEl"};
   TTreeReaderValue<ULong64_t> TrigHltElMu = {fReader, "TrigHltElMu"};
   TTreeReaderArray<unsigned int> TrigHltPhot_prescale = {fReader, "TrigHltPhot_prescale"};
   TTreeReaderArray<unsigned int> TrigHltMu_prescale = {fReader, "TrigHltMu_prescale"};
   TTreeReaderArray<unsigned int> TrigHltDiMu_prescale = {fReader, "TrigHltDiMu_prescale"};
   TTreeReaderArray<unsigned int> TrigHltEl_prescale = {fReader, "TrigHltEl_prescale"};
   TTreeReaderArray<unsigned int> TrigHltDiEl_prescale = {fReader, "TrigHltDiEl_prescale"};
   TTreeReaderArray<unsigned int> TrigHltElMu_prescale = {fReader, "TrigHltElMu_prescale"};
   TTreeReaderArray<float> METPtType1XY = {fReader, "METPtType1XY"};
   TTreeReaderArray<float> METPhiType1XY = {fReader, "METPhiType1XY"};
   TTreeReaderArray<float> METsigx2 = {fReader, "METsigx2"};
   TTreeReaderArray<float> METsigxy = {fReader, "METsigxy"};
   TTreeReaderArray<float> METsigy2 = {fReader, "METsigy2"};
   TTreeReaderArray<float> METsig = {fReader, "METsig"};
   TTreeReaderArray<float> GLepBarePt = {fReader, "GLepBarePt"};
   TTreeReaderArray<float> GLepBareEta = {fReader, "GLepBareEta"};
   TTreeReaderArray<float> GLepBarePhi = {fReader, "GLepBarePhi"};
   TTreeReaderArray<float> GLepBareE = {fReader, "GLepBareE"};
   TTreeReaderArray<int> GLepBareId = {fReader, "GLepBareId"};
   TTreeReaderArray<int> GLepBareMomId = {fReader, "GLepBareMomId"};
   TTreeReaderArray<float> GPhotPt = {fReader, "GPhotPt"};
   TTreeReaderArray<float> GPhotPrompt = {fReader, "GPhotPrompt"};
   TTreeReaderArray<float> GJetAk04Pt = {fReader, "GJetAk04Pt"};
   TTreeReaderArray<float> GJetAk04Eta = {fReader, "GJetAk04Eta"};
   TTreeReaderArray<float> GJetAk04Phi = {fReader, "GJetAk04Phi"};
   TTreeReaderArray<float> GJetAk04E = {fReader, "GJetAk04E"};
   TTreeReaderArray<float> MuPt = {fReader, "MuPt"};
   TTreeReaderArray<float> MuEta = {fReader, "MuEta"};
   TTreeReaderArray<float> MuPhi = {fReader, "MuPhi"};
   TTreeReaderArray<float> MuE = {fReader, "MuE"};
   TTreeReaderArray<unsigned int> MuId = {fReader, "MuId"};
   TTreeReaderArray<unsigned int> MuIdTight = {fReader, "MuIdTight"};
   TTreeReaderArray<float> MuCh = {fReader, "MuCh"};
   TTreeReaderArray<float> MuPfIso = {fReader, "MuPfIso"};
   TTreeReaderArray<unsigned int> MuHltMatch = {fReader, "MuHltMatch"};
   TTreeReaderArray<float> ElPt = {fReader, "ElPt"};
   TTreeReaderArray<float> ElEta = {fReader, "ElEta"};
   TTreeReaderArray<float> ElEtaSc = {fReader, "ElEtaSc"};
   TTreeReaderArray<float> ElPhi = {fReader, "ElPhi"};
   TTreeReaderArray<float> ElE = {fReader, "ElE"};
   TTreeReaderArray<unsigned int> ElId = {fReader, "ElId"};
   TTreeReaderArray<float> ElCh = {fReader, "ElCh"};
   TTreeReaderArray<float> ElPfIsoRho = {fReader, "ElPfIsoRho"};
   TTreeReaderArray<float> PhotPt = {fReader, "PhotPt"};
   TTreeReaderArray<float> PhotEta = {fReader, "PhotEta"};
   TTreeReaderArray<float> PhotPhi = {fReader, "PhotPhi"};
   TTreeReaderArray<float> PhotScEta = {fReader, "PhotScEta"};
   TTreeReaderArray<float> PhotPfIsoChHad = {fReader, "PhotPfIsoChHad"};
   TTreeReaderArray<float> PhotPfIsoNeutralHad = {fReader, "PhotPfIsoNeutralHad"};
   TTreeReaderArray<float> PhotPfIsoPhot = {fReader, "PhotPfIsoPhot"};
   TTreeReaderArray<float> PhotSigmaIetaIeta = {fReader, "PhotSigmaIetaIeta"};
   TTreeReaderArray<float> PhotSigmaIphiIphi = {fReader, "PhotSigmaIphiIphi"};
   TTreeReaderArray<float> PhotHoE = {fReader, "PhotHoE"};
   TTreeReaderArray<float> PhotR9 = {fReader, "PhotR9"};
   TTreeReaderValue<vector<bool>> PhotHasPixelSeed = {fReader, "PhotHasPixelSeed"};

   LooperMain(Options const &options);
   virtual ~LooperMain();
   virtual void     Loop();
   virtual void     Loop_InstrMET();
   virtual void     Loop_TnP();
   virtual void     Loop_NRB();
   virtual bool     passTrigger(int triggerType);
   virtual void     FillNbEntries(TChain *);
   virtual void     FillTheTChain(TChain *, TString, int, int);

  /**
   * \brief Fills histograms with jets passing b-tagging selection
   */
  void FillBTagEfficiency(std::vector<Jet> const &jets, double weight,
                          SmartSelectionMonitor_hzz &mon) const;

private :
   TRandom3 randomGenerator_;
};

#if defined(HZZ2l2nuLooper_cxx) || defined(InstrMETLooper_cxx) || defined(TnPLooper_cxx)
LooperMain::LooperMain(Options const &options)
    : options_(options), fChain(0),
      randomGenerator_(options.GetAs<unsigned>("seed")) {

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
#endif // #if defined(HZZ2l2nuLooper_cxx) || defined(InstrMETLooper_cxx) || defined(TnPLooper_cxx)

#endif
