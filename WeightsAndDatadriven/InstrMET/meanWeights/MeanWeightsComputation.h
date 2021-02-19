//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb  7 10:33:03 2020 by ROOT version 6.18/00
// from TTree Vars/
// found on file: ../../../OUTPUTS/PhotonTrees_v4/merged/Data.root
//////////////////////////////////////////////////////////

#ifndef MeanWeightsComputation_h
#define MeanWeightsComputation_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TLorentzVector.h"

class MeanWeightsComputation {
  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.

    std::vector<float> pT_thresholds_ = {0., 55., 82.5, 99., 132., 181.5, 250, 9999};
    //std::vector<float> mT_thresholds_ = {150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 2100, 3000, 9999};
    std::vector<float> mT_thresholds_ = {100, 200, 300, 350, 400, 450, 500, 550, 600, 700, 850, 1000, 1250, 1500, 3000, 9999};
    //Double_t mT_binning_[15] = {150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 2100, 3000};
    Double_t mT_binning_[15] = {100, 200, 300, 350, 400, 450, 500, 550, 600, 700, 850, 1000, 1250, 1500, 3000};
    Int_t n_mT_binning_ = sizeof(mT_binning_)/sizeof(Double_t);

    // Declaration of leaf types
    Int_t           jet_cat;
    Float_t         photon_pt;
    Float_t         ptmiss;
    Float_t         mT;
    Float_t         sm_DjjVBF;
    Int_t           num_pv_good;
    Float_t         trigger_weight;
    Float_t         photon_reweighting;
    Float_t         photon_nvtx_reweighting;

    // List of branches
    TBranch        *b_jet_cat;   //!
    TBranch        *b_photon_pt;   //!
    TBranch        *b_ptmiss;   //!
    TBranch        *b_mT;   //!
    TBranch        *b_sm_DjjVBF;   //!
    TBranch        *b_num_pv_good;   //!
    TBranch        *b_trigger_weight;   //!
    TBranch        *b_photon_reweighting;   //!
    TBranch        *b_photon_nvtx_reweighting;   //!

    MeanWeightsComputation(TTree *tree=0);
    virtual ~MeanWeightsComputation();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
    std::pair<int, int> find_thresholds_binning(float pT, float mT);
};

#endif

#ifdef MeanWeightsComputation_cxx
MeanWeightsComputation::MeanWeightsComputation(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TString base_path = std::string(getenv("HZZ2L2NU_BASE")) + "/";
    TString tree_path = base_path+ "OUTPUTS/PhotonTrees_2017_DataOnly_v10/merged/"; // Path to be updated
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(tree_path+"Data.root");
    if (!f || !f->IsOpen()) {
      f = new TFile(tree_path+"Data.root");
    }
    f->GetObject("Vars",tree);

  }
  Init(tree);
}

MeanWeightsComputation::~MeanWeightsComputation()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t MeanWeightsComputation::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t MeanWeightsComputation::LoadTree(Long64_t entry)
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

void MeanWeightsComputation::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("jet_cat", &jet_cat, &b_jet_cat);
  fChain->SetBranchAddress("photon_pt", &photon_pt, &b_photon_pt);
  fChain->SetBranchAddress("ptmiss", &ptmiss, &b_ptmiss);
  fChain->SetBranchAddress("mT", &mT, &b_mT);
  fChain->SetBranchAddress("sm_DjjVBF", &sm_DjjVBF, &b_sm_DjjVBF);
  fChain->SetBranchAddress("num_pv_good", &num_pv_good, &b_num_pv_good);
  fChain->SetBranchAddress("trigger_weight", &trigger_weight, &b_trigger_weight);
  fChain->SetBranchAddress("photon_reweighting", &photon_reweighting, &b_photon_reweighting);
  fChain->SetBranchAddress("photon_nvtx_reweighting", &photon_nvtx_reweighting, &b_photon_nvtx_reweighting);
  Notify();
}

Bool_t MeanWeightsComputation::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void MeanWeightsComputation::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t MeanWeightsComputation::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef MeanWeightsComputation_cxx
