#define HZZ2l2nuLooper_cxx
#include "../Includes/HZZ2l2nuLooper.h"
#include "../Includes/SmartSelectionMonitor.h"
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>

void HZZ2l2nuLooper::Loop()
{
   if (fChain == 0) return;

   SmartSelectionMonitor mon;
   mon.addHistogram(new TH1F("pile-up",";Number of PU events;Events",50,0,50));

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   cout << "nb of entries in the input file =" << nentries << endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      double weight = 1.;

      //cout << "nb of interactions=" << EvtPuCnt << endl;
      mon.fillHisto("pile-up","all",EvtPuCnt,weight);
   }

   TFile* outFile=TFile::Open("out.root","recreate");
   mon.Write();
   outFile->Close();

}
