#define HZZ2l2nuLooper_cxx
#include "../Includes/HZZ2l2nuLooper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void HZZ2l2nuLooper::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   cout << "nb of entries in the input file =" << nentries << endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      cout << "nb of interactions=" << EvtPuCnt << endl;
   }
}
