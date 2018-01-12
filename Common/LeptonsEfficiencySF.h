#ifndef LeptonsEfficiencySF_h
#define LeptonsEfficiencySF_h

#include <iostream>
#include <string>
#include <TMath.h>

namespace llvvElecIdIso { enum ElecIdIso  {Veto, Loose, Medium, Tight, TightIso}; }
//namespace llvvMuonId { enum MuonId  {Loose, Soft, Tight, tkHighPT, TightAndTlkHighPt}; }
/*namespace llvvPhotonId { enum PhotonId  {Loose, Medium, Tight}; }
namespace llvvElecIso{ enum ElecIso {Veto, Loose, Medium, Tight}; }
namespace llvvMuonIso{ enum MuonIso {Loose,Tight, H4lWP, TightBoosted, TightAndTkRelatBoosted}; }*/
namespace CutVersion { enum CutSet {Moriond17Cut}; }


namespace trigAndIDsfs
{
    std::pair<float,float> leptonEffSF(int, float, float, int, int );

}


#endif
