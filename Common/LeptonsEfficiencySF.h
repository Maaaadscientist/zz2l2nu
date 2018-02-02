#ifndef LeptonsEfficiencySF_h
#define LeptonsEfficiencySF_h

#include <iostream>
#include <string>
#include <TMath.h>

namespace llvvElecRecoIdIso { enum ElecRecoIdIso  {Reco, Veto, Loose, Medium, Tight, TightIso}; }
namespace llvvRecoMuonIdIso { enum MuonRecoIdIso  {Tracking, Loose, Soft, Tight, tkHighPT, TightAndTlkHighPt, TightIso}; }
/*namespace llvvPhotonId { enum PhotonId  {Loose, Medium, Tight}; }
namespace llvvElecIso{ enum ElecIso {Veto, Loose, Medium, Tight}; }
namespace llvvMuonIso{ enum MuonIso {Loose,Tight, H4lWP, TightBoosted, TightAndTkRelatBoosted}; }*/
namespace CutVersion { enum CutSet {Moriond17Cut, Moriond17CutRunGH}; }


namespace trigAndIDsfs
{
    std::pair<float,float> leptonEffSF(int, float, float, int, int );

}


#endif
