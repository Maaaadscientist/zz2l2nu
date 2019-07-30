#include <ObjectSelection.h>

#include <algorithm>
#include <limits>

#include <Utils.h>

namespace objectSelection
{

  bool cleanPathologicEventsInPhotons(TString datasetName, float run, float luminosityBlock, float event){
    bool eventShouldBeCleaned = false;
    //We remove events in the spike for those samples 
    if(datasetName.Contains("QCD_HT100to200")){
      //Spike at 190~200GeV (cut if pt > 190GeV)
      if( run == 1 && luminosityBlock == 21997  && event == 32986438 ) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 123682 && event == 185472705) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 133696 && event == 200489234) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 301998 && event == 452875030) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 237717 && event == 356480214) eventShouldBeCleaned = true;
      
      //Spike at ~330-340GeV (cut if pt > 330GeV)
      if( run == 1 && luminosityBlock == 405615 && event == 608258936) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 238040 && event == 356963627) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 192575 && event == 288784917) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 405110 && event == 607502440) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 398170 && event == 597094584) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 242217 && event == 363227739) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 175934 && event == 263829468) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 336765 && event == 505011533) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 252239 && event == 378257731) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 239877 && event == 359718643) eventShouldBeCleaned = true;
    }
    if(datasetName.Contains("QCD_Pt-30to50_EMEnriched")){
      //Spike at ~330-340GeV (cut if pt > 330GeV)
      if( run == 1 && luminosityBlock == 5960 && event == 16066353) eventShouldBeCleaned = true;
    }
    if(datasetName.Contains("QCD_Pt-20toInf_MuEnrichedPt15")){
      //Spike at ~330-340GeV (cut if pt > 330GeV)
      if( run == 1 && luminosityBlock == 181694 && event == 2097414398) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 40564  && event == 3384544677) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 49551  && event == 1742764771) eventShouldBeCleaned = true;
      if( run == 1 && luminosityBlock == 94563  && event == 2145642832) eventShouldBeCleaned = true;
    }
    return eventShouldBeCleaned;
  }

}
