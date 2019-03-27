#ifndef objectselection_h
#define objectselection_h

#include <string>
#include <vector>

#include <TTreeReaderArray.h>

#include <PhysicsObjects.h>
#include <TLorentzVectorWithIndex.h>


namespace objectSelection
{
  bool cleanPathologicEventsInPhotons(TString datasetName, float EvtRunNum, float EvtLumiNum, float EvtNum);
}

#endif
