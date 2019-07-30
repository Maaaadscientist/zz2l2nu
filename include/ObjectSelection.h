#ifndef objectselection_h
#define objectselection_h

#include <string>
#include <vector>

#include <TTreeReaderArray.h>

#include <PhysicsObjects.h>


namespace objectSelection
{
  bool cleanPathologicEventsInPhotons(TString datasetName, float run, float luminosityBlock, float event);
}

#endif
