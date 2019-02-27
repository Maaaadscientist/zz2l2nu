#ifndef trigger_h
#define trigger_h

#include <vector>

#include <TTreeReaderArray.h>


namespace trigger
{
  enum {DoubleMu, SingleMu, DoubleE, HighPtE, SingleE, EMu, SinglePhoton, MC_DiLepton, MC_Photon}; //List of triggers used for our analysis

  using Prescales = TTreeReaderArray<unsigned>;

  int passTrigger(int trig, ULong64_t TrigHltDiMu, ULong64_t TrigHltMu, ULong64_t TrigHltDiEl, ULong64_t TrigHltEl, ULong64_t TrigHltElMu, ULong64_t TrigHltPhot, Prescales const &TrigHltDiMu_prescale, Prescales const &TrigHltMu_prescale, Prescales const &TrigHltDiEl_prescale, Prescales const &TrigHltEl_prescale, Prescales const &TrigHltElMu_prescale, Prescales const &TrigHltPhot_prescale, double selectedPhotonPt = 0);
}

#endif
