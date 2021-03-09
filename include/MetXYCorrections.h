#ifndef HZZ2L2NU_INCLUDE_METXYCORRECTIONS_H_
#define HZZ2L2NU_INCLUDE_METXYCORRECTIONS_H_

#include <utility>

#include <TMath.h>


namespace metXYCorrections {

/// XY corrections to Pt-miss, code from https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/XYMETCorrection.h
std::pair<double,double> METXYCorr_Met_MetPhi(double originalMet, double originalMet_phi, int runnb, int year, bool isMC, int npv);

}  // namespace metXYCorrections

#endif  // HZZ2L2NU_INCLUDE_METXYCORRECTIONS_H_
