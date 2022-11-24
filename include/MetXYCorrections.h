#ifndef HZZ2L2NU_INCLUDE_METXYCORRECTIONS_H_
#define HZZ2L2NU_INCLUDE_METXYCORRECTIONS_H_

#include <utility>

#include <TMath.h>
#include <TString.h>


namespace metXYCorrections {

// See: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections?rev=72#xy_Shift_Correction_MET_phi_modu

/// XY corrections to Pt-miss, code from https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/XYMETCorrection_withUL17andUL18andUL16.h
std::pair<double,double> METXYCorr_Met_MetPhi(double uncormet, double uncormet_phi, int runnb, TString year, bool isMC, int npv, bool isUL, bool ispuppi);

}  // namespace metXYCorrections

#endif  // HZZ2L2NU_INCLUDE_METXYCORRECTIONS_H_
