#!/usr/bin/env bash

suffix="InstrMET_CRStudyFinal_10April_GJetsLONoCorr"
path="$CMSSW_BASE/src/shears/HZZ2l2nu/WeightsAndDatadriven/"

rm -rf pT_Z_SpikesRemovedAndSigmaCut_passQCDveto_ZMetPhi*
root -l -q -b "${path}computeGJetsNLOWeight.C(\"$suffix\")"
cp ${path}pT_Z_SpikesRemovedAndSigmaCut_passQCDveto_ZMetPhi* ~/public_html/computeGJetsWeights/.
