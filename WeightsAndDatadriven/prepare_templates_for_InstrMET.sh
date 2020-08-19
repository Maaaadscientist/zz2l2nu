#!/bin/bash

rm template_dilepton_2017_withInstrMET_asimov.root
rm template_photon_2017_bin_*.root
root -l -q 'makePhotonTemplates.C(true)'
root -l -q 'makePhotonTemplates.C(false)'
hadd -f template_dilepton_2017_withInstrMET_asimov.root temporary_SR_bin_*.root
rm temporary_SR_bin*.root
