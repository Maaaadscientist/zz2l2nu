# JERC

This subpackage provides means to compute corrections to jet momentum and jet momentum resolution, as well as the corresponding uncertainties.
Inputs are described in the form of standard text files used in the JME POG (file format for [JEC](https://hypernews.cern.ch/HyperNews/CMS/get/jes/423/1.html) and [JER factors](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution?rev=14#Scale_factors_file_format)).


## Jet momentum scale

Source files have been copied from [this directory](https://github.com/miquork/jecsys/tree/194510cedf65259bc4b58092120df2b87e6a3b24/CondFormats/JetMETObjects), with only technical modifications. The latest commit included from the source repository is from 2014-01-15.


## Jet momentum resolution

Source files with code to access JER factors have been copied from `CMSSW_10_2_22` (latest commit on 2020-03-31), from packages `CondFormats/JetMETObjects` (class `JetResolutionObject`) and `JetMETCorrections/Modules` (class `JetResolution`).
