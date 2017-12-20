#!/usr/bin/env bash

# Colors
BLUE='\033[1;34m'
RED='\033[1;31m'
GREEN='\033[1;32m'
YEL='\033[1;33m'
MAG='\033[1;35m'
DEF='\033[0;m'

# Nice printing
I="$GREEN[INFO] $DEF"
E="$RED[ERROR] $DEF"
W="$YEL[WARN] $DEF"


if [[ $# -eq 0 ]]; then 
  printf "$BLUE NAME $DEF \n\tlaunchAnalysis.sh - Launcher of the analysis for the HZZ2l2nu group\n"
  printf "\n$BLUE SYNOPSIS $DEF\n"
  printf "\n\t%-5b\n"         "./launchAnalysis.sh $GREEN [OPTION] $DEF $YEL [ANALYSIS_TYPE] $DEF $MAG [LOCAL_COPY] $DEF $RED [EXPRESS] $DEF" 
  printf "\n$GREEN OPTIONS $DEF\n" 
  printf "\n\t%-5b  %-40b\n"  "$GREEN 0 $DEF"  "completely clean up the directory (inputs AND outputs)" 
  printf "\n\t%-5b  %-40b\n"  "$GREEN 1 $DEF"  "run $YEL 'ANALYSIS_TYPE' $DEF on all samples" 
  printf "\n\t%-5b  %-40b\n"  "$GREEN 2 $DEF"  "harvest all samples from $YEL 'ANALYSIS_TYPE' $DEF from step 1" 
  printf "\n\t%-5b  %-40b\n"  "$GREEN 3 $DEF"  "perform data MC comparison for $YEL 'ANALYSIS_TYPE' $DEF" 
  printf "\n$YEL ANALYSIS_TYPE $DEF\n"
  printf "\n\t%-5b  %-40b\n"  "$YEL HZZanalysis $DEF (default)"  "perform the actions described above for the analysis 'HZZanalysis' (default option if no arguments)" 
  printf "\n\t%-5b  %-40b\n"  "$YEL InstrMET $DEF"               "perform the actions described above for the analysis 'InstrMET'" 
  printf "\n\t%-5b  %-40b\n"  "$YEL TnP $DEF"                    "perform the actions described above for the analysis 'Tag and Probe'" 
  printf "\n$MAG LOCAL_COPY $DEF\n"
  printf "\n\t%-5b  %-40b\n"  "$MAG 0 $DEF (default)"  "jobs will do a local copy on their node first. This makes them less sensitive to bandwidth issue (default option if no arguments)" 
  printf "\n\t%-5b  %-40b\n"  "$MAG 1 $DEF"            "jobs will read in streaming their ROOT files" 
  printf "\n$RED EXPRESS $DEF\n"
  printf "\n\t%-5b  %-40b\n"  "$RED 0 $DEF (default)"  "launch jobs on the localgrid (i.e. the normal) queue (default option if no arguments)" 
  printf "\n\t%-5b  %-40b\n"  "$RED 1 $DEF"            "launch jobs on the express queue" 
fi


step=$1   #first variable stores the analysis step to run

if [ $# -gt 1 ] #second variable is for analysis type
then
  analysisType=$2
else
  analysisType="HZZanalysis"
fi

if [ $# -gt 2 ] && [ $3 == "1" ] #third variable is to set the local copy mode on
then
  doLocalCopy=""
else
  doLocalCopy="--localCopy" #we want the default (i.e. value 0) to be the local copy
fi

if [ $# -gt 3 ] && [ $4 == "1" ] #fourth variable is to set the express mode on
then
  doExpress="--express"
else
  doExpress="" #default: launch on cream02, not on express
fi


#####################
### Configuration ###
#####################
#HZZanalysis
listDataset="listSamplesToRun.txt"
suffix="firstTest"

#InstrMET
listDataset_InstrMET="listSamplesToRun_InstrMET.txt"
suffix_InstrMET="firstTest_InstrMET_newBaobabs_withGJetQCDcleaning"

#TnP
listDataset_TnP="TO_BE_ADDED_WHEN_IT_WILL_EXIST"
suffix_TnP="firstTest_TnP"

if [ $analysisType == "HZZanalysis" ];then
  analysis="" #default option
elif [  $analysisType == "InstrMET" ];then
  listDataset=$listDataset_InstrMET
  suffix=$suffix_InstrMET
  analysis="--doInstrMETAnalysis"
elif [  $analysisType == "TnP" ];then
  listDataset=$listDataset_TnP
  suffix=$suffix_TnP
  analysis="--doTnPTree"
else
  echo -e "$E The analysis '$analysisType' does not exist"
  step=-1 #small trick to just go out of the function
fi

#give access to the suffix to the "doFullAnalysisScript"
if [ -f tmp_shared_variables.txt ]; then
  echo $suffix > tmp_shared_variables.txt
fi

if [ "$CMSSW_BASE" == "" ]; then
    echo -e "$W Setting CMSSW environment here (if you don't want to see this all the time, either source the script or do a cmsenv !"
    eval `scramv1 runtime -sh`
    echo "Done!"
fi

##############
### STEP 0 ###
##############
if [[ $step == 0 ]]; then 
  #analysis cleanup
  echo -e "$W ALL DATA WILL BE LOST (answer 'a' if you also want to recompile)! [N/y/a]?"
  read answer
  if [[ $answer == "y" || $answer == "a" ]];
  then
    echo "CLEANING UP..."
    rm -rf JOBS/ OUTPUTS_$suffix runOnBatch_* sendJobs* big-submission-* merged_$suffix plots_$suffix *.sh.o* ~/public_html/SHEARS_PLOTS/plots_$suffix
    if [[ $answer == "a" ]]; then make mrproper; fi
  fi
  echo "Done."
fi #end of step0

##############
### STEP 1 ###
##############
if [[ $step == 1 ]]; then
  echo -e "$W Do you want to launch the analysis $RED'$analysisType'$DEF? [N/y]?"
  read answer
  if [[ $answer == "y" ]];
  then
    if ! [ -f runHZZanalysis ]; then
      echo -e "$I runHZZanalysis was not found, I'm going to compile it..."
      make clean
      make -j4
    fi
    if ! [ -f runHZZanalysis ]; then
      echo -e "$E The compilation failed! Exiting..."
      return 0
    else
      ./prepareAllJobs.py --listDataset $listDataset --suffix $suffix $analysis $doLocalCopy $doExpress
      big-submission sendJobs_${suffix}.cmd
      return 1 #Those lines will complain when using this script with 'sh' but that's not an issue and they are needed for the 'doFullAnalysis'
    fi
  fi

fi

##############
### STEP 2 ###
##############
if [[ $step == 2 ]]; then
  echo -e "$W Do you want to merge all the plots for $RED'$analysisType'$DEF? The previous ones may be deleted!!! [N/y]?"
  read answer
  if [[ $answer == "y" ]];
  then
  ./prepareAllJobs.py --listDataset $listDataset --suffix $suffix --harvest

  fi
fi

##############
### STEP 3 ###
##############
if [[ $step == 3 ]]; then
  echo -e "$W Do you want to perform a data-MC comparison for $RED'$analysisType'$DEF? The previous plots in the plots/ folder will be deleted!!! [N/y]?"
  read answer
  if [[ $answer == "y" ]];
  then
    rm -rf plots_${suffix}
    mkdir plots_${suffix}
    root -l -q -b "dataMCcomparison.C(\"$analysisType\",\"$suffix\")"

  fi
fi
