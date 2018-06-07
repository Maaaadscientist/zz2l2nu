#!/usr/bin/env bash

###############################################
##########    /!\ PLEASE READ /!\    ##########
###############################################
#The function user_configuration (below) should be the only function you should touch.
#In this function you can choose on which dataset list you want to run and which suffix you want to give to your submission
function user_configuration() {
  #HZZanalysis
  listDataset="listSamplesToRun.txt"
  suffix="firstTest"
  
  #HZZdatadriven
  #listDataset_HZZ=$listDataset
  #listDataset_Photon=$listDataset_InstrMET
  suffix_HZZdatadriven="HZZdatadriven"
  
  #InstrMET
  listDataset_InstrMET="listSamplesToRun_InstrMET.txt"
  suffix_InstrMET="InstrMET_ALL_REWEIGHTING_APPLIED"
  
  #TnP
  listDataset_TnP="TO_BE_ADDED_WHEN_IT_WILL_EXIST"
  suffix_TnP="firstTest_TnP"
}
###############################################
##########      /!\ The End /!\      ##########
###############################################

function load_options() {
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

  # To be on he safe side, always start by setting the cmsenv. Even if already set, this ensures that it is set at the right place!
  eval `scramv1 runtime -sh`

  # Paths definition
  path="$CMSSW_BASE/src/shears/HZZ2l2nu/"

  # Running options (default)
  doExpress=""   # set to "--express" to run on express all the time, or use the --express option
  doLocalCopy="--localCopy" # set to "--localCopy" to NOT do a local copy and stream the ROOT files
  step="printHelpAndExit" #default: one has to select a step. if one forgot, then just print help and exit
  analysisType="default" #default: run the HZZanalysis with MC (no datadriven)

  #launch suffix
  user_configuration
}

function usage(){
  printf "$BLUE NAME $DEF \n\tlaunchAnalysis.sh - Launcher of the analysis for the HZZ2l2nu group\n"
  printf "\n$BLUE SYNOPSIS $DEF\n"
  printf "\n\t%-5b\n"         "./launchAnalysis.sh $GREEN [ANALYSIS_TYPE] $DEF $RED [STEPS] $DEF $MAG [OPTIONS] $DEF" 
  printf "\n$RED STEPS $DEF\n" 
  printf "\n\t%-5b  %-40b\n"  "$RED 0 $DEF"  "completely clean up the directory (inputs AND outputs)" 
  printf "\n\t%-5b  %-40b\n"  "$RED 1 $DEF"  "run $GREEN 'ANALYSIS_TYPE' $DEF on all samples" 
  printf "\n\t%-5b  %-40b\n"  "$RED 2 $DEF"  "harvest all samples from $GREEN 'ANALYSIS_TYPE' $DEF from step 1" 
  printf "\n\t%-5b  %-40b\n"  "$RED 3 $DEF"  "perform data MC comparison for $GREEN 'ANALYSIS_TYPE' $DEF" 
  printf "\n\t%-5b  %-40b\n"  "$RED 4 $DEF"  "publish the plots for $GREEN 'ANALYSIS_TYPE' $DEF on your website" 
  printf "\n$GREEN ANALYSIS_TYPE $DEF\n"
  printf "\n\t%-5b  %-40b\n"  "$GREEN HZZdatadriven $DEF (default)"  "perform the actions described above for the analysis 'HZZdatadriven' (default option if no arguments)" 
  printf "\n\t%-5b  %-40b\n"  "$GREEN HZZanalysis $DEF"  "perform the actions described above for the analysis 'HZZanalysis'" 
  printf "\n\t%-5b  %-40b\n"  "$GREEN InstrMET $DEF"               "perform the actions described above for the analysis 'InstrMET'" 
  printf "\n\t%-5b  %-40b\n"  "$GREEN TnP $DEF"                    "perform the actions described above for the analysis 'Tag and Probe'" 
  printf "\n$MAG OPTIONS $DEF\n"
  printf "\n\t%-5b  %-40b\n"  "$MAG -nlc/-noLocalCopy/--noLocalCopy $DEF"  "jobs will NOT be copied locally on their node first. This makes them more sensitive to bandwidth issue (default option is to perform a local copy to avoid streaming)" 
  printf "\n\t%-5b  %-40b\n"  "$RED -e/-express/--express $DEF (default)"  "launch jobs on the express queue" 
}

#
##
### MAIN
##
#

load_options
for arg in "$@"
do
  case $arg in -h|-help|--help) usage  ; exit 0 ;; esac
  case $arg in -e|-express|--express) doExpress="--express"  ;; esac #default: launch on cream02, not on express.
  case $arg in -nlc|-noLocalCopy|--noLocalCopy) doLocalCopy=""  ;; esac #default: do a local copy, don't stream the ROOT files
  case $arg in 0|1|2|3|4) step="$arg"  ;; esac
  case $arg in HZZanalysis|HZZdatadriven|InstrMET|TnP) analysisType="$arg"  ;; esac
done

if [ $step == "printHelpAndExit" ]; then
  usage
  exit 0
fi
if [ $analysisType == "default" ]; then
  analysisType="HZZdatadriven"
  echo -e "$W You did not enter a known analysis name, by default the $RED'$analysisType'$DEF will be launched."
fi

#####################
### Configuration ###
#####################

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
elif [ $analysisType ==  "HZZdatadriven" ]; then
  listDataset_HZZ=$listDataset
  listDataset_Photon=$listDataset_InstrMET
  suffix=$suffix_HZZdatadriven
else
  echo -e "$E The analysis '$analysisType' does not exist"
  step=-1 #small trick to just go out of the function
fi

pathAndSuffix=${path}OUTPUTS/${suffix}/

##############
### STEP 0 ###
##############
if [[ $step == 0 ]]; then 
  #analysis cleanup
  echo -e "$W ALL DATA WILL BE LOST IN $RED'OUTPUTS/${suffix}'$DEF and $RED'~/public_html/SHEARS_PLOTS/plots_$suffix'$DEF (answer 'a' if you also want to recompile)! [N/y/a]?"
  read answer
  if [[ $answer == "y" || $answer == "a" ]];
  then
    echo "CLEANING UP..."
    rm -rf ${pathAndSuffix} ~/public_html/SHEARS_PLOTS/plots_$suffix
    if [[ $answer == "a" ]]; then make -C $path mrproper; fi
  fi
  echo "Done."
fi #end of step0

##############
### STEP 1 ###
##############
if [[ $step == 1 ]]; then
  echo -e "$W Do you want to launch the analysis $RED'$analysisType'$DEF with suffix $RED'${suffix}'$DEF? [N/y]?"
  read answer
  if [[ $answer == "y" ]];
  then
    if ! [ -f ${path}runHZZanalysis ]; then
      echo -e "$I runHZZanalysis was not found, I'm going to compile it..."
      make -C $path clean
      make -C $path -j4
    fi
    if ! [ -f ${path}runHZZanalysis ]; then
      echo -e "$E The compilation failed! Exiting..."
      return 0
    else
      mkdir -p ${pathAndSuffix}
      if [ $analysisType ==  "HZZdatadriven" ]; then
        if ! [ -f ${path}WeightsAndDatadriven/InstrMET/InstrMET_weight_NVtx.root ] || ! [ -f ${path}WeightsAndDatadriven/InstrMET/InstrMET_weight_pt.root ] || ! [ -f ${path}WeightsAndDatadriven/InstrMET/InstrMET_lineshape_mass.root ] ; then
          echo -e "$E The Weight Files used for Instr. MET don't exist! Please produce them (see ${path}WeightsAndDatadriven/InstrMET/computeInstrMETWeights.sh). Exiting..."
          return 0
        else
          sed '/^Bonzais.*-DYJets.*$/d' ${path}${listDataset_HZZ} > ${pathAndSuffix}$(basename ${listDataset_HZZ}) #Copy HZZ list without DYJets MC (since we are datadriven)
          sed '/^Bonzais.*-GJets_.*$/d' ${path}${listDataset_Photon} | sed '/^Bonzais.*-QCD_.*$/d' > ${pathAndSuffix}$(basename ${listDataset_Photon}) #Copy Photon withoug GJets and QCD
          ${path}Tools/prepareAllJobs.py --listDataset ${pathAndSuffix}$(basename ${listDataset_HZZ}) --suffix $suffix $analysis $doLocalCopy $doExpress
          ${path}Tools/prepareAllJobs.py --listDataset ${pathAndSuffix}$(basename ${listDataset_Photon}) --suffix $suffix --isPhotonDatadriven $analysis $doLocalCopy $doExpress
        fi
      else
        cp ${path}${listDataset} ${pathAndSuffix}$(basename ${listDataset})
        ${path}Tools/prepareAllJobs.py --listDataset ${pathAndSuffix}$(basename ${listDataset}) --suffix $suffix $analysis $doLocalCopy $doExpress
      fi
      cd ${path}OUTPUTS/${suffix}/
      big-submission sendJobs_${suffix}.cmd
      cd - > /dev/null
      return 1 #Those lines will complain when using this script with 'sh' but that's not an issue and they are needed for the 'doFullAnalysis'
    fi
  fi

fi

##############
### STEP 2 ###
##############
if [[ $step == 2 ]]; then
  echo -e "$W Do you want to merge all the plots for $RED'$analysisType'$DEF with suffix $RED'${suffix}'$DEF? The previous ones may be deleted!!! [N/y]?"
  read answer
  if [[ $answer == "y" ]];
  then
    if [ $analysisType ==  "HZZdatadriven" ]; then
      ${path}Tools/prepareAllJobs.py --listDataset ${pathAndSuffix}$(basename ${listDataset_HZZ}) --suffix $suffix $analysis --harvest
      ${path}Tools/prepareAllJobs.py --listDataset ${pathAndSuffix}$(basename ${listDataset_Photon}) --suffix $suffix --isPhotonDatadriven $analysis --harvest
      root -l -q -b "Tools/harvestInstrMET.C(\"$suffix\")" #Harvest Instr.MET
    else
      ${path}Tools/prepareAllJobs.py --listDataset ${pathAndSuffix}$(basename ${listDataset}) --suffix $suffix $analysis --harvest
    fi
  fi
fi

##############
### STEP 3 ###
##############
if [[ $step == 3 ]]; then
  echo -e "$W Do you want to perform a data-MC comparison for $RED'$analysisType'$DEF with suffix $RED'${suffix}'$DEF? The previous plots in the plots/ folder will be deleted!!! [N/y]?"
  read answer
  if [[ $answer == "y" ]];
  then
    rm -rf ${path}OUTPUTS/${suffix}/PLOTS
    mkdir -p ${path}OUTPUTS/${suffix}/PLOTS
    root -l -q -b "${path}dataMCcomparison.C(\"$analysisType\",\"$suffix\")"

  fi
fi

##############
### STEP 4 ###
##############
if [[ $step == 4 ]]; then
  echo -e "$W Do you want to publish plots for $RED'$analysisType'$DEF with suffix $RED'${suffix}'$DEF? [N/y]?"
  read answer
  if [[ $answer == "y" ]];
  then
    mkdir -p ~/public_html
    chmod 755 ~/public_html
    mkdir -p ~/public_html/SHEARS_PLOTS
    rm -rf ~/public_html/SHEARS_PLOTS/$suffix
    mkdir -p ~/public_html/SHEARS_PLOTS/$suffix
    ln -s  ${path}OUTPUTS/${suffix}/PLOTS/* ~/public_html/SHEARS_PLOTS/$suffix/.
    cp ${base}Tools/index.php ~/public_html/SHEARS_PLOTS/$suffix/.
    echo -e "$I Your plots are available in ~/public_html/SHEARS_PLOTS/$suffix/, i.e. on http://homepage.iihe.ac.be/~$USER/SHEARS_PLOTS/$suffix/"
  fi
fi

