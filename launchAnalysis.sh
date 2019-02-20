#!/usr/bin/bash


if [ $(ps -p $$ | grep bash |wc -l) -ne 1 ]
then
  #echo "Not using bash, relaunching with bash!"
  exec bash "$0" "$@"
fi

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
  
  #NRB
  listDataset_NRB="listSamplesToRun_NRB.txt"
  suffix_NRB="NRB"
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

  # Paths definition
  path="$HZZ2L2NU_BASE/"
  build_dir_path="$path/build"

  # Running options (default)
  doExpress=""   # set to "--express" to run on express all the time, or use the --express option
  doLocalCopy="--localCopy" # set to "--localCopy" to NOT do a local copy and stream the ROOT files
  step="printHelpAndExit" #default: one has to select a step. if one forgot, then just print help and exit
  analysisType="default" #default: run the HZZanalysis with MC (no datadriven)
  systType="no"

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
  printf "\n\t%-5b  %-40b\n"  "$RED 5 $DEF"  "Compute yields, prepare datacards for $GREEN 'ANALYSIS_TYPE' $DEF and publish them" 
  printf "\n$GREEN ANALYSIS_TYPE $DEF\n"
  printf "\n\t%-5b  %-40b\n"  "$GREEN HZZdatadriven $DEF (default)"  "perform the actions described above for the analysis 'HZZdatadriven' (default option if no arguments)" 
  printf "\n\t%-5b  %-40b\n"  "$GREEN HZZanalysis $DEF"  "perform the actions described above for the analysis 'HZZanalysis'" 
  printf "\n\t%-5b  %-40b\n"  "$GREEN InstrMET $DEF"               "perform the actions described above for the analysis 'InstrMET'" 
  printf "\n\t%-5b  %-40b\n"  "$GREEN TnP $DEF"                    "perform the actions described above for the analysis 'Tag and Probe'" 
  printf "\n\t%-5b  %-40b\n"  "$GREEN NRB $DEF"                    "perform the actions described above for the analysis 'Non-resonant Bkg.'" 
  printf "\n$MAG OPTIONS $DEF\n"
  printf "\n\t%-5b  %-40b\n"  "$MAG --syst YOUR_SYST $DEF "  "Run the analysis on YOUR_SYST (see systList.txt for the names; 'all' to run on all systs in this file)" 
  printf "\n\t%-5b  %-40b\n"  "$MAG -nlc/-noLocalCopy/--noLocalCopy $DEF"  "jobs will NOT be copied locally on their node first. This makes them more sensitive to bandwidth issue (default option is to perform a local copy to avoid streaming)" 
  printf "\n\t%-5b  %-40b\n"  "$MAG -e/-express/--express $DEF (default)"  "launch jobs on the express queue" 
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
  case $arg in -e|-express|--express) doExpress="--express"; shift  ;; esac #default: launch on cream02, not on express.
  case $arg in -nlc|-noLocalCopy|--noLocalCopy) doLocalCopy=""; shift  ;; esac #default: do a local copy, don't stream the ROOT files
  case $arg in 0|1|2|3|4|5) step="$arg" ;shift ;; esac
  case $arg in HZZanalysis|HZZdatadriven|InstrMET|TnP|NRB) analysisType="$arg"; shift  ;; esac
  case $arg in --syst) systType="$2"; shift;shift  ;; esac
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
elif [  $analysisType == "NRB" ];then
  listDataset=$listDataset_NRB
  suffix=$suffix_NRB
  analysis="--doNRBAnalysis"
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
    if [[ $answer == "a" ]]; then rm -rf "$build_dir_path"; fi
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
    if ! [ -x "$(command -v runHZZanalysis)" ]; then
      echo -e "$I runHZZanalysis was not found, I'm going to compile it..."
      mkdir -p "$build_dir_path"
      cd "$build_dir_path"
      cmake ..
      make -j 4
      make install
      cd -
    fi
    if ! [ -x "$(command -v runHZZanalysis)" ]; then
      echo -e "$E The compilation failed! Exiting..."
      exit 0
    else
      mkdir -p ${pathAndSuffix}
      if [ $analysisType ==  "HZZdatadriven" ]; then
        if ! [ -f ${path}WeightsAndDatadriven/InstrMET/InstrMET_weight_NVtx.root ] || ! [ -f ${path}WeightsAndDatadriven/InstrMET/InstrMET_weight_pt.root ] || ! [ -f ${path}WeightsAndDatadriven/InstrMET/InstrMET_lineshape_mass.root ] ; then
          echo -e "$E The Weight Files used for Instr. MET don't exist! Please produce them (see ${path}WeightsAndDatadriven/InstrMET/computeInstrMETWeights.sh). Exiting..."
          return 0
        else
          sed '/^Bonzais.*-DYJets.*$/d' ${path}${listDataset_HZZ} > ${pathAndSuffix}$(basename ${listDataset_HZZ}) #Copy HZZ list without DYJets MC (since we are datadriven)
          sed '/^Bonzais.*-GJets_.*$/d' ${path}${listDataset_Photon} | sed '/^Bonzais.*-QCD_.*$/d' > ${pathAndSuffix}$(basename ${listDataset_Photon}) #Copy Photon withoug GJets and QCD
          ${path}Tools/prepareAllJobs.py --listDataset ${pathAndSuffix}$(basename ${listDataset_HZZ}) --suffix $suffix $analysis $doLocalCopy $doExpress --syst $systType
          ${path}Tools/prepareAllJobs.py --listDataset ${pathAndSuffix}$(basename ${listDataset_Photon}) --suffix $suffix --isPhotonDatadriven $analysis $doLocalCopy $doExpress --syst $systType
        fi
      else
        cp ${path}${listDataset} ${pathAndSuffix}$(basename ${listDataset})
        ${path}Tools/prepareAllJobs.py --listDataset ${pathAndSuffix}$(basename ${listDataset}) --syst $systType --suffix $suffix $analysis $doLocalCopy $doExpress
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
    jobsSent=$(eval ls -1 ${path}OUTPUTS/${suffix}/JOBS/scripts | wc -l)
    jobsSucceeded=$(eval ls -1 ${path}OUTPUTS/${suffix}/OUTPUTS | wc -l)
    jobsFailed=$(($jobsSent-$jobsSucceeded))
    if [ $jobsSent == 0 ] || [ $jobsFailed -ne 0 ]; then
      echo -e "$E $RED$jobsFailed jobs failed!$DEF Here is the list:"
      diff \
        <(ls ${path}OUTPUTS/${suffix}/OUTPUTS | cut  -d "." -f -1) \
        <(ls ${path}OUTPUTS/${suffix}/JOBS/scripts | cut  -d "." -f -1 | cut -d "_" -f 2-) \
        | grep '>' | cut -c3-
      echo -e "$W Do you want to merge plots anyway? [N/y]"
      read mergeAnyway
      if [[ $mergeAnyway == "y" ]]; then
        echo -e "$W OK. I hope you know what you're doing."
      else
        echo -e "$I Stopping here."
        exit 0
      fi
    fi
    if [ $analysisType ==  "HZZdatadriven" ]; then
      ${path}Tools/prepareAllJobs.py --listDataset ${pathAndSuffix}$(basename ${listDataset_HZZ}) --suffix $suffix $analysis --syst $systType --harvest 
      ${path}Tools/prepareAllJobs.py --listDataset ${pathAndSuffix}$(basename ${listDataset_Photon}) --suffix $suffix --isPhotonDatadriven $analysis --syst $systType --harvest
      root -l -q -b "Tools/harvestInstrMET.C(\"$suffix\",\"$systType\")" #Harvest Instr.MET
    else
      ${path}Tools/prepareAllJobs.py --listDataset ${pathAndSuffix}$(basename ${listDataset}) --suffix $suffix $analysis --syst $systType --harvest
    fi
    if [ $systType == "all" ]; then
      echo -e "$I Merging is done. Removing all temporary files and renaming them to have a clean output."
      cd ${pathAndSuffix}/MERGED/
      ls | grep -v '_final' | xargs rm
      rename '_final' '' *
      cd -
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
    #filesToMerge=$(eval ls -1 ${path}OUTPUTS/${suffix}/OUTPUTS | grep _0.root | grep -v _up.*root | grep -v _down.*root | wc -l)
    #if [ $analysisType == "HZZdatadriven" ]; then filesToMerge=$(($filesToMerge+1)); fi #There is one more file to merge in the datadriven case: the Instr.MET
    #filesSucceeded=$(eval ls -1 ${path}OUTPUTS/${suffix}/MERGED | grep -v output.*_Data.root | grep -v output.*_final.root | wc -l)
    #filesFailed=$(($filesToMerge-$filesSucceeded))
    #if [ $filesToMerge == 0 ] || [ $filesFailed -ne 0  ]; then
    #  echo -e "$E $RED$filesFailed dataset types were not merged!$DEF Do you want to perform a data/MC comparison anyway? [N/y]"
    #  read doItAnyway
    #  if [[ $doItAnyway == "y" ]]; then
    #    echo -e "$W OK. I hope you know what you're doing."
    #  else
    #    echo -e "$I Stopping here."
    #    exit 0
    #  fi
    #fi

    rm -rf ${path}OUTPUTS/${suffix}/PLOTS
    mkdir -p ${path}OUTPUTS/${suffix}/PLOTS
    root -l -q -b "${path}dataMCcomparison.C(\"$analysisType\",\"$suffix\")"
    if [ $systType != "no" ]; then root -l -q -b "Tools/InstrMET_syst_study.C(\"$suffix\")"; fi;
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
    cp ${base}Tools/index.php ~/public_html/SHEARS_PLOTS/$suffix/*/.
    echo -e "$I Your plots are available in ~/public_html/SHEARS_PLOTS/$suffix/, i.e. on http://homepage.iihe.ac.be/~$USER/SHEARS_PLOTS/$suffix/"
  fi
fi

##############
### STEP 5 ###
##############
if [[ $step == 5 ]]; then
  echo -e "$W Do you want to compute yields (and publish them) and prepare datacards for $RED'$analysisType'$DEF with suffix $RED'${suffix}'$DEF? [N/y]?"
  read answer
  if [[ $answer == "y" ]];
  then
    if [ $systType != "all" ]; then
      echo -e "$E You can compute yields and prepare datacards only if all syst are produced. So please launch it that way (--syst all)."
      exit 0
    fi
    if [ $analysisType ==  "HZZdatadriven" ]; then
      python ${base}Tools/prepareDataCards.py --suffix $suffix --dataDriven
    else
      python ${base}Tools/prepareDataCards.py --suffix $suffix
    fi
    mkdir -p ~/public_html
    chmod 755 ~/public_html
    mkdir -p ~/public_html/SHEARS_PLOTS
    mkdir -p ~/public_html/SHEARS_PLOTS/$suffix
    rm -rf ~/public_html/SHEARS_PLOTS/$suffix/YIELDS
    mkdir -p ~/public_html/SHEARS_PLOTS/$suffix/YIELDS
    ln -s  ${path}OUTPUTS/${suffix}/PLOTS/YIELDS/* ~/public_html/SHEARS_PLOTS/$suffix/YIELDS/.
    cp ${base}Tools/index.php ~/public_html/SHEARS_PLOTS/$suffix/.
    cp ${base}Tools/index.php ~/public_html/SHEARS_PLOTS/$suffix/*/.
    echo -e "$I Yields table and datacards are available in ~/public_html/SHEARS_PLOTS/$suffix/YIELDS/, i.e. on http://homepage.iihe.ac.be/~$USER/SHEARS_PLOTS/$suffix/YIELDS"
  fi
fi

