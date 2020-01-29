#!/bin/bash


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
  # Different versions of lists of input datasets
  listDataset_Main="$HZZ2L2NU_BASE/listSamplesToRun_${year}.txt"
  listDataset_MELA="$HZZ2L2NU_BASE/listSamplesToRun_MELA_${year}.txt"
  listDataset_InstrMET="$HZZ2L2NU_BASE/listSamplesToRun_InstrMET_${year}.txt"
  listDataset_NRB="$HZZ2L2NU_BASE/listSamplesToRun_NRB_${year}.txt"
  master_config="${year}.yaml"
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

  # Running options (default)
  doLocalCopy="--local-copy" # set to "--local-copy" to NOT do a local copy and stream the ROOT files
  step="printHelpAndExit" #default: one has to select a step. if one forgot, then just print help and exit
  analysisType="default" #default: run the HZZanalysis with MC (no datadriven)
  systType="no"
  doMelaReweight="" #default: signal samples will NOT be reweighted using MELA
  year="2016" #default
  #launch suffix
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
  printf "\n\t%-5b  %-40b\n"  "$GREEN NRB $DEF"                    "perform the actions described above for the analysis 'Non-resonant Bkg.'" 
  printf "\n\t%-5b  %-40b\n"  "$GREEN DileptonTrees $DEF" "MC-based analysis that produces trees" 
  printf "\n$MAG OPTIONS $DEF\n"
  printf "\n\t%-5b  %-40b\n"  "$MAG -d/--task-dir DIR $DEF "  "Directory for the task. Defaults to OUTPUTS/<suffix>"
  printf "\n\t%-5b  %-40b\n"  "$MAG --syst YOUR_SYST $DEF "  "Run the analysis on YOUR_SYST (see config/syst.yaml for the names; 'all' to run on all systs in this file)" 
  printf "\n\t%-5b  %-40b\n"  "$MAG -nlc/-noLocalCopy/--noLocalCopy $DEF"  "jobs will NOT be copied locally on their node first. This makes them more sensitive to bandwidth issue (default option is to perform a local copy to avoid streaming)" 
  printf "\n\t%-5b  %-40b\n"  "$MAG --mela $DEF"  "Run on MELA-weighted samples"
  printf "\n\t%-5b  %-40b\n"  "$MAG --year YEAR $DEF "  "The dataset year (default option is 2016)"
  printf "\n\t%-5b  %-40b\n"  "$MAG --add_options OPTIONS $DEF "  "Additional options to be forwarded to the main executable"
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
  case $arg in -nlc|-noLocalCopy|--noLocalCopy) doLocalCopy=""; shift  ;; esac #default: do a local copy, don't stream the ROOT files
  case $arg in 0|1|2|3|4|5) step="$arg" ;shift ;; esac
  case $arg in HZZanalysis|HZZdatadriven|InstrMET|NRB|DileptonTrees) analysisType="$arg"; shift ;; esac
  case $arg in -d|--task-dir) task_dir="$2"; shift; shift ;; esac
  case $arg in --mela) doMelaReweight="--doMelaReweight"; shift  ;; esac #default: no mela reweight 
  case $arg in --syst) systType="$2"; shift;shift  ;; esac
  case $arg in --add-options) add_options="$2"; shift; shift ;; esac
  case $arg in --year) year="$2"; shift; shift ;; esac
done
user_configuration

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
  analysis="Main"
  if [ -n "$doMelaReweight" ]; then
    listDataset=$listDataset_MELA
  else
    listDataset=$listDataset_Main
  fi
  suffix="HZZ"
elif [  $analysisType == "InstrMET" ];then
  analysis="InstrMET"
  listDataset=$listDataset_InstrMET
  suffix="InstrMET_ALL_REWEIGHTING_APPLIED"
elif [  $analysisType == "NRB" ];then
  analysis="NRB"
  listDataset=$listDataset_NRB
  suffix="NRB"
elif [ $analysisType ==  "HZZdatadriven" ]; then
  analysis="Main"
  if [ -n "$doMelaReweight" ]; then
    listDataset_HZZ=$listDataset_MELA
  else
    listDataset_HZZ=$listDataset_Main
  fi
  listDataset_Photon=$listDataset_InstrMET
  suffix="HZZdatadriven"
elif [ $analysisType == "DileptonTrees" ]; then
  analysis="DileptonTrees"
  listDataset=$listDataset_NRB
  suffix="DileptonTrees"
else
  echo -e "$E The analysis '$analysisType' does not exist"
  step=-1 #small trick to just go out of the function
fi

if [ -z "$task_dir" ]; then
  task_dir=OUTPUTS/$suffix
fi



##############
### STEP 0 ###
##############
if [[ $step == 0 ]]; then 
  #analysis cleanup
  echo -e "$W ALL DATA WILL BE LOST IN $RED'$task_dir'$DEF and $RED'~/public_html/SHEARS_PLOTS/plots_$suffix'$DEF! [N/y]?"
  read answer
  if [[ $answer == "y" ]];
  then
    echo "CLEANING UP..."
    rm -rf ${task_dir} ~/public_html/SHEARS_PLOTS/plots_$suffix
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
      echo -e "$E Executable not found! Exiting..."
      exit 0
    else
      mkdir -p ${task_dir}
      if [ $analysisType ==  "HZZdatadriven" ]; then
        if ! [ -f ${HZZ2L2NU_BASE}/WeightsAndDatadriven/InstrMET/InstrMET_weight_NVtx.root ] || ! [ -f ${HZZ2L2NU_BASE}/WeightsAndDatadriven/InstrMET/InstrMET_weight_pt.root ] || ! [ -f ${HZZ2L2NU_BASE}/WeightsAndDatadriven/InstrMET/InstrMET_lineshape_mass.root ] ; then
          echo -e "$E The Weight Files used for Instr. MET don't exist! Please produce them (see ${HZZ2L2NU_BASE}/WeightsAndDatadriven/InstrMET/computeInstrMETWeights.sh). Exiting..."
          return 0
        else
          sed '/^Bonzais.*-DYJets.*$/d' ${listDataset_HZZ} > ${task_dir}/$(basename ${listDataset_HZZ}) #Copy HZZ list without DYJets MC (since we are datadriven)
          sed '/^Bonzais.*-GJets_.*$/d' ${listDataset_Photon} | sed '/^Bonzais.*-QCD_.*$/d' > ${task_dir}/$(basename ${listDataset_Photon}) #Copy Photon withoug GJets and QCD
          prepare_jobs.py ${task_dir}/$(basename ${listDataset_HZZ}) -d $task_dir -a $analysis $doLocalCopy --syst $systType "--add-options=$add_options" --config $master_config
          prepare_jobs.py ${task_dir}/$(basename ${listDataset_Photon}) -d $task_dir -a $analysis --dd-photon $doLocalCopy --syst $systType "--add-options=$add_options" --config $master_config
        fi
      else
        cp ${listDataset} ${task_dir}/$(basename ${listDataset})
        prepare_jobs.py ${task_dir}/$(basename ${listDataset}) -d $task_dir -a $analysis $doLocalCopy --syst $systType "--add-options=$add_options" --config $master_config
      fi
      cd $task_dir
      big-submission send_jobs.sh
      cd - > /dev/null
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
    jobsSent=$(eval ls -1 "${task_dir}/jobs/scripts" | wc -l)
    jobsSucceeded=$(eval ls -1 "${task_dir}/output" | wc -l)
    jobsFailed=$(($jobsSent-$jobsSucceeded))
    if [ $jobsSent == 0 ] || [ $jobsFailed -ne 0 ]; then
      echo -e "$E $RED$jobsFailed jobs failed!$DEF Here is the list:"
      diff \
        <(ls "${task_dir}/output" | cut  -d "." -f -1) \
        <(ls "${task_dir}/jobs/scripts" | cut  -d "." -f -1 | cut -d "_" -f 2-) \
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
      harvest.py ${task_dir}/$(basename ${listDataset_HZZ}) -d $task_dir -a $analysis --syst $systType --config $master_config
      harvest.py ${task_dir}/$(basename ${listDataset_Photon}) -d $task_dir -a $analysis --dd-photon --syst $systType --config $master_config
      root -l -q -b "$HZZ2L2NU_BASE/Tools/harvestInstrMET.C(\"$task_dir\",\"$systType\")" #Harvest Instr.MET
    else
      harvest.py ${task_dir}/$(basename ${listDataset}) -d $task_dir -a $analysis --syst $systType --config $master_config
    fi
    if [ "$systType" == "all" ] && [ "$analysisType" != "DileptonTrees" ]; then
      echo -e "$I Merging is done. Removing all temporary files and renaming them to have a clean output."
      cd "${task_dir}/merged" && ls | grep -v '_final' | xargs rm
      rename '_final' '' *
      cd - > /dev/null
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
    rm -rf "${task_dir}/plots"
    mkdir -p "${task_dir}/plots"
    root -l -q -b "${HZZ2L2NU_BASE}/dataMCcomparison.C(\"$analysisType\",\"$task_dir\",\"$doMelaReweight\")"
    if [ $systType != "no" ]; then root -l -q -b "${HZZ2L2NU_BASE}/Tools/InstrMET_syst_study.C(\"$task_dir\")"; fi;
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
    ln -s  ${task_dir}/plots/* ~/public_html/SHEARS_PLOTS/$suffix/.
    cp $HZZ2L2NU_BASE/Tools/index.php ~/public_html/SHEARS_PLOTS/$suffix/.
    cp $HZZ2L2NU_BASE/Tools/index.php ~/public_html/SHEARS_PLOTS/$suffix/*/.
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
      if [ $doMelaReweight != "" ]; then
        python $HZZ2L2NU_BASE/Tools/prepareDataCards.py $task_dir/merged --InstrMETdataDrivenOnly --signalType ggH --massPoint 800 --higgsWidth 100 --yields --datacard
      else
        python $HZZ2L2NU_BASE/Tools/prepareDataCards.py $task_dir/merged --InstrMETdataDrivenOnly --yields
      fi
    else
      python $HZZ2L2NU_BASE/Tools/prepareDataCards.py $task_dir/merged --yields
    fi
    mkdir -p ~/public_html
    chmod 755 ~/public_html
    mkdir -p ~/public_html/SHEARS_PLOTS
    mkdir -p ~/public_html/SHEARS_PLOTS/$suffix
    rm -rf ~/public_html/SHEARS_PLOTS/$suffix/YIELDS
    mkdir -p ~/public_html/SHEARS_PLOTS/$suffix/YIELDS
    ln -s  ${task_dir}/plots/yields/* ~/public_html/SHEARS_PLOTS/$suffix/YIELDS/.
    cp $HZZ2L2NU_BASE/Tools/index.php ~/public_html/SHEARS_PLOTS/$suffix/.
    cp $HZZ2L2NU_BASE/Tools/index.php ~/public_html/SHEARS_PLOTS/$suffix/*/.
    echo -e "$I Yields table and datacards are available in ~/public_html/SHEARS_PLOTS/$suffix/YIELDS/, i.e. on http://homepage.iihe.ac.be/~$USER/SHEARS_PLOTS/$suffix/YIELDS"
  fi
fi

