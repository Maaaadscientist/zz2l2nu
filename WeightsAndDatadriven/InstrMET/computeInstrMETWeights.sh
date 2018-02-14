#!/bin/bash


####################################################################
############################   Options   ###########################
####################################################################
function load_options() {
  #pipefail: the return value of a pipeline is the status of the last command to exit with a non-zero status, or zero if no command exited with a non-zero status
  #this is needed for retry function to work properly
  set -o pipefail
  
  SLEEP_TIME_QSTAT=60 # in seconds
  SLEEP_TIME=60 #in seconds
  CHECK_TIME=6 #check qstat every CHECK_TIME*SLEEP_TIME_QSTAT
  
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
  base_path="$CMSSW_BASE/src/shears/HZZ2l2nu/"
  instrMET_path="WeightsAndDatadriven/InstrMET/"
  full_path="${base_path}${instrMET_path}"
  
  launchAnalysis_step1="${full_path}launchAnalysis_computeInstrMETWeights_step1.sh"
  launchAnalysis_step3="${full_path}launchAnalysis_computeInstrMETWeights_step3.sh"
  listSamplesToRun_HZZ="${full_path}listSamplesToRun_dataOnly.txt"
  listSamplesToRun_Photon="${full_path}listSamplesToRun_InstrMET_dataOnly.txt"

  # The various suffixes
  suffix_step1_HZZ="computeInstrMET_DiLeptonData"
  suffix_step1_InstrMET="computeInstrMET_PhotonData_NoWeight"
  suffix_step3_InstrMET="computeInstrMET_PhotonData_NVtx_WeightApplied"

  # Log file
  logFile="${base_path}OUTPUTS/fullAnalysis_computeInstrMETWeights.$(current_time).log"

  # Running options (default)
  express="0"   # set to 1 to run on express all the time, or use the --express option
  noLocalCopy="0" # set to 1 to NOT do a local copy and stream the ROOT files
  step="all" # by default, run on all 5 steps of the analysis (0, 1, 2, 3, 4)

}


####################################################################
##########################   Help (usage)   ########################
####################################################################
function usage() {
  printf "$BLUE NAME $DEF \n\tcomputeInstrMETWeights.sh - All-in-one launcher of the analysis for the HZZ2l2nu group\n"
  printf "\n\t%-5b  %-40b\n"  "$MAG -h/-help/--help $DEF"            "print this help"
  printf "\n\t%-5b  %-40b\n"  "$MAG -p/-publish/--publish $DEF"      "publish plots depending on the analysis type"
  printf "\n$BLUE SYNOPSIS $DEF\n"
  printf "\n\t%-5b\n"         "./doFullAnalysis.sh $YEL [ANALYSIS_TYPE] $DEF $MAG [LOCAL_COPY] $DEF $RED [EXPRESS] $DEF"
  printf "\n$YEL ANALYSIS_TYPE $DEF\n"
  printf "\n\t%-5b  %-40b\n"  "$YEL HZZanalysis $DEF (default)"  "launch the full (step 1, 2 and 3) 'HZZanalysis' analysis (default option if no arguments)"
  printf "\n\t%-5b  %-40b\n"  "$YEL InstrMET $DEF"               "launch the full (step 1, 2 and 3) 'InstrMET' analysis"
  printf "\n$MAG LOCAL_COPY $DEF\n"
  printf "\n\t%-5b  %-40b\n"  "$MAG 0 $DEF (default)"  "jobs will do a local copy on their node first. This makes them less sensitive to bandwidth issue (default option if no arguments)"
  printf "\n\t%-5b  %-40b\n"  "$MAG 1 $DEF"            "jobs will read in streaming their ROOT files"
  printf "\n$RED EXPRESS $DEF\n"
  printf "\n\t%-5b  %-40b\n"  "$RED 0 $DEF (default)"  "launch jobs on the localgrid (i.e. the normal) queue (default option if no arguments)"
  printf "\n\t%-5b  %-40b\n"  "$RED 1 $DEF"            "launch jobs on the express queue"
}

####################################################################
####################   Miscellaneous Functions   ###################
####################################################################

function current_time() {
  local datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
  echo "$datestamp"
}

function send_mail(){
  export LD_LIBRARY_PATH="/lib64/:$LD_LIBRARY_PATH"
  mailAddress=$(grep $USER ${base_path}Tools/userInfo.db | awk  '{print $2}')
  if [[ $mailAddress == *"@"* ]]; then
    echo -e "$I Mail address found in the db : $mailAddress"
  elif [[ $mailAddress == "no" ]]; then
    echo -e "$I User asked to not received emails"
    return 0
  else
    mailAddress=$(ldapsearch -LLL -x uid=$USER mail | sed -n 's/^[ \t]*mail:[\t]*\(.*\)/\1/p')
  fi
  sed -i "1s/^/Subject: Jobs for SHEARS $(current_time)\n/" $logFile
  sendmail $mailAddress < $logFile
}

function check_running_process() {
  #Kill the process if it is running already...
  me=$(basename $0);
  for pid in $(pidof -x $me); do
      if [ $pid != $$ ]; then
        echo -e "$E An instance of this script is still running in the background. Should I kill it before launching the script? [N/y]?"
        read answer
        if [[ $answer == "y" ]]; then
          kill -9 $pid
          echo -e "$I Instance $pid has been killed"
        fi
      fi
  done

}

# Retries a command on failure.
# $1 - the max number of attempts
# $2... - the command to run
function retry() {
    baseWaitingTime=10
    local -r -i max_attempts="$1"; shift
    local -r cmd="$@"
    local -i attempt_num=1

    until $cmd
    do
        if (( attempt_num == max_attempts ))
        then
            echo -e "$E Attempt $attempt_num failed and there are no more attempts left!" >&2
            exit 5
        else
            echo -e "$W Attempt $attempt_num failed! Trying again in $(( baseWaitingTime * attempt_num )) seconds..." >&2
            sleep $(( baseWaitingTime * attempt_num ))
            (( attempt_num++ ))
        fi
    done
}


####################################################################
#######################   Introduction Text   ######################
####################################################################
function print_introduction_text(){
  # Steps name
  declare -A stepName
  stepName[all]='Clean - run on DiLepton/Photon Data - compute first weight for Nvtx - rerun on photon data with this weight - compute final weight for Pt'
  stepName[0]='Cleaning'
  stepName[1]='Run On DiLepton and Photon Data Only'
  stepName[2]='Compute Weight NVtx'
  stepName[3]='Re Run On Photon Data With NVtx Weights'
  stepName[4]='Compute Final Weights Pt'

  # The text
  echo -e "$I For $YEL help $DEF just add $YEL -h $DEF"
  echo -e "$I Don't forget only one instance of this script can run at the same time."
  echo -e "$I If you need to kill the script, just relaunch this script: it will ask you automatically if you want to kill the running script."
  echo -e "$I You have chosen the following options:"
  echo -e "\t $MAG step: ${step}, i.e. ${stepName[$step]} $DEF"
  echo -e "\t $MAG express: $express $DEF"
  echo -e "\t $MAG noLocalCopy: $noLocalCopy $DEF"

}

####################################################################
#######################   Prepare Scripts   ########################
####################################################################
function prepare_scripts() {
  if [[ $step == "all" || $step == "0" || $step == "1" || $step == "3" ]]; then
    echo -e "$I $(current_time) Preparing scripts..."
    echo -e "$I $(current_time) Removing old scripts..."
    rm -f $launchAnalysis_step1 $launchAnalysis_step3 $listSamplesToRun_HZZ $listSamplesToRun_Photon

    echo -e "$I $(current_time) Copying scripts and dataset lists..."
    cp ../../launchAnalysis.sh $launchAnalysis_step1
    cp ../../listSamplesToRun.txt $listSamplesToRun_HZZ
    cp ../../listSamplesToRun_InstrMET.txt $listSamplesToRun_Photon
  
    echo -e "$I $(current_time) Modifying scripts and dataset lists..."
    #First dataset list names
    sed -i "s|^listDataset=\".*\"$|listDataset=\"${instrMET_path}listSamplesToRun_dataOnly.txt\"|" $launchAnalysis_step1
    sed -i "s|^listDataset_InstrMET=\".*\"$|listDataset_InstrMET=\"${instrMET_path}listSamplesToRun_InstrMET_dataOnly.txt\"|" $launchAnalysis_step1
    #Then we make a copy of the launchAnalysis because the suffix will need to be different between step1 and step3
    if [ $step != "1" ]; then cp $launchAnalysis_step1 $launchAnalysis_step3; fi
    #Then we change the suffix according to the step
    sed -i "s|^suffix=\".*\"$|suffix=\"${suffix_step1_HZZ}\"|" $launchAnalysis_step1
    sed -i "s|^suffix_InstrMET=\".*\"$|suffix_InstrMET=\"${suffix_step1_InstrMET}\"|" $launchAnalysis_step1
    if [ $step != "1" ]; then sed -i "s|^suffix_InstrMET=\".*\"$|suffix_InstrMET=\"${suffix_step3_InstrMET}\"|" $launchAnalysis_step3; fi

    sed -i '/^Bonzais-.*MC.*$/d' $listSamplesToRun_HZZ
    sed -i '/^Bonzais-.*MC.*$/d' $listSamplesToRun_Photon
    echo -e "$I $(current_time) Done."
  fi

}

function backup_previousWeights() {
  if [[ (-f "${full_path}InstrMET_weight_NVtx_vs_pt.root") || (-f "${full_path}InstrMET_weight_pt.root") ]]; then
    if ! [ "$step" == "0" ]; then
      backupForWeight="${full_path}backupForWeight_$(current_time)"
      mkdir -p $backupForWeight
      echo -e "$W The previous weight files (if not needed by the current step) have been put in a backup folder here: $backupForWeight"
    fi
    
    if [[ $step == "all" || $step == "1" || $step == "2" ]]; then #for the first steps we don't need weight files
      mv InstrMET_weight_NVtx_vs_pt.root $backupForWeight
      mv InstrMET_weight_pt.root $backupForWeight

    elif [[ $step == "3" || $step == "4" ]]; then #but for those steps, we need the 1st weight file, so we just move the second weight file
      mkdir -p $backupForWeight
      mv InstrMET_weight_pt.root $backupForWeight

    fi
  fi

}


####################################################################
######################   Waiting Functions   #######################
####################################################################
function check_number_of_remaining_jobs_to_send_from_bigSubmission() { #this is here in case a job failed to be sent in the first steps of big-submission... and so we wait for it
  #This function takes one argument: the suffix of the job where we want to do the check
  theSuffix=$1
  if [ -z "$theSuffix" ]; then
    echo -e "$E Error, you used the function 'check_number_of_remaining_jobs_to_send_from_bigSubmission' without specifying an argument!"
    return 0
  fi

  while [ $(wc -l < ${base_path}OUTPUTS/${theSuffix}/sendJobs_${theSuffix}.cmd) -gt 0 ]
  do
    echo -e "$I $(current_time) There are still $(wc -l < ${base_path}OUTPUTS/${theSuffix}/sendJobs_${theSuffix}.cmd) jobs to send for ${theSuffix}"
    sleep 60
  done
  if [ $(grep -c -e '^qsub ' $(ls -Art ${base_path}OUTPUTS/${theSuffix}/big-submission-*.err | tail -n 1) ) -gt 0 ]; then
    retryCounter=0
    while [ $(grep -c -e '^qsub ' $(ls -Art ${base_path}OUTPUTS/${theSuffix}/big-submission-*.err | tail -n 1) ) -gt 0 ]
    do
      echo -e "$W There are jobs that failed to be submitted in ${base_path}OUTPUTS/${theSuffix}. Let's wait a bit to see if they manage to be submitted"
      sleep 60
      if [ $retryCounter  == 7 ]; then
        echo -e "$E Big-submission didn't manage to send all jobs. We stop here!"
        send_mail
        return 0
      fi
      retryCounter=$((retryCounter+1))
    done
  fi

}

function check_if_jobs_are_done() {
  #This function takes three arguments: the suffix of the job where we want to do the check; the folder where we want to check the output; and the folder that we should use for comparison to know if there are remaining jobs or not.
  theSuffix="$1"
  outputFolderToCheck="$2"
  inputFolderToCompareTo="$3"
  if [[ (-z "$theSuffix") || (-z "$outputFolderToCheck") || (-z "$inputFolderToCompareTo") ]]; then
    echo -e "$E Error, you used the function 'check_if_jobs_are_done' without specifying three arguments!"
    return 0
  fi
  
  totalJobs=$(retry 5 eval ls -1 ${base_path}OUTPUTS/${theSuffix}/${inputFolderToCompareTo} | wc -l)
  if [ $? == 5 ]; then
    send_mail
    return 1
  fi
  sleptTime=1  #don't make it start at 0
  while [ $(getRemainingJobs $theSuffix $outputFolderToCheck $totalJobs) -gt 0 ]
  do
    echo -e "$I $(current_time) There are $(getRemainingJobs $theSuffix $outputFolderToCheck $totalJobs) jobs remaining" 

    if (( $sleptTime % ($SLEEP_TIME_QSTAT*$CHECK_TIME/$SLEEP_TIME) == 0 ))
    then
      echo "Checking with qstat to see if there are still running/pending jobs, or if they crashed..."
      nofJobs=$(getNumJobsOnCE)
      echo -e "$I $(current_time) There are $nofJobs jobs running/pending on the cluster"
      if (( $nofJobs == 0 ))
      then
        echo -e "No jobs are running or pending in the grid, I guess some jobs failed!"
        if (( $(getRemainingJobs $theSuffix $outputFolderToCheck $totalJobs) == 0 ))
        then
          echo "No, it's fine."
        else
          echo -e "$E Some jobs have failed! Exiting"
          exit 3
        fi
      fi
    fi

    sleep $SLEEP_TIME
    sleptTime=$((sleptTime + 1))
  done

}

function getRemainingJobs() {
  theSuffix=$1
  outputFolderToCheck=$2
  totalJobs=$3
  jobsDone=$(retry 5 eval ls -1 ${base_path}OUTPUTS/${theSuffix}/${outputFolderToCheck} | wc -l)
  if [ $? == 5 ]; then
    send_mail
    return 1
  fi
  echo $(($totalJobs-$jobsDone)) #number of running/remaining jobs
}

function fileIsFresh(){
  file=$1
  fileTimeInSeconds=$(timeout -s 9 3 ls -l --time-style=+%s $file|awk '{print $6}')
  currentTimeInSeconds=$(date +%s)
  let currentTimeInSeconds=$currentTimeInSeconds-$SLEEP_TIME_QSTAT

  if [ $fileTimeInSeconds -gt $currentTimeInSeconds ];then
    return 0
  else
    return 1
  fi
}

function getNumJobsOnCE(){
  nameOfJob=$1
  # Setting default for number of jobs
  nJobs=-1

  # NFS file that is used first
  qstatFile='/group/log/dumpOfFullQstat'

  # Checking if NFS file is accessible and not empty
  timeout -s 9 3 ls $qstatFile &> /dev/null
  if [ $? -eq 0 ] && [ $(wc -l $qstatFile|awk '{print $1}') -gt 0 ];then

    # Checking if NFS file has been updated recently
    if fileIsFresh $qstatFile;then
      if [ -z "$nameOfJob" ]; then
        nJobs=$(cat $qstatFile|grep -e Job_Owner|grep -c $USER)
      else
        nJobs=$(cat $qstatFile|grep -e Job_Owner|grep $nameOfJob|grep -c $USER)
      fi
    fi
  fi

  # If NFS is not accessible, or file is not fresh
  if [ $nJobs -eq -1 ];then

    # Doing a standard qstat
    if [ -z "$nameOfJob" ]; then
      nJobs=$(qstat -u $USER |grep $USER|wc -l)
    else
      nJobs=$(qstat -u $USER |grep $USER|grep $nameOfJob|wc -l)
    fi

    # If qstat command fails, sleep then start again
    if [ $? -ne 0 ];then
      sleep 20
      if [ -z "$nameOfJob" ]; then
        nJobs=$(qstat -u $USER |grep $USER|wc -l)
      else
        nJobs=$(qstat -u $USER |grep $USER|grep $nameOfJob|wc -l)
      fi

      # If qstat fails again, set it to 9999 so we know it failed
      if [ $? -ne 0 ];then
        let nJobs=9999
      fi

    fi
  fi

  echo $nJobs
}


####################################################################
#######################   Step 0 - cleaning   ######################
####################################################################
function cleaning() {
  echo -e "$I $(current_time) Starting step 0: full cleaning..."
  echo "a" | source $launchAnalysis_step1 0 #clean HZZanalysis with step 1 suffixes
  echo "a" | source $launchAnalysis_step1 0 InstrMET #clean InstrMET with step 1 suffixes
  echo "a" | source $launchAnalysis_step3 0 InstrMET #clean InstrMET with step 3 suffixes
  echo -e "$I $(current_time) Cleaning done."
}

####################################################################
#########   Step 1 - Run On DiLepton and Photon Data Only   ########
####################################################################
function runOnDiLeptonAndPhotonDataOnly() {
  echo -e "$I $(current_time) Starting step 1: Run On DiLepton and Photon Data Only..."
 
  echo -e "$I $(current_time) Launching HZZanalysis analysis..."
  yes | source $launchAnalysis_step1 1 HZZanalysis $noLocalCopy $express
  if [ $? -eq 0 ]; then
    echo -e "$E Step 1 failed for HZZanalysis, exiting."
    send_mail
    return 0
  fi
  
  echo -e "$I $(current_time) Launching InstrMET analysis..."
  yes | source $launchAnalysis_step1 1 InstrMET $noLocalCopy $express
  if [ $? -eq 0 ]; then
    echo -e "$E Step 1 failed for InstrMET, exiting."
    send_mail
    return 0
  fi
 
  echo -e "$I $(current_time) Waiting for the jobs from step 1 to be over..."
  sleep 60
  check_number_of_remaining_jobs_to_send_from_bigSubmission ${suffix_step1_HZZ}
  check_number_of_remaining_jobs_to_send_from_bigSubmission ${suffix_step1_InstrMET}

  outputFolderToCheck="OUTPUTS"
  inputFolderToCompareTo="JOBS/scripts"
  check_if_jobs_are_done "${suffix_step1_HZZ}" "$outputFolderToCheck" "$inputFolderToCompareTo"
  echo -e "$I $(current_time) Jobs for ${suffix_step1_HZZ} are done, sending harvesting..."
  send_harvesting ${suffix_step1_HZZ} HZZanalysis
  check_if_jobs_are_done "${suffix_step1_InstrMET}" "$outputFolderToCheck" "$inputFolderToCompareTo"
  echo -e "$I $(current_time) Jobs for ${suffix_step1_InstrMET} are done, sending harvesting..."
  send_harvesting ${suffix_step1_InstrMET} InstrMET

  echo -e "$I $(current_time) Waiting for the the harvesting to be over for both samples..."
  sleep 60
  outputFolderToCheck="MERGED | grep -v output_Data.root "
  inputFolderToCompareTo="OUTPUTS | grep _0.root "
  check_if_jobs_are_done "${suffix_step1_HZZ}" "$outputFolderToCheck" "$inputFolderToCompareTo"
  check_if_jobs_are_done "${suffix_step1_InstrMET}" "$outputFolderToCheck" "$inputFolderToCompareTo"

  echo -e "$I $(current_time) The harvesting is done for both samples and can be found in ${suffix_step1_HZZ} and ${suffix_step1_InstrMET}." 
}


function send_harvesting() {
  #This function takes two arguments: the suffix of the job where we want to do the check; and the folder where we want to check the output
  theSuffix=$1
  analysisType=$2
  if [[ (-z "$theSuffix") || (-z "$analysisType") ]]; then
    echo -e "$E Error, you used the function 'send_harvesting' without specifying two arguments!"
    return 0
  fi


  prepare_jobs_for_express > ${base_path}OUTPUTS/$theSuffix/harvesting.sh
  echo "echo \"y\" | sh $launchAnalysis_step1 2 $analysisType" >> ${base_path}OUTPUTS/$theSuffix/harvesting.sh
  retryCounter=0
  while [ $retryCounter -lt 3 ]; do
    if qsub -q express -l walltime=00:30:00 -j oe -o ${base_path}OUTPUTS/$theSuffix/ ${base_path}OUTPUTS/$theSuffix/harvesting.sh 2>&1 | grep -q 'qsub'; then
      echo -e "$W Failed to submit to the grid, retry in 30s"
      retryCounter=$((retryCounter+1))
      sleep 30
    else
      echo -e "$I Harvesting for $analysisType submitted"
      retryCounter=10
    fi
  done
  if [ $retryCounter  == 3 ]; then
    echo -e "$E Failed 3 times to send jobs, exiting"
    send_mail
    return 0
  fi

}

function prepare_jobs_for_express() {
  echo "source \$VO_CMS_SW_DIR/cmsset_default.sh"
  echo "export SCRAM_ARCH=slc6_amd64_gcc530"
  echo "export INITDIR=$CMSSW_BASE/src/shears/HZZ2l2nu"
  echo "cd \$INITDIR"
  echo "hostname ;"
  echo "date;"

}

####################################################################
#################   Step 2 - Compute Weight NVtx   #################
####################################################################
function computeWeightNVtx() {
  echo -e "$I $(current_time) Starting step 2: Compute Weight NVtx..."
  root -l -q 'macroToComputeInstrMETWeights.C++(1)'
  echo -e "$I $(current_time) Step 2 is done."

}


####################################################################
#######   Step 3 - Re Run On Photon Data With NVtx Weights   #######
####################################################################
function reRunOnPhotonDataWithNVtxWeights() {
  echo "test"
#Penser a faire un backup des anciens poids. QUand on lance cette fonction, on backup les anciens poids et on supprime l'ancienne position des poids. Le script n'en trouve plus et est donc content


}


####################################################################
###############   Step 4 - Compute Final Weights Pt   ##############
####################################################################
function computeFinalWeightsPt() {
  echo "test"
#Penser a faire un backup des anciens poids. QUand on lance cette fonction, on backup les anciens poids et on supprime l'ancienne position des poids. Le script n'en trouve plus et est donc content


}


####################################################################
####################################################################
####################################################################
###############  #####  #######   #######  ###  ####  ##############
###############   ###   ######  #  ######  ###   ###  ##############
###############  # # #  #####  ###  #####  ###  # ##  ##############
###############  ## ##  ###           ###  ###  ## #  ##############
###############  #####  ####  #####  ####  ###  ###   ##############
###############  #####  ####  #####  ####  ###  ####  ##############
####################################################################
####################################################################
####################################################################

function main() {
  prepare_scripts
  backup_previousWeights
  if [[ $step == "all" || $step == "0" ]]; then cleaning; fi
  if [[ $step == "all" || $step == "1" ]]; then runOnDiLeptonAndPhotonDataOnly; fi
  if [[ $step == "all" || $step == "2" ]]; then computeWeightNVtx; fi
  if [[ $step == "all" || $step == "3" ]]; then reRunOnPhotonDataWithNVtxWeights; fi
  if [[ $step == "all" || $step == "4" ]]; then computeFinalWeightsPt; fi
  send_mail
}



load_options
for arg in "$@"
do
  case $arg in -h|-help|--help) usage  ; exit 0 ;; esac
  case $arg in -e|-express|--express) express="1"  ;; esac
  case $arg in -nlc|-noLocalCopy|--noLocalCopy) noLocalCopy="1"  ;; esac
  case $arg in all|0|1|2|3|4) step="$arg"  ;; esac
done

check_running_process

print_introduction_text
echo -e "$W Do you wish to compute Instr. MET weights with those options? [${RED}N${DEF}/${GREEN}y${DEF}]"
read answer
if [[ $answer == "y" ]];
then
  mkdir -p ${base_path}OUTPUTS 

  echo -e "$I Script launched! The log are available here: tail -f $YEL ${logFile} $DEF"
  echo -e "$I Open it with 'tail -f' for realtime update or with 'less -R' to benefit from the colour output."

  main &> $logFile &


fi
