#!/bin/bash
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

function usage(){
  printf "$BLUE NAME $DEF \n\tdoFullAnalysis.sh - All-in-one launcher of the analysis for the HZZ2l2nu group\n"
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
  # Setting default for number of jobs
  nJobs=-1

  # NFS file that is used first
  qstatFile='/group/log/dumpOfFullQstat'

  # Checking if NFS file is accessible and not empty
  timeout -s 9 3 ls $qstatFile &> /dev/null
  if [ $? -eq 0 ] && [ $(wc -l $qstatFile|awk '{print $1}') -gt 0 ];then

    # Checking if NFS file has been updated recently
    if fileIsFresh $qstatFile;then
      nJobs=$(cat $qstatFile|grep -e Job_Owner|grep -c $USER)
    fi
  fi

  # If NFS is not accessible, or file is not fresh
  if [ $nJobs -eq -1 ];then

    # Doing a standard qstat
    nJobs=$(qstat -u $USER |grep $USER|wc -l)

    # If qstat command fails, sleep then start again
    if [ $? -ne 0 ];then
      sleep 20
      nJobs=$(qstat -u $USER |grep $USER|wc -l)

      # If qstat fails again, set it to 9999 so we know it failed
      if [ $? -ne 0 ];then
        let nJobs=9999
      fi

    fi
  fi

  echo $nJobs
}

function getRemainingJobs(){
  theSuffix=$1
  folder_output=$2
  totalJobs=$3
  jobsDone=$(retry 5 ls -1 ${CMSSW_BASE}/src/shears/HZZ2l2nu/OUTPUTS/${theSuffix}/${folder_output} | wc -l) #exportedSuffix is an exported variable from the script launchAnalysis.sh
  if [ $? == 5 ]; then 
    send_mail
    return 1
  fi
  echo $(($totalJobs-$jobsDone)) #number of running/remaining jobs
}

function publish_plots(){
  theSuffix=$1
  doNotSendMailIfFail=$2 #don't give a value this argument for normal behaviour, i.e send a mail even if an error occurs
  datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
  echo -e "$I Creating symbolic link to your public_html folder..."
  plots_to_publish=$(retry 5 ls -1 ${CMSSW_BASE}/src/shears/HZZ2l2nu/OUTPUTS/${theSuffix}/PLOTS/ |wc -l)
  if [ $plots_to_publish -eq 0 ]; then
    echo -e "$E No plots to publish"
    if [ -z "$doNotSendMailIfFail" ]; then send_mail; fi
    exit 4
  else
    mkdir -p ~/public_html
    chmod 755 ~/public_html
    mkdir -p ~/public_html/SHEARS_PLOTS
    rm -rf ~/public_html/SHEARS_PLOTS/plots_${theSuffix}
    mkdir -p ~/public_html/SHEARS_PLOTS/plots_${theSuffix}
    ln -s ${CMSSW_BASE}/src/shears/HZZ2l2nu/OUTPUTS/${theSuffix}/PLOTS/* ~/public_html/SHEARS_PLOTS/plots_${theSuffix}/.
    cp ${CMSSW_BASE}/src/shears/HZZ2l2nu/Tools/index.php ~/public_html/SHEARS_PLOTS/plots_${theSuffix}/.
    echo -e "$I Your plots are available in ~/public_html/SHEARS_PLOTS/plots_${theSuffix}/, i.e. on http://homepage.iihe.ac.be/~$USER/SHEARS_PLOTS/plots_${theSuffix}/"
  fi
}

function send_mail(){
  export LD_LIBRARY_PATH="/lib64/:$LD_LIBRARY_PATH"
  mailAddress=$(grep $USER Tools/userInfo.db | awk  '{print $2}')
  if [[ $mailAddress == *"@"* ]]; then
    echo -e "$I Mail address found in the db : $mailAddress"
  elif [[ $mailAddress == "no" ]]; then
    echo -e "$I User asked to not received emails"
    return 0
  else
    mailAddress=$(ldapsearch -LLL -x uid=$USER mail | sed -n 's/^[ \t]*mail:[\t]*\(.*\)/\1/p')
  fi
  datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
  sed -i "1s/^/Subject: Jobs for SHEARS $datestamp\n/" $logFile
  sendmail $mailAddress < $logFile
}

function main(){
  rm -f OUTPUTS/$theSuffix/prepare_tmp.sh
  rm -f OUTPUTS/$theSuffix/step2_tmp.sh
  rm -f OUTPUTS/$theSuffix/step3_tmp.sh
  
  #0) Full cleaning
  datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
  echo $datestamp
  echo "Starting full cleaning..."
  echo "a" | source launchAnalysis.sh 0 $analysisType

  #1) Launch analysis on cluster
  echo "Starting step 1..."
  yes | source launchAnalysis.sh 1 $analysisType $localCopy $express
  if [ $? -eq 0 ]; then
    echo -e "$E Step 1 failed, exiting."
    send_mail
    return 0
  fi
  sleep 60

  #2) Harvest on express queue
  echo "Waiting for step 1 to be over..."
  #a) prepare environment
  echo "source \$VO_CMS_SW_DIR/cmsset_default.sh" >> OUTPUTS/$theSuffix/prepare_tmp.sh
  echo "export SCRAM_ARCH=slc6_amd64_gcc530" >> OUTPUTS/$theSuffix/prepare_tmp.sh
  echo "export INITDIR=$CMSSW_BASE/src/shears/HZZ2l2nu" >> OUTPUTS/$theSuffix/prepare_tmp.sh
  echo "cd \$INITDIR" >> OUTPUTS/$theSuffix/prepare_tmp.sh
  echo "hostname ;" >> OUTPUTS/$theSuffix/prepare_tmp.sh
  echo "date;" >> OUTPUTS/$theSuffix/prepare_tmp.sh
  #b) check number of jobs
  while [ $(wc -l < OUTPUTS/${theSuffix}/sendJobs_${theSuffix}.cmd) -gt 0 ] #this is here in case a job failed to be sent... and so we wait for him.
  do
    datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
    echo -e "$I [$datestamp] There are still $(wc -l < OUTPUTS/${theSuffix}/sendJobs_${theSuffix}.cmd) jobs to send"
    sleep 60
  done
  if [ $(grep -c -e '^qsub ' $(ls -Art OUTPUTS/${theSuffix}/big-submission-*.err | tail -n 1) ) -gt 0 ]; then
    retryCounter=0
    while [ $(grep -c -e '^qsub ' $(ls -Art OUTPUTS/${theSuffix}/big-submission-*.err | tail -n 1) ) -gt 0 ]
    do
      echo -e "$W There are jobs that failed to be submitted. Let's wait a bit to see if they manage to be submitted"
      sleep 60
      if [ $retryCounter  == 7 ]; then
        echo -e "$E Big-submission didn't manage to send all jobs. We stop here!"
        send_mail
        return 0
      fi
      retryCounter=$((retryCounter+1))
    done
  fi
  folder="OUTPUTS"
  totalJobs=$(retry 5 ls -1 ${CMSSW_BASE}/src/shears/HZZ2l2nu/OUTPUTS/${theSuffix}/JOBS/scripts | wc -l)
  if [ $? == 5 ]; then
    send_mail
    return 1
  fi
  sleptTime=1  #don't make it start at 0
  while [ $(getRemainingJobs $theSuffix $folder $totalJobs) -gt 0 ]
  do
    datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
    echo -e "$I [$datestamp] There are $(getRemainingJobs $theSuffix $folder $totalJobs) jobs remaining" 

    if (( $sleptTime % ($SLEEP_TIME_QSTAT*$CHECK_TIME/$SLEEP_TIME) == 0 ))
    then
      echo "Checking with qstat to see if there are still running/pending jobs, or if they crashed..."
      echo -e "$I [$datestamp] There are $(getNumJobsOnCE) jobs running/pending on the cluster"
      if (( $(getNumJobsOnCE) == 0 ))
      then
        echo -e "No jobs are running or pending in the grid, I guess some jobs failed!"
        if (( $(getRemainingJobs $theSuffix $folder $totalJobs) == 0 ))
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
  echo "All jobs are done, launch step 2 !"  
  cp OUTPUTS/$theSuffix/prepare_tmp.sh OUTPUTS/$theSuffix/step2_tmp.sh
  echo "echo \"y\" | sh launchAnalysis.sh 2 $analysisType" >> OUTPUTS/$theSuffix/step2_tmp.sh
  
  retryCounter=0
  while [ $retryCounter -lt 3 ]; do
    if qsub -q express -l walltime=00:30:00 -j oe -o OUTPUTS/$theSuffix/ OUTPUTS/$theSuffix/step2_tmp.sh 2>&1 | grep -q 'qsub'; then
      echo -e "$W Failed to submit to the grid, retry in 30s"
      retryCounter=$((retryCounter+1))
      sleep 30
    else 
      echo -e "$I Step 2 submitted"
      retryCounter=10
    fi
  done
  if [ $retryCounter  == 3 ]; then
    echo -e "$E Failed 3 times to send jobs, exiting"
    send_mail
    return 0
  fi
  sleep 60

  #3) Do data-MC comparison
  echo "Waiting for step 2 to be over..." 
  folder="MERGED"
  totalJobs=$(retry 5 ls -1 ${CMSSW_BASE}/src/shears/HZZ2l2nu/OUTPUTS/${theSuffix}/OUTPUTS |grep _0.root | wc -l)
  if [ $? == 5 ]; then
    send_mail
    return 1
  fi
  sleptTime=1 #don't make it start at 0
  while [ $(getRemainingJobs $theSuffix $folder $totalJobs) -gt 0 ]
  do
    datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
    echo -e "$I [$datestamp] There are $(getRemainingJobs $theSuffix $folder $totalJobs) datasets to merge remaining" 

    if (( $sleptTime % ($SLEEP_TIME_QSTAT*$CHECK_TIME/$SLEEP_TIME) == 0 ))
    then
      echo "Checking with qstat to see if there are still running/pending jobs, or if they crashed..."
      echo -e "$I [$datestamp] There are $(getNumJobsOnCE) jobs running/pending on the cluster"
      if (( $(getNumJobsOnCE) == 0 ))
      then
        echo -e "No jobs are running or pending in the grid, I guess some jobs failed!"
        if (( $(getRemainingJobs $theSuffix $folder $totalJobs) == 0 ))
        then
          echo "No, it's fine."
        else
          echo -e "$E Some jobs have failed! Exiting"
          send_mail
          exit 3
        fi
      fi
    fi

    sleep $SLEEP_TIME
    sleptTime=$((sleptTime + 1))
  done
  echo "All jobs are done, launch step 3 !" 
  cp OUTPUTS/$theSuffix/prepare_tmp.sh OUTPUTS/$theSuffix/step3_tmp.sh
  echo "echo \"y\" | sh launchAnalysis.sh 3 $analysisType" >> OUTPUTS/$theSuffix/step3_tmp.sh
 
  retryCounter=0
  while [ $retryCounter -lt 3 ]; do
    if qsub -q express -l walltime=00:30:00 -j oe -o OUTPUTS/$theSuffix/ OUTPUTS/$theSuffix/step3_tmp.sh 2>&1 | grep -q 'qsub'; then
      echo -e "$W Failed to submit to the grid, retry in 30s"
      retryCounter=$((retryCounter+1))
      sleep 30
    else 
      echo -e "$I Step 3 submitted"
      retryCounter=10
    fi
  done
  if [ $retryCounter  == 3 ]; then
    echo -e "$E Failed 3 times to send jobs, exiting"
    send_mail
    return 0
  fi
  sleep 60
  if [ $(qstat -u $USER |grep $USER|grep step3|wc -l) -gt 0 ]; then
    while [ $(getNumJobsOnCE) -gt 0 ]
    do
      datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
      echo -e "$I [$datestamp] There are $(getNumJobsOnCE) jobs remaining for step 3"
      sleep $SLEEP_TIME
    done
  fi
  echo -e "$I Step 3 done."
 
  publish_plots $theSuffix #comment this line if you don't want to publish the plots
  send_mail

  rm -f OUTPUTS/$theSuffix/prepare_tmp.sh
  rm -f OUTPUTS/$theSuffix/step2_tmp.sh
  rm -f OUTPUTS/$theSuffix/step3_tmp.sh
  
}


# Printing help if argument looks like it
onlyPublishPlots=0
for arg in "$@"
do
  case $arg in -h|-help|--help) usage  ; exit 0 ;; esac
  case $arg in -p|-publish|--publish) onlyPublishPlots=1  ;; esac #if option publish is on, just publish plots
done



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

analysisType=$1
localCopy=$2
express=$3

if [[ (-z "$analysisType") || ("$analysisType" == "-p") || ("$analysisType" == "-publish") || ("$analysisType" == "--publish") ]]; then analysisType="HZZanalysis"; fi
if [ -z "$localCopy" ]; then localCopy="0"; fi
if [ -z "$express" ]; then express="0"; fi

if ! [ "$analysisType" == "HZZanalysis" ] && ! [ "$analysisType" == "InstrMET" ] && ! [ "$analysisType" == "TnP" ]
then
  echo "$analysisType is not a known analysis"
  exit 0
fi

#Find the suffix according to the analysis type
suffixType="suffix"
if [ "$analysisType" == "HZZanalysis" ]; then
  suffixType="suffix"
elif [ "$analysisType" == "InstrMET" ]; then
  suffixType="suffix_InstrMET"
elif [ "$analysisType" == "TnP" ]; then
  suffixType="suffix_TnP"
fi
theSuffix=$(grep -oP "(?<=${suffixType}\=\").*" launchAnalysis.sh | tr -d '"')

#Possibility to publish_plots only
if [ "$onlyPublishPlots" == 1 ]; then
  publish_plots $theSuffix 1
  exit 0
fi

#Express option
queue="on the $GREEN localgrid $DEF queue"
if [ "$express" == "0" ]; then queue="on the $GREEN localgrid $DEF queue";
elif [ "$express" == "1" ]; then queue="on the $GREEN express $DEF queue";
else
  echo "$express is not a valid option for express (should be 0 or 1)"
  exit 0
fi

#Local copy option
localCopyText="to $MAG copy $DEF the data $MAG locally $DEF on the running node"
if [ "$localCopy" == "0" ]; then localCopyText="to $MAG copy $DEF the data $MAG locally $DEF on the running node";
elif [ "$localCopy" == "1" ]; then localCopyText="to read the data in $MAG streaming $DEF on the running node";
else
  echo "$localCopy is not a valid option for localCopy (should be 0 or 1)"
  exit 0
fi

#Main script
echo -e "$I Please perform some $RED tests $DEF before using this script! (PS: for $YEL help just do -h $DEF)\n"
echo -e "$W Do you wish to launch $queue the $RED FULL '$analysisType' $DEF analysis? [N/y]\n You have also asked for the jobs ${localCopyText}."
read answer
if [[ $answer == "y" ]];
then
  #if [ "$CMSSW_BASE" == "" ]; then
  #  echo -e "$W Setting CMSSW environment here (if you don't want to see this all the time, either source the script or do a cmsenv !"
  #  eval `scramv1 runtime -sh`
  #  echo "Done!"
  #fi
  #To be on he safe side, always start by setting the cmsenv. Even if already set, this ensures that it is set at the right place!
  eval `scramv1 runtime -sh`
  #Create directories
  mkdir -p OUTPUTS

  datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
  logFile="${CMSSW_BASE}/src/shears/HZZ2l2nu/OUTPUTS/fullAnalysis_${theSuffix}.${datestamp}.log"
  echo -e "$I Script launched! The log are available here: $YEL ${logFile} $DEF"
  echo -e "$I Open it with 'tail -f' for realtime update or with 'less -R' to benefit from the colour output."
  main &> $logFile &
    

fi

