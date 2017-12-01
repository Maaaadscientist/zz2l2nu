#!/bin/bash

SLEEP_TIME_QSTAT=900 # in seconds. 6min is the min
SLEEP_TIME=60 #in seconds

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
  sharedSuffix=$1
  folder_output=$2
  totalJobs=$3
  jobsDone=$(ls -1 ${CMSSW_BASE}/src/shears/HZZ2l2nu/${folder_output}_${sharedSuffix} | wc -l) #exportedSuffix is an exported variable from the script launchAnalysis.sh
  echo $(($totalJobs-$jobsDone)) #number of running/remaining jobs
}

function publish_plots(){
  sharedSuffix=$1
  datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
  echo -e "$I Creating symbolic link to your public_html folder..."
  if [ $(ls -1 ${CMSSW_BASE}/src/shears/HZZ2l2nu/plots_${sharedSuffix}/ |wc -l) -eq 0 ]; then
    echo -e "$E No plots to publish"
    exit 4
  else
    mkdir -p ~/public_html
    chmod 755 ~/public_html
    mkdir -p ~/public_html/SHEARS_PLOTS
    mkdir -p ~/public_html/SHEARS_PLOTS/plots_${sharedSuffix}
    ln -s ${CMSSW_BASE}/src/shears/HZZ2l2nu/plots_${sharedSuffix}/* ~/public_html/SHEARS_PLOTS/plots_${sharedSuffix}/.
    cp ${CMSSW_BASE}/src/shears/HZZ2l2nu/Tools/index.php ~/public_html/SHEARS_PLOTS/plots_${sharedSuffix}/.
    echo -e "$I Your plots are available in ~/public_html/SHEARS_PLOTS/plots_${sharedSuffix}/, i.e. on http://homepage.iihe.ac.be/~$USER/SHEARS_PLOTS/plots_${sharedSuffix}/"
  fi
}

function send_mail(){
  export LD_LIBRARY_PATH="/lib64/:$LD_LIBRARY_PATH"
  mailAddress=$(ldapsearch -LLL -x uid=$USER mail | sed -n 's/^[ \t]*mail:[\t]*\(.*\)/\1/p')
  datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
  sed -i "1s/^/Subject: Jobs for SHEARS $datestamp\n/" $logFile
  sendmail $mailAddress < $logFile
}

function main(){
  rm -f tmp_shared_variables.txt
  touch tmp_shared_variables.txt #create the shared variables file to interact with launchAnalysis.sh to share the suffix it uses
  rm -f prepare_tmp.sh
  rm -f step2_tmp.sh
  rm -f step3_tmp.sh
  
  #0) Full cleaning
  datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
  echo $datestamp
  echo "Starting full cleaning..."
  echo "a" | source launchAnalysis.sh 0
  
  #1) Launch analysis on cluster
  echo "Starting step 1..."
  yes | source launchAnalysis.sh 1 $analysisType $localCopy $express
  if [ $? -eq 0 ]; then
    echo -e "$E Step 1 failed, exiting."
    return 0
  fi
  sleep 60

  #2) Harvest on express queue
  echo "Waiting for step 1 to be over..."
  #a) prepare environment
  echo "source \$VO_CMS_SW_DIR/cmsset_default.sh" >> prepare_tmp.sh
  echo "export SCRAM_ARCH=slc6_amd64_gcc530" >> prepare_tmp.sh
  echo "export INITDIR=$CMSSW_BASE/src/shears/HZZ2l2nu" >> prepare_tmp.sh
  echo "cd \$INITDIR" >> prepare_tmp.sh
  echo "hostname ;" >> prepare_tmp.sh
  echo "date;" >> prepare_tmp.sh
  #b) check number of jobs
  while read -r line
  do
    sharedSuffix="$line"
  done < tmp_shared_variables.txt
  while [ $(wc -l < sendJobs_${sharedSuffix}.cmd) -gt 0 ] #this is here in case a job failed to be sent... and so we wait for him.
  do
    datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
    echo -e "$I [$datestamp] There are still $(wc -l < sendJobs_${sharedSuffix}.cmd) jobs to send"
    sleep 60
  done
  if [ $(grep -c -e '^qsub ' $(ls -Art big-submission-*.err | tail -n 1) ) -gt 0 ]; then
    echo -e "$E There are jobs not submitted. To be safe, let's stop here"
    return 0
  fi
  folder="OUTPUTS"
  totalJobs=$(ls -1 ${CMSSW_BASE}/src/shears/HZZ2l2nu/JOBS | wc -l)
  sleptTime=1  #don't make it start at 0
  while [ $(getRemainingJobs $sharedSuffix $folder $totalJobs) -gt 0 ]
  do
    datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
    echo -e "$I [$datestamp] There are $(getRemainingJobs $sharedSuffix $folder $totalJobs) jobs remaining" 

    if (( $sleptTime % ($SLEEP_TIME_QSTAT/$SLEEP_TIME) == 0 ))
    then
      echo "Checking with qstat to see if there are still running/pending jobs, or if they crashed..."
      echo -e "$I [$datestamp] There are $(getNumJobsOnCE) jobs running/pending on the cluster"
      if (( $(getNumJobsOnCE) == 0 ))
      then
        echo -e "No jobs are running or pending in the grid, I guess some jobs failed!"
        if (( $(getRemainingJobs $sharedSuffix $folder $totalJobs) == 0 ))
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
  cp prepare_tmp.sh step2_tmp.sh
  echo "echo \"y\" | sh launchAnalysis.sh 2 $analysisType" >> step2_tmp.sh
  qsub -q express -l walltime=00:30:00 -j oe step2_tmp.sh 
  sleep 60

  #3) Do data-MC comparison
  echo "Waiting for step 2 to be over..." 
  folder="merged"
  totalJobs=$(ls -1 ${CMSSW_BASE}/src/shears/HZZ2l2nu/OUTPUTS_${sharedSuffix} |grep _0.root | wc -l)
  sleptTime=1 #don't make it start at 0
  while [ $(getRemainingJobs $sharedSuffix $folder $totalJobs) -gt 0 ]
  do
    datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
    echo -e "$I [$datestamp] There are $(getRemainingJobs $sharedSuffix $folder $totalJobs) jobs remaining" 

    if (( $sleptTime % ($SLEEP_TIME_QSTAT/$SLEEP_TIME) == 0 ))
    then
      echo "Checking with qstat to see if there are still running/pending jobs, or if they crashed..."
      echo -e "$I [$datestamp] There are $(getNumJobsOnCE) jobs running/pending on the cluster"
      if (( $(getNumJobsOnCE) == 0 ))
      then
        echo -e "No jobs are running or pending in the grid, I guess some jobs failed!"
        if (( $(getRemainingJobs $sharedSuffix $folder $totalJobs) == 0 ))
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
  echo "All jobs are done, launch step 3 !" 
  cp prepare_tmp.sh step3_tmp.sh
  echo "echo \"y\" | sh launchAnalysis.sh 3 $analysisType" >> step3_tmp.sh
  qsub -q express -l walltime=00:30:00 -j oe step3_tmp.sh 
  sleep 60
  if [ $(qstat -u $USER |grep $USER|wc -l) -gt 0 ]; then
    while [ $(getNumJobsOnCE) -gt 0 ]
    do
      echo -e "$I [$datestamp] There are $(getNumJobsOnCE) jobs remaining for step 3"
      sleep $SLEEP_TIME
    done
  fi
  echo -e "$I Step 3 done."
 
  publish_plots $sharedSuffix #comment this line if you don't want to publish the plots
  send_mail #Comment this line if you don't like to receive emails when jobs done :)

  rm -f tmp_shared_variables.txt
  rm -f prepare_tmp.sh
  rm -f step2_tmp.sh
  rm -f step3_tmp.sh
  
}


# Printing help if argument looks like it
case $1 in -h|-help|--help) usage  ; exit 0 ;; esac

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

if [ -z "$analysisType" ]; then analysisType="HZZanalysis"; fi
if [ -z "$localCopy" ]; then localCopy="0"; fi
if [ -z "$express" ]; then express="0"; fi

if ! [ "$analysisType" == "HZZanalysis" ] && ! [ "$analysisType" == "InstrMET" ]
then
  echo "$analysisType is not a known analysis"
  exit 0
fi
queue="on the $GREEN localgrid $DEF queue"
if [ "$express" == "0" ]; then queue="on the $GREEN localgrid $DEF queue";
elif [ "$express" == "1" ]; then queue="on the $GREEN express $DEF queue";
else
  echo "$express is not a valid option for express (should be 0 or 1)"
  exit 0
fi

localCopyText="to $MAG copy $DEF the data $MAG locally $DEF on the running node"
if [ "$localCopy" == "0" ]; then localCopyText="to $MAG copy $DEF the data $MAG locally $DEF on the running node";
elif [ "$localCopy" == "1" ]; then localCopyText="to read the data in $MAG streaming $DEF on the running node";
else
  echo "$localCopy is not a valid option for localCopy (should be 0 or 1)"
  exit 0
fi

echo -e "$I Please perform some $RED tests $DEF before using this script! (PS: for $YEL help just do -h $DEF)\n"
echo -e "$W Do you wish to launch $queue the $RED FULL '$analysisType' $DEF analysis? [N/y]\n You have also asked for the jobs ${localCopyText}."
read answer
if [[ $answer == "y" ]];
then
  if [ "$CMSSW_BASE" == "" ]; then
    echo -e "$W Setting CMSSW environment here (if you don't want to see this all the time, either source the script or do a cmsenv !"
    eval `scramv1 runtime -sh`
    echo "Done!"
  fi
  datestamp=$(date  +%Y-%m-%d-%H:%M:%S)
  logFile="${CMSSW_BASE}/src/shears/HZZ2l2nu/fullAnalysis.$datestamp.log"
  echo -e "$I Script launched! The log are available here: $YEL fullAnalysis.${datestamp}.log $DEF"
  echo -e "$I Open it with 'tail -f' for realtime update or with 'less -R' to benefit from the colour output."
  main &> $logFile &
    

fi

