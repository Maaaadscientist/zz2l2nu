#!/bin/bash
#SBATCH --job-name=parallel_job
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null


# Directory containing .sh scripts
SCRIPT_DIR="jobs/scripts"
LOG_DIR="jobs/logs"

# Ensure the log directory exists
mkdir -p $LOG_DIR

# Execute scripts in parallel
for SCRIPT in ${SCRIPT_DIR}/*.sh; do
    # Extract the script name without path for use in log filenames
    SCRIPT_NAME=$(basename $SCRIPT)
    # Redirect output and errors of each script to separate log files
    bash $SCRIPT > "${LOG_DIR}/${SCRIPT_NAME}_$$.out" &
    
    # Manage job submission to limit to the number of ntasks
    # Ensure that the number of background jobs does not exceed the number of ntasks
    while [ $(jobs -p | wc -l) -ge $SLURM_NTASKS ]; do
        sleep 1
        wait -n
    done
done

# Wait for all remaining background jobs to complete
wait
