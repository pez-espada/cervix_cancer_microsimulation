#!/bin/bash

## Describe requirements for computing ----
#SBATCH --job-name=R-jobarray
#SBATCH --mail-type=ALL
#SBATCH --mail-user=netid@illinois.edu
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5gb

## Standard output and error log (to a generic file in submission dir)
#SBATCH --output=slurm_%A-%a.out

# Array range
#SBATCH --array=1-5

## Setup computing environment for job ----
JOB_ID="${SLURM_ARRAY_JOB_ID}"
TASK_ID="${SLURM_ARRAY_TASK_ID}"
OUTPUT_DIR="${SLURM_SUBMIT_DIR}/${JOB_ID}"
SLURM_OUTPUT_FILE="slurm_${SLURM_JOB_ID}-${TASK_ID}.out"
TARGET_SLURM_OUTPUT="${OUTPUT_DIR}/${SLURM_OUTPUT_FILE}"

## Create a directory for the data output based on the SLURM_ARRAY_JOB_ID
mkdir -p "${OUTPUT_DIR}"

### Switch directory into job ID (puts all output here)
#cd "${OUTPUT_DIR}"

## Run simulation ----
module load apps/R

#export PARAMS=`cat ${HOME}/microSim/cervix_cancer_microsimulation/data/params_job_array/inputs.txt | sed -n "${TASK_ID}"p`

# Read the TASK_ID’th line of your inputs file and split it directly into p1 p2 p3
read -r p1 p2 p3 < <(sed -n "${TASK_ID}p" "${HOME}/microSim/cervix_cancer_microsimulation/data/params_job_array/inputs.txt")


Rscript "$HOME/microSim/cervix_cancer_microsimulation/cervix_microSim_ARRAY_stacked_list_v.02_original_20250514.R" \
    --args "$p1" "$p2" "$p3" \
    > "${OUTPUT_DIR}/output_${TASK_ID}.txt"

#Rscript "$HOME/microSim/cervix_cancer_microsimulation/cervix_microSim_stacked_list_v.02_original_20250514.R" --args "$PARAMS" > "output_${TASK_ID}.txt"
#Rscript "$HOME/microSim/cervix_cancer_microsimulation/cervix_microSim_ARRAY_stacked_list_v.02_original_20250514.R" --args "$PARAMS" > "output_${TASK_ID}.txt"
#Rscript "$HOME/microSim/cervix_cancer_microsimulation/cervix_microSim_ARRAY_stacked_list_v.02_original_20250514.R" --args "$PARAMS" > "${OUTPUT_DIR}/output_${TASK_ID}.txt"

## Move the Slurm output file to the job directory
mv "${SLURM_SUBMIT_DIR}/${SLURM_OUTPUT_FILE}" "${TARGET_SLURM_OUTPUT}"
