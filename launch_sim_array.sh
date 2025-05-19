#!/bin/bash
#SBATCH --job-name=R-jobarray
#SBATCH --mail-type=ALL
#SBATCH --mail-user=netid@illinois.edu
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5gb

#SBATCH --output=slurm_%A-%a.out
#SBATCH --array=1-5

# Fail fast on any error, and echo every command
set -e
set -x

echo "## SLURM env:"
echo "   SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}"
echo "   SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "   SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}"

TASK_ID=${SLURM_ARRAY_TASK_ID}
OUTPUT_DIR="${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}"
echo "Creating output directory: ${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Path to your params file
PARAM_FILE="${HOME}/microSim/cervix_cancer_microsimulation/data/params_job_array/inputs.txt"

# Grab the TASK_ID'th line, echo it, split into p1,p2,p3
LINE=$(sed -n "${TASK_ID}p" "${PARAM_FILE}")
echo "Line ${TASK_ID} from inputs.txt: '${LINE}'"
read -r p1 p2 p3 <<< "${LINE}"
echo "Parsed p1='${p1}'  p2='${p2}'  p3='${p3}'"

# Load R and run, passing each coverage as its own arg
module load apps/R

Rscript "$HOME/microSim/cervix_cancer_microsimulation/cervix_microSim_ARRAY_stacked_list_v.02_original_20250514.R" \
    --args "$p1" "$p2" "$p3" \
    > "${OUTPUT_DIR}/output_${TASK_ID}.txt"

# Move the SLURM log into the same folder
mv "${SLURM_SUBMIT_DIR}/slurm_${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}.out" \
   "${OUTPUT_DIR}/"
