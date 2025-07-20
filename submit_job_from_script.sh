#!/bin/bash

#SBATCH --job-name=R-simulation       # Job name
#SBATCH --cpus-per-task=20            # CPUs
#SBATCH --mem=256gb                   # Job memory request
#SBATCH --time=48:0:0                 # Time limit
#SBATCH --output=slurm-%j.out

#Rscript -e "rmarkdown::render('~/carlos-code/cervix_cancer_microsimulation--20240617T075924Z-001/cervix_cancer_microsimulation/Cervix_MicroSim_RMarkdown_v.071.Rmd')"


module load apps/R

## Generate a unique tag (timestamp + random ID)
#unique_tag=$(date +%Y%m%d_%H%M%S)_$RANDOM

# Call your R script and pass the tag
#Rscript PARALLEL_cervix_microSim_stacked_list_v.01_original_20241010.R "$unique_tag"

# Pass the SLURM job ID to the R script:
#Rscript PARALLEL_cervix_microSim_stacked_list_v.01_original_20241010.R "$SLURM_JOB_ID"
#Rscript cervix_microSim_ARRAY_stacked_list_v.02_original_20250514.R "$SLURM_JOB_ID"
#Rscript cervix_microSim_stacked_list_v.02_original_20250514.R "$SLURM_JOB_ID"
Rscript cervix_microSim_stacked_list_v.02_original_20250514_FOR_SUBMISSION.R "$SLURM_JOB_ID"
#Rscript cervix_microSim_stacked_list_v.02_B_original_20250514.R "$SLURM_JOB_ID"
#Rscript PARALLEL_cervix_microSim_stacked_Natural_History_REPRODUCIBLE_20250328.R
#Rscript PARALLEL_cervix_microSim_stacked_Natural_History_20250211.R # for natural history simulation
