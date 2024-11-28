#!/bin/bash

#SBATCH --job-name=R-simulation       # Job name
#SBATCH --cpus-per-task=20            # CPUs
#SBATCH --mem=256gb                    # Job memory request
#SBATCH --time=48:0:0                 # Time limit
#SBATCH --output=slurm-%j.out

#Rscript -e "rmarkdown::render('~/carlos-code/cervix_cancer_microsimulation--20240617T075924Z-001/cervix_cancer_microsimulation/Cervix_MicroSim_RMarkdown_v.071.Rmd')"

#Rscript /home/07075107P/carlos-code/cervix_cancer_microsimulation--20240617T075924Z-001/cervix_cancer_microsimulation/cervix_microSim_stacked_list.R      # Run your R script
#Rscript /home/07075110P/carlos-code/cervix_cancer_microsimulation--20240617T075924Z-001/cervix_cancer_microsimulation/cervix_microSim_stacked_list_v.01.R    # Run your R script
#Rscript /home/07075110P/carlos-code/cervix_cancer_microsimulation--20240617T075924Z-001/cervix_cancer_microsimulation/parallel_cervix_microSim_stacked_list_v.02_DeBug.R   # Run your R script
#Rscript /home/07075107P/microSim/cervix_cancer_microsimulation/parallel_cervix_microSim_stacked_list_v.02_DeBug.R
#Rscript /home/07075107P/microSim/cervix_cancer_microsimulation/parallel_cervix_microSim_stacked_list_v.02_DeBug.R
#Rscript /home/07075107P/microSim/cervix_cancer_microsimulation/PARALLEL_cervix_microSim_stacked_list_v.01_original_20241010.R

module load apps/R
Rscript PARALLEL_cervix_microSim_stacked_list_v.01_original_20241010.R
