#!/bin/bash

#SBATCH --job-name=R-simulation       # Job name
#SBATCH --cpus-per-task=8             # CPUs
#SBATCH --mem=128gb                   # Job memory request
#SBATCH --time=48:0:0                 # Time limit

#Rscript -e "rmarkdown::render('~/carlos-code/cervix_cancer_microsimulation--20240617T075924Z-001/cervix_cancer_microsimulation/Cervix_MicroSim_RMarkdown_v.071.Rmd')"


#Rscript /home/07075107P/carlos-code/cervix_cancer_microsimulation--20240617T075924Z-001/cervix_cancer_microsimulation/cervix_microSim_stacked_list.R      # Run your R script
Rscript /home/07075107P/carlos-code/cervix_cancer_microsimulation--20240617T075924Z-001/cervix_cancer_microsimulation/cervix_microSim_stacked_list_v.01.R    # Run your R script
