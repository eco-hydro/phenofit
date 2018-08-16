#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH --array=0-5
#SBATCH --job-name=phenofit
#SBATCH --error=result/log/%x_err_%a.err
#SBATCH --output=result/log/%x_out_%a.out
#SBATCH --share
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=1Gb
#
#SBATCH --mail-type=END
#SBATCH --mail-user=Dongdong.Kong@csiro.au

module load R
# Rscript -e 'source("test/GEE/02_whit_lambda_main.R")'
# Rscript -e 'source("test/07_whit/01_phenofit_main_test_flux&cam.R")'
# Rscript -e 'source("test/07_whit/whit_lambda/02_whit_lambda_main.R")'
# Rscript -e 'source("test/07_whit/05_compare_with_wHANTS&wSG.R")'
Rscript -e 'source("test/07_whit/07_val_phenofit_methods.R")'


# ssh kon055@bracewell-login.hpc.csiro.au 
# sinteractive -n1 -c4 -t03:30:00 -m8GB

# # cd /OSM/CBR/CoRE/working/timeseries/Climate/kon055

# module load matlab
# module load matlab/R2015b
# matlab -nojvm -nosplash -singleCompThread -r mycode

# squeue -u kon055
# scancel -u kon055 -n flux112_
# sacct -u kon055 --job 13127327libr
