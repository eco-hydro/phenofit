#!/bin/bash

#SBATCH --array=0-3
#SBATCH --job-name=phenofit
#SBATCH --error=log/%x_err_%a.err
#SBATCH --output=log/%x_out_%a.out
#SBATCH --share
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1Gb
#
#SBATCH --mail-type=END
#SBATCH --mail-user=Dongdong.Kong@csiro.au

module load R
# Rscript -e 'source("test/GEE/02_whit_lambda_main.R")'
Rscript -e 'source("test/07_whit/01_phenofit_main_test_flux&cam.R")'
