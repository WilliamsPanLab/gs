#!/bin/bash
#
#SBATCH --job-name=Fig2Fin
#SBATCH --time=2-00:00:00
#SBATCH -p normal,leanew1  # Queue names you can submit to
#SBATCH -n 1
#SBATCH --mem=45G
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
module load R/4.1
Rscript Fig_2_povNonPov_sRE_fancy_numericPov.R
