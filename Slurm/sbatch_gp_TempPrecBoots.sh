#!/bin/bash
#
#SBATCH --job-name=TmpPrec
#SBATCH --time=7-00:00:00
#SBATCH -p leanew1  # Queue names you can submit to
#SBATCH -n 1
#SBATCH --mem=110G
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
module load R/4.1

#Rscript gp_TempPrecBoots.R

echo intermission

Rscript pg_TemPrecBoots.R
