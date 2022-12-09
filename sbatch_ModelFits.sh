#!/bin/bash
#
#SBATCH --job-name=BigModels
#SBATCH --time=6-00:00:00
#SBATCH -p leanew1  # Queue names you can submit to
#SBATCH -n 1
#SBATCH --mem=80G
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
module load R/4.2.0
Rscript ModelFits.R
