#!/bin/bash
#
#SBATCH --job-name=Fig3_Boots1k.R
#SBATCH --time=1-10:00:00
#SBATCH -p normal,leanew1  # Queue names you can submit to
#SBATCH -n 1
#SBATCH --mem=70G
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
module load R/4.1

Rscript Fig3_Boots1k.R
