#!/bin/bash
#
#SBATCH --job-name=DevExplainedBoots
#SBATCH --time=2-00:00:00
#SBATCH -p leanew1  # Queue names you can submit to
#SBATCH -n 1
#SBATCH --mem=25G
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
module load R/4.1
Rscript DevExplBoots.R 
