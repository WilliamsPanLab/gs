#!/bin/bash
#
#SBATCH --job-name=BigModels
#SBATCH --time=7-00:00:00
#SBATCH -p leanew1  # Queue names you can submit to
#SBATCH -n 1
#SBATCH --mem=90G
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
module load R/4.1
#Rscript gp_DevExplainedBoots.R 
#Rscript gExt_Fits.R
#Rscript gInt_Fits.R
#Rscript gp_Fits.R
Rscript gp_gausVnb_qq.R
