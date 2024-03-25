#!/bin/bash
#
#SBATCH --job-name=Fig1Fin
#SBATCH --time=7-00:00:00
#SBATCH -p leanew1  # Queue names you can submit to
#SBATCH -n 1
#SBATCH --mem=50G
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
module load R/4.1

<<<<<<< Updated upstream
Rscript Fig1_Boots.R
=======
Rscript Fig1_Boots_Poly_sRE.R
>>>>>>> Stashed changes
