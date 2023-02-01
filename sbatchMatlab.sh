#!/bin/bash
#
#SBATCH --job-name=matlab
#SBATCH --time=1-00:00:00
#SBATCH -p leanew1,normal  # Queue names you can submit to
#SBATCH -n 1
#SBATCH --mem=40G
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
ml aws-cli
ml python3
matlab -nodisplay -r "OpFl_abcd('${1}','${2}','${3}')"
