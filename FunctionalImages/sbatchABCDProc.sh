#!/bin/bash
#
#SBATCH --job-name=ABCDParticipant
#SBATCH --time=6:00:00
#SBATCH -p normal,leanew1  # Queue names you can submit to
#SBATCH -n 1
#SBATCH --mem=30G
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
ml biology
ml matlab
ml workbench
ml freesurfer
ml python/3
# run matlab pipeline
matlab -nodisplay -r "preProc_SSP_OpFl_Delete('${1}')"
