# record input argument as subj
subj=$1

# load modules needed
ml biology
ml workbench

# set rsfrmi fp
rsfp=/scratch/users/apines/abcd_images/fmriresults01/derivatives/abcd-hcp-pipeline/${subj}/ses-baselineYear1Arm1/func/${subj}_ses-baselineYear1Arm1_task-rest_p2mm_masked.dtseries.nii

# circuit ROI fps
# DMN
D1=~/2021-masks/Medial_amPFC_DefaultModeNetwork_n2_50_n6.nii.gz.shape.gii
D2=~/2021-masks/Left_AG_DefaultModeNetwork_n46_n70_32.nii.gz.shape.gii
D3=~/2021-masks/Right_AG_DefaultModeNetwork_50_n62_26.nii.gz.shape.gii
D4=~/2021-masks/Medial_PCC_DefaultModeNetwork_0_n50_28.nii.gz.shape.gii
# Salience
S1=~/2021-masks/Left_antInsula_Salience_n38_14_n6.nii.gz.shape.gii
S2=~/2021-masks/Right_antInsula_Salience_38_18_2.nii.gz.shape.gii
####
# use coarse Tian to extract amyg
####
# Attention
A1=~/2021-masks/Medial_msPFC_Attention_n2_14_52.nii.gz.shape.gii
A2=~/2021-masks/Left_lPFC_Attention_n44_6_32.nii.gz.shape.gii
A3=~/2021-masks/Right_lPFC_Attention_50_10_28.nii.gz.shape.gii
A4=~/2021-masks/Left_aIPL_Attention_n30_n54_40.nii.gz.shape.gii
A5=~/2021-masks/Right_aIPL_Attention_38_n56_48.nii.gz.shape.gii
A6=~/2021-masks/Left_precuneus_Attention_n14_n66_52.nii.gz.shape.gii
A7=~/2021-masks/Right_precuneus_Attention_18_n68_52.nii.gz.shape.gii


# convert alff to pt series
wb_command -cifti-parcellate ${alff_rs1xcpSubcort_fp} ${scale3fp} COLUMN ${alff_p_SubcortTS1}
wb_command -cifti-parcellate ${alff_rs2xcpSubcort_fp} ${scale3fp} COLUMN ${alff_p_SubcortTS2}
wb_command -cifti-parcellate ${alff_wmxcpSubcort_fp} ${scale3fp} COLUMN ${alff_p_SubcortTS3}
wb_command -cifti-parcellate ${alff_gambxcpSubcort_fp} ${scale3fp} COLUMN ${alff_p_SubcortTS5}
wb_command -cifti-parcellate ${alff_emoxcpSubcort_fp} ${scale3fp} COLUMN ${alff_p_SubcortTS6}

# convert ptseries to text
wb_command -cifti-convert -to-text $alff_p_SubcortTS1 $alff_SubcortTS1
wb_command -cifti-convert -to-text $alff_p_SubcortTS2 $alff_SubcortTS2
wb_command -cifti-convert -to-text $alff_p_SubcortTS3 $alff_SubcortTS3
wb_command -cifti-convert -to-text $alff_p_SubcortTS5 $alff_SubcortTS5
wb_command -cifti-convert -to-text $alff_p_SubcortTS6 $alff_SubcortTS6

