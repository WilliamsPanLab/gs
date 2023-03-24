# convert MNI-coordinate surface-circuit ROIs to cifti-space
ml biology
ml workbench
ml ants/2.4.0

# DMN
D1=~/2021-masks/Medial_amPFC_DefaultModeNetwork_n2_50_n6.nii.gz
D2=~/2021-masks/Left_AG_DefaultModeNetwork_n46_n70_32.nii.gz
D3=~/2021-masks/Right_AG_DefaultModeNetwork_50_n62_26.nii.gz
D4=~/2021-masks/Medial_PCC_DefaultModeNetwork_0_n50_28.nii.gz

# Salience
S1=~/2021-masks/Left_antInsula_Salience_n38_14_n6.nii.gz
S2=~/2021-masks/Right_antInsula_Salience_38_18_2.nii.gz
# amyg already in cifti-space

# Attention
A1=~/2021-masks/Medial_msPFC_Attention_n2_14_52.nii.gz
A2=~/2021-masks/Left_lPFC_Attention_n44_6_32.nii.gz
A3=~/2021-masks/Right_lPFC_Attention_50_10_28.nii.gz
A4=~/2021-masks/Left_aIPL_Attention_n30_n54_40.nii.gz
A5=~/2021-masks/Right_aIPL_Attention_38_n56_48.nii.gz
A6=~/2021-masks/Left_precuneus_Attention_n14_n66_52.nii.gz
A7=~/2021-masks/Right_precuneus_Attention_18_n68_52.nii.gz

# load modules for resampling
ml biology
ml workbench

# A4 needs dilation
ImageMath 3 ${A4}_dil.nii.gz GD ${A4} 1

### resample each

# sink some more time into seeing if there is a better approach for MNI->FSLR

# DMN
wb_command -volume-to-surface-mapping ${D1} ~/S1200_MSMAll3T1071.L.inflated_MSMAll.32k_fs_LR.surf.gii ${D1}_L.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${D1} ~/S1200_MSMAll3T1071.R.inflated_MSMAll.32k_fs_LR.surf.gii ${D1}_R.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${D2} ~/S1200_MSMAll3T1071.L.inflated_MSMAll.32k_fs_LR.surf.gii ${D2}.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${D3} ~/S1200_MSMAll3T1071.R.inflated_MSMAll.32k_fs_LR.surf.gii ${D3}.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${D4} ~/S1200_MSMAll3T1071.L.inflated_MSMAll.32k_fs_LR.surf.gii ${D4}_L.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${D4} ~/S1200_MSMAll3T1071.R.inflated_MSMAll.32k_fs_LR.surf.gii ${D4}_R.shape.gii -trilinear
# total of 6 rois for L + R DMN

# Salience
wb_command -volume-to-surface-mapping ${S1} ~/S1200_MSMAll3T1071.L.inflated_MSMAll.32k_fs_LR.surf.gii ${S1}.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${S2} ~/S1200_MSMAll3T1071.R.inflated_MSMAll.32k_fs_LR.surf.gii ${S2}.shape.gii -trilinear
# amygdalae to be added straight from fslr later for total of 4

# Attention
wb_command -volume-to-surface-mapping ${A1} ~/S1200_MSMAll3T1071.L.inflated_MSMAll.32k_fs_LR.surf.gii ${A1}_L.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${A1} ~/S1200_MSMAll3T1071.R.inflated_MSMAll.32k_fs_LR.surf.gii ${A1}_R.shape.gii -trilinear

wb_command -volume-to-surface-mapping ${A2} ~/S1200_MSMAll3T1071.L.inflated_MSMAll.32k_fs_LR.surf.gii ${A2}.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${A3} ~/S1200_MSMAll3T1071.R.inflated_MSMAll.32k_fs_LR.surf.gii ${A3}.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${A4}_dil.nii.gz ~/S1200_MSMAll3T1071.L.inflated_MSMAll.32k_fs_LR.surf.gii ${A4}.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${A5} ~/S1200_MSMAll3T1071.R.inflated_MSMAll.32k_fs_LR.surf.gii ${A5}.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${A6} ~/S1200_MSMAll3T1071.L.inflated_MSMAll.32k_fs_LR.surf.gii ${A6}.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${A7} ~/S1200_MSMAll3T1071.R.inflated_MSMAll.32k_fs_LR.surf.gii ${A7}.shape.gii -trilinear

### great, now lets make them parcellations so we can extract from em
# convert shape giis to dlabels
#DMN
wb_command -label-convert -from-nifti ${A1}_L.shape.gii -to-label lh.dlabel.nii
wb_command -label-convert -from-nifti lh.shape.gii -to-label lh.dlabel.nii
# 

# concatenate shape giis



# convert dscalar subcortex to dlabel




