# convert MNI-coordinate surface-circuit ROIs to cifti-space


# DMN
D1=~/2021-masks/Medial_amPFC_DefaultModeNetwork_n2_50_n6.nii.gz
D2=~/rois/Defaultmode_AG_L.nii
D3=~/rois/Defaultmode_AG_R.nii
D4=~/rois/Defaultmode_PCC_M.nii

# Salience
S1=~/rois/Salience_AI_L.nii
S2=~/rois/Salience_AI_R.nii
# amyg already in cifti-space

# Attention
A1=~/rois/Attention_msPFC_M.nii
A2=~/rois/Attention_LPFC_L.nii
A3=~/rois/Attention_LPFC_R.nii
A4=~/rois/Attention_aIPL_L.nii
A5=~/rois/Attention_aIPL_R.nii
A6=~/rois/Attention_PCUN_L.nii
A7=~/rois/Attention_PCUN_R.nii

# load modules for resampling
ml biology
ml workbench

### resample each
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
wb_command -volume-to-surface-mapping ${A4} ~/S1200_MSMAll3T1071.L.inflated_MSMAll.32k_fs_LR.surf.gii ${A4}.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${A5} ~/S1200_MSMAll3T1071.R.inflated_MSMAll.32k_fs_LR.surf.gii ${A5}.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${A6} ~/S1200_MSMAll3T1071.L.inflated_MSMAll.32k_fs_LR.surf.gii ${A6}.shape.gii -trilinear
wb_command -volume-to-surface-mapping ${A7} ~/S1200_MSMAll3T1071.R.inflated_MSMAll.32k_fs_LR.surf.gii ${A7}.shape.gii -trilinear

