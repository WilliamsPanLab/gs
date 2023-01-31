subj=$1
sesh=ses-baselineYear1Arm1
task1=task-rest
# can add in others once dl working

# set freesurfer dir
export FREESURFER_HOME=/share/software/user/open/freesurfer/6.0.0
# set subjs dir
export SUBJECTS_DIR=/share/software/user/open/freesurfer/6.0.0/subjects
# set freesurfer license
export FS_LICENSE=/oak/stanford/groups/leanew1/users/apines/license.txt
# file fp
parentfp=/scratch/users/apines/derivatives/abcd-hcp-pipeline/${subj}/${sesh}/func
childfp=/scratch/users/apines/derivatives/abcd-hcp-pipeline/${subj}/${sesh}/func
# make directory if it does not exist
mkdir ${childfp}

# subject's aggregated time series
AgTS=${childfp}/${subj}_${sesh}_${task1}_bold_desc-filtered_timeseries.dtseries.nii

# separate hemispheres - left
wb_command -cifti-separate $AgTS COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_${sesh}_L_AggTS.func.gii 

# right hemi
wb_command -cifti-separate $AgTS COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_${sesh}_R_AggTS.func.gii

### resample both hemis to 10k vertices
# left hemisphere
wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_AggTS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_AggTS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii

# right hemisphere
wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_AggTS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_AggTS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii

# convert to mgh for reading individual hemisphere time series into matlab
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_AggTS_10k.func.gii ${childfp}/${subj}_${sesh}_L_AggTS_10k.mgh
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_AggTS_10k.func.gii ${childfp}/${subj}_${sesh}_R_AggTS_10k.mgh


