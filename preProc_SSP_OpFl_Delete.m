function preProc_SSP_Delete(subj)
%%% This function will take a single subject's NDAR name, download their fMRI data and motion masks, concatenate the fMRI data, mask according to Robert's instructions (.2mm FD, power outliers out), derive a single-subject parcellation based on Pines et al. 2022's group templates, and delete the input fMRI data.

% print subject being ran
subj

% add matlab path for used functions
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));

% tell it where AWS tools are for downloads
%system('export PATH=/cbica/projects/abcdfnets/aws/dist/:$PATH')

% echo subj name into a .txt
%subjTxtCommand=['echo ' subj ' >> /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/' subj '.txt'];
%system(subjTxtCommand)

% download that one subject's data

% ∆∆∆∆∆∆
% keep looking at github for updates on authent. fix
% ∆∆∆∆∆∆

% subjDlCommand=['python3 /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/download.py -i /scratch/users/apines/datastructure_manifest.txt -o /scratch/users/apines/ -s /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/' subj '.txt -l /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/dl_logs -d /oak/stanford/groups/leanew1/users/apines/data_subsets_3_9_21_from9620.txt &']

% note: downloader tool does not seem to communicate when it is done to matlab
% added '&' and 'pause' so that matlab waits 5 minutes to proceed rather than getting caught up indefinitely
system(subjDlCommand)
pause(300)



% now the matlab portions. Apply the motion mask to the downloaded data
%%% probably will lose psychopathology instances if we do this
%%% apply_motion_mask(subj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ∆∆∆∆ EITHER RUN THIS FOR EACH TASK, OR CONCATENATE PRIOR TO LOADING THIS IN ∆∆∆∆∆∆∆∆∆∆∆∆
% time series filepath
TS=
% make output dir
mkdir(['/scratch/users/apines/' subj '/subcort.ptseries.nii'])
% cifti-parcellate to get subcortical time series
ParcCommand=['wb_command -cifti-parcellate ' TS '/oak/stanford/groups/leanew1/users/apines/maps/Tian_Subcortex_S2_3T_32k.dlabel.nii COLUMN /scratch/users/apines/' subj '/subcort.ptseries.nii'];
system(ParcCommand)

% downsample to fsaverage5
subjDSCommand=['DS_surf_ts_fs5.sh ' subj ' &']
system(subjDSCommand)
pause(45)

% this will become concatenate only as download works
% concatenate masked time series and isolate the cortex (for cortical surface only SSP)
%concat_TS_and_IsoCort(subj)

% derive an indivudalized parcellation (just parent filepath for inputTS_fp, g_ls will find 10k mgh extension within folder
for k=2:30
	PersonalizeNetworks(inputTS_fp,k,subj)
end

% for "cifti_read", interferes if added earlier
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/cifti-matlab/'));

% calculate multiscale FC

% saveout low memory version

% run optical flow pipeline ∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆
%%% calculate optical flow
OpFl_abcd(tsIn_L,tsIn_R,tsOut)
% set output filepaths to calculcate angular distances
rsOut=${childfp}/${subj}_${sesh}_PGGDist_rs_fs5.mat

% needed?
mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj}
% calculate angular distances (change output filepath to something meaningful)
AngDCalcCmd=['/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/run_Extract_BUTD_ResultantVecs_Gran_fs5.sh /share/software/user/restricted/matlab/R2018a/ ' tsOut '$childfp/OpFl_timeseries_L_fs5.mat $childfp/OpFl_timeseries_R_fs5.mat'];
system(AngDCalcCmd)

% saveout low memory version

% delete input data
Delete_input_data(subj)
