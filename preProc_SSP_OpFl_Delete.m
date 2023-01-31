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

% keep looking at github for updates on authent. fix

% subjDlCommand=['python3 /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/download.py -i /scratch/users/apines/datastructure_manifest.txt -o /scratch/users/apines/ -s /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/' subj '.txt -l /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/dl_logs -d /oak/stanford/groups/leanew1/users/apines/data_subsets_3_9_21_from9620.txt &']

% note: downloader tool does not seem to communicate when it is done to matlab
% added '&' and 'pause' so that matlab waits 5 minutes to proceed rather than getting caught up indefinitely
system(subjDlCommand)
pause(300)

% now the matlab portions. Apply the motion mask to the downloaded data
apply_motion_mask(subj)

% consider downsampling
subjDSCommand=['DS_surf_ts_fs5.sh ' subj ' &']
system(subjDSCommand)
pause(45)

% concatenate masked time series and isolate the cortex (for cortical surface only SSP)
concat_TS_and_IsoCort(subj)

% derive an indivudalized parcellation
Individualize_ciftiSurf_resampledGroCon(subj)

% for "cifti_read", interferes if added earlier
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/cifti-matlab/'));

% convert to dscalar hard parcel
mat_to_dlabel(subj)

% calculate multiscale FC
% saveout low memory version

% run optical flow pipeline
%
%
%
%
% saveout low memory version

% delete input data
Delete_input_data(subj)
