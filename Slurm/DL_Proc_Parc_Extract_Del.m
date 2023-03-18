function preProc_SSP_Delete(subj)
%%% This function will take a single subject's NDAR name, download their fMRI data and motion masks, concatenate the fMRI data, mask according to Robert's instructions (.2mm FD, power outliers out), derive a single-subject parcellation based on Pines et al. 2022's group template, and delete the input fMRI data.

% print subject being ran
subj

% add matlab path for used functions
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));

% echo subj name into a .txt
subjTxtCommand=['echo ' subj ' >> /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/' subj '.txt'];
system(subjTxtCommand)

% download that one subject's data
subjDlCommand=['python3 /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/download.py -dp 1210784 -m ~/datastructure_manifest.txt -o /scratch/users/apines/dl_folder/dl_folder/ -s /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/' subj '.txt -b /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/data_subsets_smr.txt -l /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/log/ &'];

% note: downloader tool does not seem to communicate when it is done to matlab
% added '&' and 'pause' so that matlab waits 5 minutes to proceed rather than getting caught up indefinitely
system(subjDlCommand)
pause(230)

% now the matlab portions. Apply the motion mask to the downloaded data
apply_motion_mask(subj)

% concatenate masked time series and isolate the cortex (for cortical surface only SSP)
concat_TS(subj)

% derive an indivudalized parcellation
Individualize_ciftiSurf_resampledGroCon(subj)

% for "cifti_read", interferes if added earlier
addpath(genpath('/cbica/projects/abcdfnets/scripts/cifti-matlab/'));

% convert to dscalar hard parcel
mat_to_dlabel(subj)

% delete input data
Delete_input_data(subj)
