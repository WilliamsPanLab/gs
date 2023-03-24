function preProc_SSP_Delete(subj)

%%% CALL ml R/4.1
% ml biology
% ml workbench
% ml python/3
%%%% IN .SH CALLER SCRIPT

%%% This function will take a single subject's NDAR name, download their fMRI data and motion masks, concatenate the fMRI data, mask according to Robert's paper (.2mm FD, power outliers out), derive a single-subject parcellation based on Pines et al. 2022's group templates,  and delete the input fMRI data.

% print subject being ran
subj

% add matlab path for used functions
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/gp/FunctionalImages');

% echo subj name into a .txt
subjTxtCommand=['echo ' subj ' >> /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/' subj '.txt'];
system(subjTxtCommand)

% download this participant's data
subjDlCommand=['python3 /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/download.py -dp 1210784 -m /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/datastructure_manifest.txt -o /scratch/users/apines/abcd_images -s /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/' subj '.txt -l /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/dl_logs -b /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/data_subsets_Final.txt']

% Define the expect script
expect_script = ['~/expect_script_' subj '.exp'];
% Create the expect script file
fid = fopen(expect_script, 'w');
fprintf(fid, '#!/usr/bin/expect -f\n');
fprintf(fid, 'set prompt "Enter your NIMH Data Archives username:"\n');
fprintf(fid, 'set timeout 300\n');
fprintf(fid, 'set pattern "Completed download:"\n');
fprintf(fid, 'spawn %s\n', subjDlCommand);
fprintf(fid, 'expect $prompt\n');
fprintf(fid, 'send "apines\\r"\n');
fprintf(fid, 'expect {\n');
fprintf(fid, '  timeout { \n');
fprintf(fid, '    puts "Timed out waiting for pattern"\n');
fprintf(fid, '    exit 1 \n');
fprintf(fid, '  }\n');
fprintf(fid, '  $pattern { \n');
fprintf(fid, '    expect eof\n');
fprintf(fid, '  }\n');
fprintf(fid, '}\n');
fclose(fid);
% Make the expect script executable
system(sprintf('chmod +x %s', expect_script));
% Call the expect script from MATLAB
system(sprintf('%s', expect_script));

%%% Motion mask
apply_motion_mask(subj)

%%% Concatenate scans
concat_TS(subj)

%%% SSP Workflow
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/gp/FunctionalImages/Networks')
% input TS
TSfp=['/scratch/users/apines/abcd_images/fmriresults01/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/' subj '_ses-baselineYear1Arm1_p2mm_masked_concat.dtseries.nii'];
%% SSP
PersonalizeNetworks(TSfp,18,subj)

% this can go in FC .m script
% re-mask with SNR to be sure: looks like 1's could be getting de-facto'ed into SNR empties after the actual NMF fit
%

% FC (personalized networks, subcortical Tian S1, Circuits) 
Extract_FC(subj)

% Sulc extract

% MM extract

%%% OpFl Workflow

% Additional Motion mask

% OpFl

% Downsample networks

% Props relative to networks


% combine into circuit scores from resting state

%%% combine features into vector, save to permanent storage

%%% delete input data

% comment out deletion for a few subjs to QC
























% tell it where AWS tools are for downloads
%system('export PATH=/cbica/projects/abcdfnets/aws/dist/:$PATH')
%subjTxtCommand=['echo ' subj ' >> /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/' subj '.txt'];

% download that one subject's data
% ∆∆∆∆∆∆
% keep looking at github for updates on authent. fix
% ∆∆∆∆∆∆
%%% it works suckkkaaaas


% subjDlCommand=['python3 /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/download.py -i /scratch/users/apines/datastructure_manifest.txt -o /scratch/users/apines/ -s /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/' subj '.txt -l /oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/dl_logs -d /oak/stanford/groups/leanew1/users/apines/data_subsets_3_9_21_from9620.txt &']

% note: downloader tool does not seem to communicate when it is done to matlab
% added '&' and 'pause' so that matlab waits 5 minutes to proceed rather than getting caught up indefinitely
system(subjDlCommand)

% set parent file path
parentFP=['/scratch/users/apines/derivatives/abcd-hcp-pipeline/ses-baselineYear1Arm1/func/' subj ];
% and child file path
childFP=['/scratch/users/apines/derivatives/pinespipe/' subj '/'];
mkdir(childFP)

% now the matlab portions. Apply the motion mask to the downloaded data
%%% probably will lose psychopathology instances if we do this

%%%%%% extract FD and TRs passing threshold

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

inputTS_fp=['/scratch/users/apines/derivatives/abcd-hcp-pipeline/ses-baselineYear1Arm1/func/' subj ];

% derive an indivudalized parcellation (just parent filepath for inputTS_fp, g_ls will find 10k mgh extension within folder
for k=2:30
	PersonalizeNetworks(inputTS_fp,k,subj)
end

% for "cifti_read", interferes if added earlier
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/cifti-matlab/'));

% calculate multiscale FC
CalcFC(['/scratch/users/apines/' subj '/subcort.ptseries.nii'],CortTS_L,CortTS_R,subj,outFP)

% run optical flow pipeline ∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆
%%% calculate optical flow
OpFl_abcd(tsIn_L,tsIn_R,tsOut)
% set output filepaths to calculcate angular distances
rsOut=[childfp '/' subj '_' '_PGGDist_fs5.mat'];

% needed?
% calculate angular distances (change output filepath to something meaningful)
AngDCalcCmd=['/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/fs_5/run_Extract_BUTD_ResultantVecs_Gran_fs5.sh /share/software/user/restricted/matlab/R2018a/ ' tsOut ' ' childFP 'OpFl_timeseries_L_fs5.mat' ' ' childFP 'OpFl_timeseries_R_fs5.mat'];
system(AngDCalcCmd)

%%%%% VECTORIZE DATA


% delete input data
Delete_input_data(subj)
