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
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/gp/FunctionalImages/Networks');
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/gp/FunctionalImages/Props');

disp('Δ Downloading data')
tic

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

toc
disp('Δ applying motion masks and concatenating scans')
tic

%%% Motion mask + flag missing tasks
apply_motion_mask(subj)

%%% Concatenate scans
concat_TS(subj)

toc
disp('Δ mapping individual functional neuroanatomy')
tic

%%% SSP Workflow
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/gp/FunctionalImages/Networks')
% input TS
TSfp=['/scratch/users/apines/abcd_images/fmriresults01/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/' subj '_ses-baselineYear1Arm1_p2mm_masked_concat.dtseries.nii'];
%% SSP
PersonalizeNetworks(TSfp,18,subj)

toc
disp('Δ Extracting structural features')
tic

% FC (personalized networks, subcortical Tian S1, Circuits) 
Extract_FC(subj)
% CT extract
Extract_CT(subj)
% MM extract
Extract_MM(subj)

toc
disp('Δ Running OpFl workflow')
tic


%%% OpFl Workflow
% downsample time series
DScommand=['./DS_surf_ts.sh ' subj];
system(DScommand)

% convert personalized networks to dscalar to downsample
mat_to_dscalar(subj)

% Downsample networks with workbench
DSCommand=['Networks/DS_surf_Networks.sh ' subj];
system(DSCommand)

% convert them to .mat so it works in compiled matlab scripts
Netgiis_2_mat(subj)

% OpFl
tasks=["rest","MID","SST","nback"];
for t=tasks
	task=tasks(t);
	% set filepaths
	LeftTS=['/scratch/users/apines/abcd_images/fmriresults01/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/' subj '/' subj '_' task '_L_AggTS_3k.func.gii'];
	RightTS=['/scratch/users/apines/abcd_images/fmriresults01/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/' subj '/' subj '_' task '_R_AggTS_3k.func.gii'];
	OpFlOut=['/scratch/users/apines/abcd_images/fmriresults01/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/' subj '/' subj '_' task '_OpFl_3k.mat'];
	% if files exist, run optical flow
	if exist(LeftTS,'file') && exist(RightTS,'file')
	% run OpFl
	OpFl_abcd(subj,t,LeftTS,RightTS,['/scratch/users/apines/abcd_images/fmriresults01/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/' subj '/' subj '_' task '_OpFl_3k.mat'])
end
% Props relative to networks
Extract_BUTD_ResultantVecs_PL(subj,OpFlOut,infileAngD,outfile_L,outfile_R)
% note this might be a compiled script
toc
disp('Δ deleting neuroimages')
tic
%%% delete input data
% comment out deletion to QC a few runs
Delete_input_data(subj)
disp('Δ done Δ')
