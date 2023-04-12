function ConcatOpFl(subj)
% concatenate optical flow derivatives

% these are the tasks
tasks={'rest','MID','SST','nback'};

% initialize a structure to use
us.vf_left={};
us.vf_right={};
for t=1:4
	task=tasks{t};
	OpFlFile=['/scratch/users/apines/abcd_images/fmriresults01/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/' subj '_' task '_OpFl_3k.mat'];
	if exist(OpFlFile,'file')
		% load file
		OpFl=load(OpFlFile)
		us.vf_left=[us.vf_left,OpFl.us.vf_left];
		us.vf_right=[us.vf_right,OpFl.us.vf_right];
	else
	end
end
% save out
OpFlFile=['/scratch/users/apines/abcd_images/fmriresults01/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/' subj '_concat_OpFl_3k.mat'];
save(OpFlFile,'us')
