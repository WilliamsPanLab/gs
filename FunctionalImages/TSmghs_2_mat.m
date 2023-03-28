function TSmghs_2_mat(ts)
% add freesurfer path
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));
% load in time series
data=MRIread(ts).vol;
% rid it of extra dimensions
data=squeeze(data);
% save out as same name except .mat at end
outfp=[ts '.mat'];
save(outfp,'data')
