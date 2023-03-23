% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));

% load in init consensus map
initConsensus=load('/oak/stanford/groups/leanew1/users/apines/data/init_Vs/18/init_UV.mat');
V=initConsensus.V;
% extract from struct
V=V{1};
%%% use snr mask to translate back to 10242 per hemisphere
% load in mask (SNR Mask)
surfML = '/oak/stanford/groups/leanew1/users/apines/data/gp/lh.Mask_SNR.label';
mwIndVec_l = read_medial_wall_label(surfML);
% get difference for valid vertices
index_l=setdiff([1:10242],mwIndVec_l);
% right hemisphere
surfMR = '/oak/stanford/groups/leanew1/users/apines/data/gp/rh.Mask_SNR.label';
mwIndVec_r = read_medial_wall_label(surfMR);
index_r=setdiff([1:10242],mwIndVec_r);

% use left and right hemi PG func.giis as templates to save parcellation over
LH=gifti('/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients_10k_L.dscalar.func.gii')
RH=gifti('/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients_10k_R.dscalar.func.gii')

% and get a dtseries template
dtsTemplateFP='/oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/nda-abcd-s3-downloader/test_out/fmriresults01/derivatives/abcd-hcp-pipeline/sub-NDARINV04CLBZAD/ses-baselineYear1Arm1/func/sub-NDARINV04CLBZAD_ses-baselineYear1Arm1_task-rest_bold_desc-filtered_timeseries.dtseries.nii'

% set number of networks
k=18
% now for each network, plop mat data into funcgii and save a copy for merging back together at the end of this script (to then be plopped into a 32k-equiv initMat)	
for n=1:k
	% put current network into LeftHemi and RightHemi
	% zero-out current vector
	LH.cdata(:,n)=0;
	RH.cdata(:,n)=0;
	% replace with indices from V
	LH.cdata(index_l,n)=V(1:length(index_l),n);
	RH.cdata(index_r,n)=V(length(index_l)+1:end,n);
	% set filepath
	LHfp=['/oak/stanford/groups/leanew1/users/apines/maps/' num2str(n) '_groL.func.gii'];
	RHfp=['/oak/stanford/groups/leanew1/users/apines/maps/' num2str(n) '_groR.func.gii'];
	% saveout
	save(LH, LHfp);
	save(RH, RHfp);	
	% sys command for resampling
	commandL=['wb_command -metric-resample ' LHfp ' /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii ADAP_BARY_AREA /oak/stanford/groups/leanew1/users/apines/resamp_32k_network' num2str(n) '_L.32k_fs_LR.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii'];
	commandR=['wb_command -metric-resample ' RHfp ' /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii ADAP_BARY_AREA /oak/stanford/groups/leanew1/users/apines/resamp_32k_network' num2str(n) '_R.32k_fs_LR.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii'];
	system(commandL)
	system(commandR)
	% downsampled and upsampled func.gii's should be done for this network
end


% next step is to convert upsampled into init.mat again, next loop
% initialize initmat
initmat=zeros(59412,k);
for n=1:k
	% gifti name
	gifNameL=['/oak/stanford/groups/leanew1/users/apines/resamp_32k_network' num2str(n) '_L.32k_fs_LR.func.gii'];
	gifNameR=['/oak/stanford/groups/leanew1/users/apines/resamp_32k_network' num2str(n) '_R.32k_fs_LR.func.gii'];
	% combine them with connectome wb
	combinedFn=['/oak/stanford/groups/leanew1/users/apines/resamp_32k_network' num2str(n) '_LR.32k_fs_LR.dscalar.nii'];
	cmd = ['wb_command -cifti-create-dense-scalar ' combinedFn ' -left-metric ' gifNameL ' -right-metric ' gifNameR];
	system(cmd);
	% now resolve to 59k cortical vertices version
	outFN= ['/oak/stanford/groups/leanew1/users/apines/resamp_32k_network' num2str(n) '_LR.32k_fs_LR.dscalar.nii']; 
	cmd2 = ['wb_command -cifti-create-dense-from-template ' dtsTemplateFP ' ' outFN ' -cifti ' combinedFn];
	system(cmd2);
	% read it back in
	combined=read_cifti(outFN);
	combinedHemis=combined.cdata;
	% populate initmat with each upsampled gifti
	initmat(:,n)=combinedHemis(1:59412,n);
end
% match initmat format
V={initmat};
% saveout initmat
save('/oak/stanford/groups/leanew1/users/apines/maps/initmat_18fslr.mat','V');

% and make a group consensus hard parcellation for reference
[~, AtlasLabel_NoMedialWall] = max(initmat, [], 2);

% cifti to replace cdata in
combined.cdata(1:59412)=AtlasLabel_NoMedialWall;
outputfile=['~/K18_Parcel.dscalar.nii'];
write_cifti(combined,outputfile)
