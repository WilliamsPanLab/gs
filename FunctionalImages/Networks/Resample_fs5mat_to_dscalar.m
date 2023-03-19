function Resample_fs5mat_to_dscalar(subj)

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
surfML = '/oak/stanford/groups/leanew1/users/apines/data/gp/lh.Mask_SNR.label'
mwIndVec_l = read_medial_wall_label(surfML);
% get difference for valid vertices
index_l=setdiff([1:10242],mwIndVec_l);
% right hemisphere
surfMR = '/oak/stanford/groups/leanew1/users/apines/maps/rh.cortex_SNR.label';
mwIndVec_r = read_medial_wall_label(surfMR);
index_r=setdiff([1:10242],mwIndVec_r);

% use left and right hemi PG func.giis as templates to save parcellation over
LH=gifti('/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients_10k_L.dscalar.func.gii')
RH=gifti('/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients_10k_R.dscalar.func.gii')

% set number of networks
k=18
% now for each network, plop mat data into funcgii and save a copy for merging back together at the end of this script (to then be plopped into a 32k-equiv initMat)	
for n=1:k
	% put current network into LeftHemi and RightHemi
	% insert that non-medial wall data from the consensus partition
	
	%%%LeftHemi(mwIndVec_l)=V(1:length(mwIndVec_l),n);
	%%%RightHemi(mwIndVec_r)=V(length(mwIndVec_l)+1:(length(mwIndVec_r)+length(mwIndVec_l)),n);

	% use left and right hemi PG func.giis as templates to save parcellation over
	LH=gifti('/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients_10k_L.dscalar.func.gii')
	RH=gifti('/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients_10k_R.dscalar.func.gii')
	% zero-out current vector
	LH.cdata(:,n)=0
	% replace with indices from V
	LH.cdata(index_l,n)=V(1:length(index_l),n);

		% ADOPT FASHION ABOVE ^

	RH.cdata(mwIndVec_r,n)=V(length(mwIndVec_l)+1:end),n);

	% set filepath
	LHfp=['/oak/stanford/groups/leanew1/users/apines/maps/' num2str(n) '_groL.func.gii'];
	RHfp=['/oak/stanford/groups/leanew1/users/apines/maps/' num2str(n) '_groR.func.gii'];
	% saveout
	save(LH, LHfp);
	save(RH, RHfp);

	
	% sys command for resampling





outputfile=['/cbica/projects/abcdfnets/results/SingleParcel_1by1/' subj '/' subj '_Parcel.dscalar.nii'];
cifti_write(HP,outputfile)

% loop over every k and save out sep. func.gii map
for i = 1:17;
	% load file
	fp=[unmasked_GC_folder 'Group_AtlasLoading_Network_' num2str(i) '.dscalar.nii'];
	NetLoads=cifti_read(fp);
	% extract LH
	lhLoad=NetLoads.cdata(1:10242);
	% initialize gifti
	LH_gif=gifti;
	LH_gif.cdata=lhLoad;
	V_lh_File = [unmasked_GC_folder 'Group_lh_Network_' num2str(i) '.func.gii'];
	% save
	save(LH_gif, V_lh_File);
	% extract RH
	rhLoad=NetLoads.cdata(10243:20484);
	% initialize gifti
	RH_gif=gifti;
	RH_gif.cdata=rhLoad;
	V_rh_File = [unmasked_GC_folder 'Group_rh_Network_' num2str(i) '.func.gii'];
	% save
	save(RH_gif, V_rh_File);
	i
end

