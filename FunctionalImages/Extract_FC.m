function Extract_FC(subj)
% add path for cifti/gifti functions
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));

% just spent about 10 hours trying to nancy-drew wb_documentation, and it looks like it's actually easier/more efficient to just take the mni->gifti metric files and manually index them in matlab rather than navigate/verify wb_labrynith 

% set rsfrmi fp
rsfp=['/scratch/users/apines/abcd_images/fmriresults01/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/' subj '_ses-baselineYear1Arm1_task-rest_p2mm_masked.dtseries.nii'];
% set concatenated fmri fp
cfp=['/scratch/users/apines/abcd_images/fmriresults01/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/' subj '_ses-baselineYear1Arm1_p2mm_masked_concat.dtseries.nii'];

% circuit ROI fps
% DMN
D1=['~/2021-masks/Medial_amPFC_DefaultModeNetwork_n2_50_n6.nii.gz.shape.gii'];
D2=['~/2021-masks/Left_AG_DefaultModeNetwork_n46_n70_32.nii.gz.shape.gii'];
D3=['~/2021-masks/Right_AG_DefaultModeNetwork_50_n62_26.nii.gz.shape.gii'];
D4=['~/2021-masks/Medial_PCC_DefaultModeNetwork_0_n50_28.nii.gz.shape.gii']'
% Salience
S1=['~/2021-masks/Left_antInsula_Salience_n38_14_n6.nii.gz.shape.gii'];
S2=['~/2021-masks/Right_antInsula_Salience_38_18_2.nii.gz.shape.gii'];
% use coarse Tian to extract amyg
% Attention
A1=['~/2021-masks/Medial_msPFC_Attention_n2_14_52.nii.gz.shape.gii'];
A2=['~/2021-masks/Left_lPFC_Attention_n44_6_32.nii.gz.shape.gii'];
A3=['~/2021-masks/Right_lPFC_Attention_50_10_28.nii.gz.shape.gii'];
A4=['~/2021-masks/Left_aIPL_Attention_n30_n54_40.nii.gz.shape.gii'];
A5=['~/2021-masks/Right_aIPL_Attention_38_n56_48.nii.gz.shape.gii'];
A6=['~/2021-masks/Left_precuneus_Attention_n14_n66_52.nii.gz.shape.gii'];
A7=['~/2021-masks/Right_precuneus_Attention_18_n68_52.nii.gz.shape.gii'];

% load in concatenated timne series
cts=read_cifti(cfp);
rts=read_cifti(rsfp);
% load in single subject parcellation
sspFP=['/scratch/users/apines/abcd_images/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/SingleParcel_1by1/k18/' subj '/IndividualParcel_Final_sbj1_comp18_alphaS21_1_alphaL2500_vxInfo1_ard0_eta0/final_UV.mat'];
ssp=load(sspFP);
% convert to hard parcellation
fV=[ssp.V{:}];
[~, HardParcel] = max(fV, [], 2);
% 0-indices
indices=sum(fV,2)==0;
% load in S1 Tian atlas
S1=read_cifti('~/Tian_Subcortex_S1_3T_32k.dscalar.nii')
% add subcort indices to network indices so we have 34 uniquely labeled ROIs
HardParcel=HardParcel+max(S1.cdata);
print('Do not forget you added 16 to all network indices, Adam.') 
% 0-indices to 0-network assignment
HardParcel(indices)=0;
% combine this participants functional networks with subcortical parcellations
S1.cdata(1:59412)=HardParcel;
% temp test
%write_cifti(S1,'~/Tian_Subcortex_S1_3T_32k_wHP.dscalar.nii');

%%%% calculate FC

%% add parallel labels like you have previously 

% because each network is treated as an ROI, we are going to treat the subcortical labels the same way for a total of 18 cortical networks + 16 subcortical rois for a 34x34 fc matrix, and then remove the diagonal/one triangles
for N=1:34
	% get indices of this network
	NetInds=S1.cdata==N;
	% get time series
	NetworkTS=cts.cdata(NetInds,:);
	% loop over every other network
	for OtherNetNumber=1:N
		% get indices for this OtherNetwork
		ONetInds=S1.cdata==OtherNetNumber;
		% get time series matrix
		OtherNetworkTS=cts.cdata(ONetInds,:);
		% get correlation matrix
		GrayordCoords=corr(OtherNetworkTS',NetworkTS');
		% return average
	end
end


% calculate circuit scores
