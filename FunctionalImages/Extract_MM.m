function Extract_MM(subj)
% extract MyelinMaps from personalized functional network ROIs 

% add path for cifti/gifti functions
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));

% load in single subject parcellation
sspFP=['/scratch/users/apines/abcd_images/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/SingleParcel_1by1/k18/' subj '/IndividualParcel_Final_sbj1_comp18_alphaS21_1_alphaL2500_vxInfo1_ard0_eta0/final_UV.mat'];
ssp=load(sspFP);
% convert to hard parcellation
fV=[ssp.V{:}];
[~, HardParcel] = max(fV, [], 2);
% 0-indices
indices=sum(fV,2)==0;
HardParcel(indices)=0;

% load in MM
MMFP=['/scratch/users/apines/abcd_images/imagingcollection01/derivatives/abcd-hcp-pipeline/sub-NDARINV1JGK90HZ/ses-baselineYear1Arm1/anat/' subj '_ses-baselineYear1Arm1_space-fsLR32k_desc-smoothed_myelinmap.dscalar.nii'];

% initialize output FC table/struct/vector
MMvec=[];
stringVec={};

% extract the 18 CT values
for N=1:18
    % get the indices for this parcel
    parcelInds=find(HardParcel==N);
    % get the CT values for this parcel
    MMvals=read_cifti(MMFP)
    MMNetVals=MMvals.cdata(parcelInds);
    % average across the parcel
    MMvals=mean(MMNetVals);
    % add to the output vector
    MMvec=[MMvec MMvals];
    % add to the output string vector
	% note add of +16 to keep numbers consistent with subcortical->cortical ROI sequence
    stringVec=[stringVec ['MM' num2str(N+16)]];
end

% save out as csv
T=table(MMvec','RowNames',stringVec);
% calc outFP
outFP=['/oak/stanford/groups/leanew1/users/apines/data/gp/anat_Feats/' subj];
% make out filepath
system(['mkdir ' outFP]);
% write out
writetable(T,[outFP '/MM_Feats.csv'],'WriteRowNames',true)
