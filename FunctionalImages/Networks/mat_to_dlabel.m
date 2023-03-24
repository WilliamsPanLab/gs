function mat_to_dlabel(subj)

% add matlab paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));

%%%
ProjectFolder = '/scratch/users/apines/abcd_images/derivatives/abcd-hcp-pipeline/';
OutputFolder = ['/scratch/users/apines/abcd_images/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/SingleParcel_1by1/k18/' subj '/'];
%%%

% resampled V of interest (group or individ - CHANGE AS FIT TO MATCH NAME)
finalUVFile=[OutputFolder 'IndividualParcel_Final_sbj1_comp18_alphaS21_1_alphaL2500_vxInfo1_ard0_eta0/final_UV.mat'];
Loading_Mat=load(finalUVFile);
V=Loading_Mat.V;
% extract from struct
V=V{:};

% apply ZC code from step 6th - subject AtlasLabel will have parcels for viz
V_Max = max(V);
trimInd = V ./ max(repmat(V_Max, size(V, 1), 1), eps) < 5e-2;
V(trimInd) = 0;
sbj_AtlasLoading_NoMedialWall = V;
[~, sbj_AtlasLabel_NoMedialWall] = max(sbj_AtlasLoading_NoMedialWall, [], 2);

% cifti to replace cdata in
HP=read_cifti('/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients.dscalar.nii');
HP.cdata(1:59412)=sbj_AtlasLabel_NoMedialWall;
outputfile=[OutputFolder 'HardParcel.dscalar.nii'];
write_cifti(HP,outputfile)
