function mat_to_dlabel(subj)

% add matlab paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));

% directories
ProjectFolder = '/scratch/users/apines/abcd_images/derivatives/abcd-hcp-pipeline/';
OutputFolder = ['/scratch/users/apines/abcd_images/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/SingleParcel_1by1/k18/' subj '/'];

% resampled V of interest (group or individ - CHANGE AS FIT TO MATCH NAME)
finalUVFile=[OutputFolder 'IndividualParcel_Final_sbj1_comp18_alphaS21_1_alphaL2500_vxInfo1_ard0_eta0/final_UV.mat'];
Loading_Mat=load(finalUVFile);
V=Loading_Mat.V;

% extract from struct
V=V{:};

% cifti to replace cdata in
HP=read_cifti('/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients.dscalar.nii');
HP.cdata(1:59412,1:18)=V;
outputfile=[OutputFolder 'SoftParcel.dscalar.nii'];
write_cifti(HP,outputfile)
