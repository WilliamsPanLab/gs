function PersonalizeNetworks(inputTS_fp,K,ID_Str)

% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));

sname=ID_Str;
% Based on the group atlas, creating each subject's individual specific atlas
ProjectFolder=['/scratch/users/apines/abcd_images/derivatives/abcd-hcp-pipeline/' sname '/ses-baselineYear1Arm1/func/'];
ResultantFolder = [ProjectFolder '/SingleParcel_1by1'];
mkdir(ResultantFolder);
% load in "prep data file"
PrepDataFile = ['/oak/stanford/groups/leanew1/users/apines/data/gp/CreatePrepData.mat'];
resId = 'IndividualParcel_Final';
% load in group consensus from pines et al 2022
initName=['/oak/stanford/groups/leanew1/users/apines/maps/initmat_18fslr.mat'];
% Use parameter in Hongming's NeuroImage paper, except for alphaL which we have expanded to account for the reduced smoothin in DCAN's pipeline
alphaS21 = 1;
alphaL = 2500;
vxI = 1;
spaR = 1;
ard = 0;
iterNum = 20;
eta = 0;
calcGrp = 0;
parforOn = 0;
% load in filepath for subject dtseries
CiftiCell = g_ls(inputTS_fp);

ResultantFolder_I = [ResultantFolder '/k' num2str(K) '/' ID_Str];
ResultFile_check = dir([ResultantFolder_I, '/**/final_UV.mat']);

% check for existing directory
if ~exist(ResultantFolder_I, 'dir')
   mkdir(ResultantFolder_I);
end
IDMatFile = [ResultantFolder_I '/ID.mat'];
save(IDMatFile, 'ID_Str');

sbjListFile = [ResultantFolder_I '/sbjListAllFile.txt'];
system(['rm -rf ' sbjListFile]);

% as two sep. lines so nmf can read both
cmd = ['echo ' CiftiCell{1} ' >> ' sbjListFile];
system(cmd);


deployFuncMvnmfL21p1_func_cifti(sbjListFile,PrepDataFile,ResultantFolder_I,resId,initName,K,alphaS21,alphaL,vxI,spaR,ard,eta,iterNum,calcGrp,parforOn)


