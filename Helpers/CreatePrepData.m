addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));
SubjectsFolder = '/share/software/user/open/freesurfer/6.0.0/subjects/fsaverage5';

% for surface data
surfL = [SubjectsFolder '/surf/lh.pial'];
surfR = [SubjectsFolder '/surf/rh.pial'];

surfML = '/oak/stanford/groups/leanew1/users/apines/data/gp/lh.Mask_SNR.label';
surfMR = '/oak/stanford/groups/leanew1/users/apines/data/gp/rh.Mask_SNR.label';

[surfStru, surfMask] = getFsSurf(surfL, surfR, surfML, surfMR);

% uncomment to implement without SNR mask.
%surfMask.l=ones(length(surfMask.l),1);
%surfMask.r=ones(length(surfMask.r),1);

gNb = createPrepData('surface', surfStru, 1, surfMask);

% save gNb into file for later use
prepDataName = ['/oak/stanford/groups/leanew1/users/apines/data/gp/CreatePrepData.mat'];
save(prepDataName, 'gNb');
