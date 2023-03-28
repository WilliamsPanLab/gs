function Extract_SCV(subj)

% assume aseg path
asFold=['/scratch/users/apines/abcd_images/imagingcollection01/derivatives/freesurfer-5.3.0-HCP/' subj '/ses-baselineYear1Arm1/stats/']
asfp=[asFold 'aseg.stats'];

% parse aseg table
asegCmd=['cat ' asfp ' | grep -v "#" > ' asFold 'aseg_stats.csv'];
system(asegCmd);

% load in parsed as table
stats=readtable([asFold 'aseg_stats.csv']);

% normal format
format long g

% initialize output SCV table/struct/vector
SCVvec=stats.Var4;
stringVec=stats.Var5;

% save as output to feature directory
T=table(SCVvec,'RowNames',stringVec);
% calc outFP
outFP=['/oak/stanford/groups/leanew1/users/apines/data/gp/anat_Feats/' subj];
% make out filepath
system(['mkdir ' outFP]);
% write out
writetable(T,[outFP '/SCV_Feats.csv'],'WriteRowNames',true)
