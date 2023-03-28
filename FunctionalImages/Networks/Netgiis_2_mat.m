function Netgiis_2_mat(subj)
% convert downsampled network .giis to .mats
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% this participant's filepath
funcgiiFolder = ['/scratch/users/apines/abcd_images/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/'];
netgiis_L=[funcgiiFolder '' subj '_L_AggNets_3k.func.gii'];
netgiis_R=[funcgiiFolder '' subj '_R_AggNets_3k.func.gii'];
% load in functional networks: Left
Lnets=gifti(netgiis_L);
% load in functional networks: Right
Rnets=gifti(netgiis_R);
% set to 18 networks
Lnets=Lnets.cdata(:,1:18);
Rnets=Rnets.cdata(:,1:18);
% save out
nets=struct;
nets.Lnets=Lnets;
nets.Rnets=Rnets;
save([funcgiiFolder '' subj '_Nets_fs4.mat'],'nets')
