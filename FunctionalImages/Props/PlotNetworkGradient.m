addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6'))
% load in fs4 surface
surfL = ['/oak/stanford/groups/leanew1/users/apines/surf/lh.sphere'];
[vx_l, faces_l] = read_surf(surfL);
faces_l = faces_l + 1;
F_L=faces_l;
V_L=vx_l;
TR_L = TriRep(F_L,V_L);
P_L = TR_L.incenters;
% load in fs4 network
testL=load('/scratch/users/apines/abcd_images/derivatives/abcd-hcp-pipeline/sub-NDARINV0CP9XGTP/ses-baselineYear1Arm1/func/sub-NDARINV0CP9XGTP_Nets_fs4.mat');
% extract 15th network
PG_LH=testL.nets.Lnets(:,15);
% get gradient of it
PGg_L = grad(faces_l, vx_l, PG_LH);
% get nan mask
g_noMW_combined_L=find(PG_LH(:,1));
% normalize vectors
u=PGg_L;
% initialize figure
figure('units','pixels','position',[0 0 2500 2500])
% plot vector fields
quiver3D([P(g_noMW_combined_L,1)./scalingfactor,P(g_noMW_combined_L,2)./scalingfactor,P(g_noMW_combined_L,3)./scalingfactor],[ret(g_noMW_combined_L,1), ret(g_noMW_combined_L,2), ret(g_noMW_combined_L,3)],[0.18 0.19 0.57])
hold on
% make surface
trisurf(faces_l, vx_l(:, 1)./scalingfactor, vx_l(:, 2)./scalingfactor, vx_l(:, 3)./scalingfactor, PG_LH, 'EdgeColor','none'); 
caxis([0,.8]) 
colormap('jet')
daspect([1 1 1])
view(60,60)
print('~/test_60_60.png','-dpng') 
