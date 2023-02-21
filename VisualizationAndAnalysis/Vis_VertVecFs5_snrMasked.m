function Vis_Vertvec(VertVec,Fn) 

addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))

%%% Load in surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/fs5surf';
surfL = [SubjectsFolder '/lh.inflated'];
surfR = [SubjectsFolder '/rh.inflated'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% faces_L
F_L=faces_l;
% vertices V
V_L=vx_l;
% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;
%%%% medial wall stuff

% use native freesurfer command for mw mask indices
surfML = '/oak/stanford/groups/leanew1/users/apines/data/gp/lh.Mask_SNR.label';
surfMR = '/oak/stanford/groups/leanew1/users/apines/data/gp/rh.Mask_SNR.label';
mwIndVec_l = read_medial_wall_label(surfML);
mwIndVec_r = read_medial_wall_label(surfMR);
Index_l = setdiff([1:10242], mwIndVec_l);
Index_r = setdiff([1:10242], mwIndVec_r);
%%%% make binary "isn't medial wall" vector for vertices
%%%%mw_L=ones(1,10242);
%%%%mw_L(mwIndVec_l)=0;
%%%%mw_R=ones(1,10242);
%%%%mw_R(mwIndVec_r)=0;

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%plotdata=zeros(1,10242);
sbj_AtlasLabel_NoMedialWall=VertVec;
sbj_AtlasLabel_lh = zeros(1, 10242);
sbj_AtlasLabel_lh(Index_l) = sbj_AtlasLabel_NoMedialWall(1:length(Index_l));
sbj_AtlasLabel_rh = zeros(1, 10242);
sbj_AtlasLabel_rh(Index_r) = sbj_AtlasLabel_NoMedialWall(length(Index_l) + 1:end);

plotdata=sbj_AtlasLabel_lh;

%%%%%%% fixed colorscale varities

%%% circular
%mincol=-10;
%maxcol=10;


%%% for red/blue 0-centered
mincol=0;
maxcol=max(VertVec);
%custommap=colormap(b2r(mincol,maxcol));
% abscense of color to gray to accom. lighting "none"
%custommap(126,:)=[.5 .5 .5];
% blue-orange color scheme
%BO_cm=inferno(9);
%BO_cm(1,:)=[49 197 244];
%BO_cm(2,:)=[71 141 203];
%BO_cm(3,:)=[61 90 168];
%BO_cm(4,:)=[64 104 178];
%BO_cm(5,:)=[126 126 126];
%BO_cm(6,:)=[240 74 35];
%BO_cm(7,:)=[243 108 33];
%BO_cm(8,:)=[252 177 11];
%BO_cm(9,:)=[247 236 31];
% scale to 1
%BO_cm=BO_cm.*(1/255);
% interpolate color gradient
%interpsteps=[0 .125 .25 .375 .5 .625 .75 .875 1];
%BO_cm=interp1(interpsteps,BO_cm,linspace(0,1,255));
%custommap=BO_cm;
custommap=colormap(parula);
custommap(1,:)=[.5 .5 .5];

%%% matches circular hist
% for 180 degree max
%roybigbl_cm=inferno(6);
%roybigbl_cm(1,:)=[0, 0, 255];
%roybigbl_cm(2,:)=[0, 255, 255];
%roybigbl_cm(3,:)=[116, 192, 68];
%roybigbl_cm(4,:)=[246, 235, 20];
%roybigbl_cm(5,:)=[255, 165, 0];
%roybigbl_cm(6,:)=[255, 0, 0];
% scale to 1
%roybigbl_cm=roybigbl_cm.*(1/255);
% interpolate color gradient
%interpsteps=[0 .2 .4 .6 .8 1];
%roybigbl_cm=interp1(interpsteps,roybigbl_cm,linspace(0,1,255));
% make circular with flipud
%custommap=vertcat(flipud(roybigbl_cm),roybigbl_cm);

% for 90 degree max
%roybigbl_cm=inferno(3);
% blue
%roybigbl_cm(1,:)=[0, 0, 255];
% cyan
%roybigbl_cm(2,:)=[0, 255, 255];
% green
%roybigbl_cm(3,:)=[116, 192, 68];
% yellow
%roybigbl_cm(4,:)=[246, 235, 20];
% scale to 1
%roybigbl_cm=roybigbl_cm.*(1/255);
% interpolate color gradient
%interpsteps=[0 .33333 .66666 1];
%interpsteps=[0 .5 1];
%roybigbl_cm=interp1(interpsteps,roybigbl_cm,linspace(0,1,255));
% make circular with flipud
%custommap=vertcat(flipud(roybigbl_cm),roybigbl_cm);
%custommap=roybigbl_cm;
%custommap=colormap(parula);
% mw to black
%custommap(1,:)=[0 0 0];

figure
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/lh.inflated']);
asub = subaxis(2,2,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);

aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),plotdata)
view([90 0]);
colormap(custommap)
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting gouraud;
shading flat;
camlight;
	alpha(1)

set(gca,'CLim',[mincol,maxcol]);
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');

asub = subaxis(2,2,4, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),plotdata)
view([90 0]);
rotate(aplot, [0 0 1], 180)
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting gouraud;
material metal %shiny %metal;
shading flat;
camlight;
alpha(1)
 pos = get(asub, 'Position');
 posnew = pos; posnew(2) = posnew(2) + 0.13; posnew(1) = posnew(1) -.11; set(asub, 'Position', posnew);
set(gcf,'Color','w')

set(gca,'CLim',[mincol,maxcol]);
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');


%%% right hemisphere
plotdata=sbj_AtlasLabel_rh;

[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/rh.inflated']);

asub = subaxis(2,2,2, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0,'Holdaxis',1);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),plotdata)
view([90 0]);
rotate(aplot, [0 0 1], 180)
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting gouraud;
material metal %shiny %metal;%shading flat;
shading flat;
camlight;
 pos = get(asub, 'Position');
 posnew = pos; posnew(1) = posnew(1) - 0.11; set(asub, 'Position', posnew);
alpha(1)


set(gca,'CLim',[mincol,maxcol]);
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');

asub = subaxis(2,2,3, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),plotdata)
view([90 0]);
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting gouraud;
material metal %shiny %metal;
shading flat;
camlight;
alpha(1)
 pos = get(asub, 'Position');
 posnew = pos; posnew(2) = posnew(2) + 0.13; set(asub, 'Position', posnew);
set(gcf,'Color','w')


set(gca,'CLim',[mincol,maxcol]);
%%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');
%colorbar
%c=colorbar
%c.Location='southoutside'

colormap(custommap)

print(Fn,'-dpng')
