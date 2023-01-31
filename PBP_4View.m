function PBP_vertWiseEffect5View(VertVec,name)% pretty picture code, AAB 4/2018 - AP 5/1/20 - Updated to take in 17734 vectors on 9/10
% data should be vectors, 10242 in length if fsaverage5 is used
% if using higher resolution, then change accordingly
% depencies include: matlab freesrufer functions, subaxis.m (matlab central), inferno color scale (matlab central - for Sam ;)


addpath(genpath('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox'));
ProjectFolder = '/cbica/projects/pinesParcels/data/SingleParcellation';
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage5';


surfML = '/cbica/projects/pinesParcels/data/H_SNR_masks/lh.Mask_SNR.label';
mwIndVec_l = read_medial_wall_label(surfML);
Index_l = setdiff([1:10242], mwIndVec_l);
surfMR = '/cbica/projects/pinesParcels/data/H_SNR_masks/rh.Mask_SNR.label';
mwIndVec_r = read_medial_wall_label(surfMR);
Index_r = setdiff([1:10242], mwIndVec_r);


% load pc age effects
%pcfn=['/cbica/projects/pinesParcels/results/EffectVecs/fSlope.mat'];
%pcstruct=load(pcfn);
%PcAgeEff=pcstruct.fSlope;

%Krange=2:30;

%for K=Krange
%K
%plot_text=[num2str(K)];
sbj_AtlasLabel_NoMedialWall=VertVec;
%sbj_AtlasLabel_NoMedialWall=PcAgeEff(:,K-1);
sbj_AtlasLabel_lh = zeros(1, 10242);
sbj_AtlasLabel_lh(Index_l) = sbj_AtlasLabel_NoMedialWall(1:length(Index_l));
sbj_AtlasLabel_rh = zeros(1, 10242);
sbj_AtlasLabel_rh(Index_r) = sbj_AtlasLabel_NoMedialWall(length(Index_l) + 1:end);
% sep. out lh and rh, merge it back with SNR mask

%%% read in data, left hemisphere
%scale=num2str(scale);
%scale=num2str(K);
%plot_text=[num2str(K)];
plot_text='';
%plot_text=strcat([scale ' communities, + only, Nodal Segreg. Avg']);
%readleft=['/cbica/projects/pinesParcels/multiscale/results/viz/k' scale '_loading_L.csv'];
% added in for SD
%readleft=['/cbica/home/pinesa/multiscale/results/viz/posseg_ROI_vertices_agedcor_SD_L.csv'];
%readright=['/cbica/home/pinesa/multiscale/results/viz/k' scale '_loading_R.csv'];
% added in for SD
%readright=['/cbica/home/pinesa/multiscale/results/viz/posseg_ROI_vertices_agecor_SD_R.csv'];
[vertices, faces] = freesurfer_read_surf('/cbica/software/external/freesurfer/scientificlinux6/6.0.0/subjects/fsaverage5/surf/rh.inflated');
%using lh.gray will make more anatomical looking plot but harder to see into sulci
%datal=importdata(readleft);
datal=sbj_AtlasLabel_lh;
datal=datal(1,1:10242)';
%datar=importdata(readright);
datar=sbj_AtlasLabel_rh;
datar=datar(1,1:10242)';

%set NaN to 0
%I generally have the midcut region set to NaN
%in the csv files that I read in
indexNaNrh = find(isnan(datar));
indexNaNlh = find(isnan(datal));
datar(indexNaNrh)=0;
datal(indexNaNlh)=0;

%%% set color scale
datalr=[datal; datar];

% set each network to it's effect size

%mincol = -.05;
%maxcol = .05;
%AP% set to make white zero on all maps
maxabs=max(abs(datalr));
mincol=0-maxabs;
maxcol=maxabs;

%%%% FIXED COLORSCALE - DR2
maxcol=.05;
mincol=-.05;

%change above to set max/min manually or by other means
%custommap=colormap('spring'); %or whatever
% for white at 0
%custommap=colormap(b2r(-3,3));

% plasma, non-0-spanning color schame
custommap=colormap('plasma');
mincol=min(datalr);
maxcol=max(datalr);


% to set 0 to black
custommap(1,:)=[.75 .75 .75];


% right hemi
data=datar;
asub = subaxis(2,2,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);
%asub = subplot(4,2,1)
% note use of subaxis is to ged rid of white space around brains
% if you don't care about that, it's faster and less likely to cause
% issues if you use subplot instead
% if so, bet rid of all of the posnew stuff below

aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
view([90 0]);
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting gouraud; %phong;
material metal %shiny %metal;
shading flat;
camlight;
alpha(1)

asub = subaxis(2,2,4, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
view([90 0]);
rotate(aplot, [0 0 1], 180)
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting gouraud; %phong;
material metal %shiny %metal;
shading flat;
camlight;
alpha(1)
 pos = get(asub, 'Position');
  posnew = pos; posnew(2) = posnew(2) + 0.13; posnew(1) = posnew(1) -.11; set(asub, 'Position', posnew);
  set(gcf,'Color','w')

  %%% left hemisphere
  data=datal;

  [vertices, faces] = freesurfer_read_surf('/cbica/software/external/freesurfer/scientificlinux6/6.0.0/subjects/fsaverage5/surf/lh.inflated');

  asub = subaxis(2,2,2, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0,'Holdaxis');
  aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
  view([90 0]);
  rotate(aplot, [0 0 1], 180)
  colormap(custommap)
  caxis([mincol; maxcol]);
  %caxis([NAval; max_data])
  daspect([1 1 1]);
  axis tight;
  axis vis3d off;
  lighting phong; %gouraud
  material metal %shiny %metal;
  shading flat;
  camlight;
   pos = get(asub, 'Position');
    posnew = pos; posnew(1) = posnew(1) - 0.11; set(asub, 'Position', posnew);
    alpha(1)
    %colormap(mycol)



    asub = subaxis(2,2,3, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);
    aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
    view([90 0]);
    colormap(custommap)
    caxis([mincol; maxcol]);
    daspect([1 1 1]);
    axis tight;
    axis vis3d off;
    lighting gouraud; %phong;
    material metal %shiny %metal;
    shading flat;
    camlight;
    alpha(1)
     pos = get(asub, 'Position');
      posnew = pos; posnew(2) = posnew(2) + 0.13; set(asub, 'Position', posnew);
      set(gcf,'Color','w')

      %%%


      acbar = colorbar('southoutside');
      acbar.FontWeight='bold';
      acbar.FontSize=13;
      acbar.TickLength=.025;
      acbar.LineWidth=3;

      % for 3 tick marks
      %acbar.TicksMode='manual';
      midval = round((maxcol-mincol)/2);
      %acbar.Ticks=[round(mincol),midval,floor(maxcol)];
      %if mincol < 0
      %       midval=((maxcol+mincol)/2);
      %       acbart.Ticks=[mincol,midval,maxcol];
      %end
      set(acbar, 'position', [0.35 0.55 0.2 0.02])

      % going lower rez for now, but giant vector rendering was beaut
      print('-dpng','-r600',['/cbica/projects/pinesParcels/results/viz/' char(name)])
      %print('-dpdf','-r600',['/cbica/home/pinesa/multiscale/results/viz/' char(name)])
      end

