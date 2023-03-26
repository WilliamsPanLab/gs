function Extract_BUTD_ResultantVecs_PL(subj,infileOpFl,infileAngD,outfile_L,outfile_R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take optical flow results, get a bottom-up and top-down resultant vector in x,y coords for each face. Measured relative to gPGG.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set assumed network input directory
NetworksFolder = ['/scratch/users/apines/abcd_images/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/SingleParcel_1by1/k18/' subj '/'];

% Load in fsav4 opflow calc
data=load(infileOpFl)
% Load in surface data
surfL = ['/oak/stanford/groups/leanew1/users/apines/surf/lh.sphere'];
surfR = ['/oak/stanford/groups/leanew1/users/apines/surf/rh.sphere'];
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
% vector fields
vfl=data.us.vf_left;
% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;
% vector fields
vfr=data.us.vf_right;
% get incenters of triangles
TR_L = TriRep(F_L,V_L);
P_L = TR_L.incenters;
TR_R = TriRep(F_R,V_R);
P_R = TR_R.incenters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use native freesurfer command for mw mask indices
surfML = '/oak/stanford/groups/leanew1/users/apines/surf/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/oak/stanford/groups/leanew1/users/apines/surf/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "is medial wall" vector for vertices
mw_L=zeros(1,2562);
mw_L(mwIndVec_l)=1;
mw_R=zeros(1,2562);
mw_R(mwIndVec_r)=1;
% convert to faces
% convert to faces
F_MW_L=sum(mw_L(faces_l),2)./3;
F_MW_R=sum(mw_R(faces_r),2)./3;
% convert "partial" medial wall to medial wall
F_MW_L=ceil(F_MW_L);
F_MW_R=ceil(F_MW_R);
% face mask indices
fmwIndVec_l=find(F_MW_L);
fmwIndVec_r=find(F_MW_R);
% make medial wall vector
g_noMW_combined_L=setdiff([1:5120],fmwIndVec_l);
g_noMW_combined_R=setdiff([1:5120],fmwIndVec_r);

% extract size of time series
vfl=data.us.vf_left;
vfr=data.us.vf_right;
NumTRs=size(vfl);
NumTRs=NumTRs(2);
lenOpFl=NumTRs;

% load in Networks
networks=load([NetworksFolder '_Nets_fs4.mat']);
nets_LH=networks.Lnets;
nets_RH=networks.Rnets;
% initialize out dataframes
Propvec=[];
stringVec={};
Thetas_L=zeros(length(g_noMW_combined_L),lenOpFl);
Thetas_R=zeros(length(g_noMW_combined_R),lenOpFl);
% load in time series, convert angles, then loop over each network for angular distances (network gradients iteratively)
% for each face
for F=g_noMW_combined_L
	% for each timepoint
	for fr=1:lenOpFl
		% current vector field
		relVf_L=vfl{fr};
		% xyz components
        	xComp_L=relVf_L(F,1);
        	yComp_L=relVf_L(F,2);
        	zComp_L=relVf_L(F,3);
		% convert to spherical coord system
        	vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(F),eld_L(F));
		% convert to spherical coordinates
       		OpFlVec_L= [vs_L(1) vs_L(2)];
        	% store in output vector (r is redundant across all vecs, only using az and el)
		[Thetas_L(F,fr),~]=cart2pol(vs_L(1),vs_L(2));
	end
end

% for each face
for F=g_noMW_combined_R
	% for each timepoint
	for fr=1:lenOpFl
		% current vector field
		relVf_R=vfr{fr};
		% xyz components
		xComp_R=relVf_R(F,1);
		yComp_R=relVf_R(F,2);
		zComp_R=relVf_R(F,3);
		% convert to spherical coord system
		vs_R=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(F),eld_R(F));
		% convert to spherical coordinates
       		OpFlVec_R= [vs_R(1) vs_R(2)];
		% store in output vector (r is redundant across all vecs, only using az and el)
		[Thetas_R(F,fr),~]=cart2pol(vs_R(1),vs_R(2));
	end
end
% use thetas to calculate angular distance
% for each network
for k=1:18
	% network of interest
	n_LH=nets_LH(:,k);
	n_RH=nets_RH(:,k);
	% calculate network gradients on sphere
	ng_L = grad(F_L, V_L, n_LH);
	ng_R = grad(F_R, V_R, gPG_RH);
	% get NA vertices and mask them out medial wall mask

	% extract face-wise vector cartesian vector components
	nx_L=ng_L(:,1);
	ny_L=ng_L(:,2);
	nz_L=ng_L(:,3);
	nx_R=ng_R(:,1);
	ny_R=ng_R(:,2);
	nz_R=ng_R(:,3);

	% translate xyz spherical coordinates to az/el/r
	[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
	[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));

	% convert from radians to degrees
	azd_L=rad2deg(az_L);
	eld_L=rad2deg(el_L);
	azd_R=rad2deg(az_R);
	eld_R=rad2deg(el_R);
	
	% translate xyz vector components at coordinates to az/el/r
	gazes_L=zeros(1,length(azd_L));
	gels_L=zeros(1,length(eld_L));
	for i=1:length(azd_L)
	    gvs_L=cart2sphvec(double([nx_L(i);ny_L(i);nz_L(i)]),azd_L(i),eld_L(i));
	    gazes_L(i)=gvs_L(1);
	    gels_L(i)=gvs_L(2);
	end
	% right hemi
	gazes_R=zeros(1,length(azd_R));
	gels_R=zeros(1,length(eld_R));
	for i=1:length(azd_R)
	    gvs_R=cart2sphvec(double([nx_R(i);ny_R(i);nz_R(i)]),azd_R(i),eld_R(i));
	    gazes_R(i)=gvs_R(1);
	    gels_R(i)=gvs_R(2);
	end

	% calculate angular distances
	Nangs_L=cart2pol(gazes_L,gels_L);
	% include masking
	NangDs_L=abs(Nangs_L-Thetas_L);
	% right
	Nangs_R=cart2pol(gazes_R,gels_R);
	% include masking
	NangDs_R=abs(Nangs_R-Thetas_R);
	% average angular distances across hemispheres
	avgD=mean([NangDs_L NangDs_R],2);
	Propvec=[Propvec avgD];
	% add label - adding 16 to be consistent with 
	stringVec=[stringVec ['AngD' num2str(N+16)]];
end
% save out as csv
T=table(Propvec','RowNames',stringVec);
% calc outFP
outFP=['/oak/stanford/groups/leanew1/users/apines/data/gp/PropFeats/' subj];
% make out filepath
system(['mkdir ' outFP]);
% write out
writetable(T,[outFP '/Prop_Feats.csv'],'WriteRowNames',true)
