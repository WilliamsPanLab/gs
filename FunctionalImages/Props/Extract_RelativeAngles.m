function Extract_RelativeAngles(subj,infileOpFl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take optical flow results, get a bottom-up and top-down resultant vector in x,y coords for each face. Measured relative to gPGG.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set assumed network input directory
NetworksFolder = ['/scratch/users/apines/abcd_images/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/func/'];

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
networks=load([NetworksFolder subj '_Nets_fs4.mat']);
nets_LH=networks.nets.Lnets;
nets_RH=networks.nets.Rnets;
% initialize out dataframes
Propvec=[];
stringVec={};
% and azez and els for opflow vectors
OpF_azes_L=zeros(length(g_noMW_combined_L),lenOpFl);
OpF_els_L=zeros(length(g_noMW_combined_L),lenOpFl);
OpF_azes_R=zeros(length(g_noMW_combined_R),lenOpFl);
OpF_els_R=zeros(length(g_noMW_combined_R),lenOpFl);

% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));

% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);

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
       		%OpFlVec_L= [vs_L(1) vs_L(2)];
		OpF_azes_L(F,fr)=vs_L(1);
		OpF_els_L(F,fr)=vs_L(2);
		% store in output vector (r is redundant across all vecs, only using az and el)
		%[Thetas_L(F,fr),~]=cart2pol(vs_L(1),vs_L(2));
	end
end
% mask mw out of angles
OpF_azes_L=OpF_azes_L(g_noMW_combined_L,:);
OpF_els_L=OpF_els_L(g_noMW_combined_L,:);

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
       		%OpFlVec_R= [vs_R(1) vs_R(2)];
		OpF_azes_R(F,fr)=vs_R(1);
		OpF_els_R(F,fr)=vs_R(2);
		% store in output vector (r is redundant across all vecs, only using az and el)
		%[Thetas_R(F,fr),~]=cart2pol(vs_R(1),vs_R(2));
	end
end
% mask mw out of angles
OpF_azes_R=OpF_azes_R(g_noMW_combined_R,:);
OpF_els_R=OpF_els_R(g_noMW_combined_R,:);

% add lukas lang functions
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6'))

% convert azs and els to only those within MW mask
az_L=az_L(g_noMW_combined_L);
el_L=el_L(g_noMW_combined_L);
az_R=az_R(g_noMW_combined_R);
el_R=el_R(g_noMW_combined_R);

% get network gradient to calculate relative angles
% for each network
for k=1:18
	% network of interest
	n_LH=nets_LH(:,k);
	n_RH=nets_RH(:,k);
	% calculate network gradients on sphere
	ng_L = grad(F_L, V_L, n_LH);
	ng_R = grad(F_R, V_R, n_RH);
	% use medial wall mask as common starting point (from which to mask both opfl vecs and net grads further)
	ng_L=ng_L(g_noMW_combined_L,:);
	ng_R=ng_R(g_noMW_combined_R,:);
	% get NA vertices
        sumLeft=sum(ng_L,2);
	sumRight=sum(ng_R,2);
	% finds 0s in left and right network gradients
	emptyLeft=find(~sumLeft);
	emptyRight=find(~sumRight);
	InclLeft=find(sumLeft);
	InclRight=find(sumRight);
	% not InclLeft and Right presume mw mask already applied!

	% mask them out of medial wall mask (medial wall mask indicates what to include, emptyLeft indicates what to exclude. setdiff excludes what should be excluded (from eL) from what should be incl. (noMW)
	%n_and_g_noMW_combined_L=setdiff(g_noMW_combined_L,emptyLeft);
	%n_and_g_noMW_combined_R=setdiff(g_noMW_combined_R,emptyRight);
	% extract face-wise vector cartesian vector components
	
	nx_L=ng_L(InclLeft,1);
	ny_L=ng_L(InclLeft,2);
	nz_L=ng_L(InclLeft,3);
	nx_R=ng_R(InclRight,1);
	ny_R=ng_R(InclRight,2);
	nz_R=ng_R(InclRight,3);

	% translate xyz spherical coordinates to az/el/r
	%[az_L,el_L,r_L]=cart2sph(P_L(n_and_g_noMW_combined_L,1),P_L(n_and_g_noMW_combined_L,2),P_L(n_and_g_noMW_combined_L,3));
	%[az_R,el_R,r_R]=cart2sph(P_R(n_and_g_noMW_combined_R,1),P_R(n_and_g_noMW_combined_R,2),P_R(n_and_g_noMW_combined_R,3));
	% get spherical coordinates (az/el/r, r equal in sphere) relevant to this network
	az_L_n=az_L(InclLeft);
	el_L_n=el_L(InclLeft);
	az_R_n=az_R(InclRight);
	el_R_n=el_R(InclRight);
	
	% now same mask for the opfl vectors
	OpF_azes_L_n=OpF_azes_L(InclLeft,:);
	OpF_els_L_n=OpF_els_L(InclLeft,:);
	OpF_azes_R_n=OpF_azes_R(InclRight,:);
	OpF_els_R_n=OpF_els_R(InclRight,:);

	% translate xyz vector components at coordinates to az/el/r
	nazes_L=zeros(length(az_L_n),1);
	nels_L=zeros(length(el_L_n),1);
	for i=1:length(az_L_n)
	    nvs_L=cart2sphvec(double([nx_L(i);ny_L(i);nz_L(i)]),az_L_n(i),el_L_n(i));
	    nazes_L(i)=nvs_L(1);
	    nels_L(i)=nvs_L(2);
	end
	% right hemi
	nazes_R=zeros(length(az_R_n),1);
	nels_R=zeros(length(el_R_n),1);
	for i=1:length(az_R_n)
	    nvs_R=cart2sphvec(double([nx_R(i);ny_R(i);nz_R(i)]),az_R_n(i),el_R_n(i));
	    nazes_R(i)=nvs_R(1);
	    nels_R(i)=nvs_R(2);
	end

	% initialize angular distance vector for each network (l and r) above
	NangDs_L=zeros(length(InclLeft),lenOpFl);
	NangDs_R=zeros(length(InclRight),lenOpFl);
	% get angular distance for each face for each timepoint
	for F=1:length(InclLeft);
		% get vector for each face (network vector)
		nVec=[nazes_L(F) nels_L(F)];
		% loop over each tp
		for fr=1:lenOpFl
			% get optical flow vector
			OpFlVec=[OpF_azes_L_n(F,fr) OpF_els_L_n(F,fr)];
			% get angular distance at that timepoint (degrees)
			a = acosd(min(1,max(-1, nVec(:).' *OpFlVec(:) / norm(nVec) / norm(OpFlVec) )));
			% populate vector
			NangDs_L(F,fr)=a;
		% end tp loop
		end
	% end each face loop
	end
	for F=1:length(InclRight);
                % get vector for each face (network vector)
                nVec=[nazes_R(F) nels_R(F)];
                % loop over each tp
                for fr=1:lenOpFl
                        % get optical flow vector
                        OpFlVec=[OpF_azes_R_n(F,fr) OpF_els_R_n(F,fr)];
                        % get angular distance at that timepoint (degrees)
                        a = acosd(min(1,max(-1, nVec(:).' *OpFlVec(:) / norm(nVec) / norm(OpFlVec) )));
                        % populate vector
                        NangDs_R(F,fr)=a;
                % end tp loop
                end
        % end each face loop
        end	
	% average for this network before proceeding to next network loop
	AllAngs=[NangDs_R(:)' NangDs_L(:)'];
	% average angular distances across hemispheres
	avgD=mean(AllAngs);
	Propvec=[Propvec avgD];
	% add label - adding 16 to be consistent with 
	stringVec=[stringVec ['AngD' num2str(k+16)]];
end
% save out as csv
T=table(Propvec','RowNames',stringVec);
% calc outFP
outFP=['/oak/stanford/groups/leanew1/users/apines/data/gp/PropFeats/' subj];
% make out filepath
system(['mkdir ' outFP]);
% write out
writetable(T,[outFP '/Prop_Feats.csv'],'WriteRowNames',true)
