function us = OpFl_abcd(subj,task,tsIn_L,tsIn_R,tsOut)
% convert nomenclature
sname=subj;
% add paths
%%%%%%%%%%%%%%%%%%%% Set parameters
N = 10; % Vector spherical harmonics degree
h = 1; % finite-difference height vector of triangles
alpha = 1; % scales Tikhonov regularization, > alpha = > spatiotemporal smoothness
s = 1; % R(u), regularizing functional, scales Tikhonov regularization more rapidly via penalization of first spatial and first temporal derivative
%%%%%%%%%%%%%%%%%%%%

% load in data
fpL=tsIn_L;
fpR=tsIn_R;
dataL=MRIread(fpL).vol;
dataR=MRIread(fpR).vol;
% squeeze to get rid of extra dimensions
TRs_l_g=squeeze(dataL);
TRs_r_g=squeeze(dataR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
% for surface data
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% normalize verts to unit sphere
% left
numV=length(vx_l);
vx_l(numV+1:end, :) = VecNormalize(vx_l(numV+1:end, :));
% right
numV=length(vx_r);
vx_r(numV+1:end, :) = VecNormalize(vx_r(numV+1:end, :));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake functional data (on surface)
% handle input data
disp('Size of input data')
sizeInDl=size(TRs_l_g)
sizeInDr=size(TRs_r_g)
disp('Number of TRs detected L and R')
TR_n=sizeInDl(2)
sizeInDr(2)
assert(sizeInDl(2) == sizeInDr(2), 'Unequal time series length between hemispheres')

% left hemi
disp('converting left hemi to struct')
fl=struct;
% populate struct
for TRP=1:TR_n;
	fl.TRs{TRP}=TRs_l_g(:,TRP);
end

% r h
disp('converting right hemi to struct')
fr=struct;
for TRP=1:TR_n;
	fr.TRs{TRP}=TRs_r_g(:,TRP);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load in continuous segment indices
% load in motin mask and use dcan/midb's string indices containing contiguous segments
fpParent=['/scratch/users/apines/abcd_images/fmriresults01/derivatives/abcd-hcp-pipeline/' sname '/ses-baselineYear1Arm1/func/'];
masfp=[fpParent sname '_ses-baselineYear1Arm1_task-' task '_desc-filtered_motion_mask.mat'];
MotMask=load(masfp);
% consistent with masking script
MotStruct=MotMask.motion_data{21};
% extract string
MotString=MotStruct.format_string;
% convert string to nx2 table where n is the number of contiguous segments, column 1 ('Var1') is the starting TR for the segment, and Var2 is the length of the contiguous segment
iterator=1;
Var1=[];
Var2=[];
% extract the good and the bad segments
segments = regexp(MotString, '\d+[x+]', 'match');
nsegs=length(segments);
% loop through each segment
for n=1:nsegs
	currSeg=segments{n};
	% if +, included segment
	if contains(currSeg,'+')
		% drop in current starting TR in masked sequence
		Var1=[Var1 iterator];
		% retrieve number by removing operator
		currSeg=currSeg(1:end-1);
		% record length of this segment 
		Var2=[Var2 str2num(currSeg)];
		% add to iterator
		iterator=iterator+str2num(currSeg);
	else % else means it was an x, meaning it was cut
	end
end
% convert variables to continuous segment indices as prior
CSI=table(Var1',Var2');
% extract numbers of TRs
segLengths=strsplit(MotString,{'x','+'});
% assure that TR count is the same between time series and valid segments txt
SegNum=height(CSI);
% trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
numTRsVS=CSI{SegNum,1}+CSI{SegNum,2}-1;
if numTRsVS ~= TR_n
        disp('TRs from Valid Segments txt and cifti do not match. Fix it.')
        return
end

% remove any segments <5 TRs that didn't get thresholded in midb pipeline
mask=CSI.Var2>4;
CSI=CSI(mask,:);
% recompute segNum
SegNum=height(CSI);

% addpath for OF to not confuse height function earlier on
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute optical flow on every pair of sequential TRs
% initialize output struct
us=struct;
disp('Computing optical flow זְרִימָה  प्रवाहः  دفق');

% initialize TRP counter: for plopping u outputs into master struct w/o/r/t their segment
% note trp = tr pair
TRPC=1;

% for each segment
for seg=1:SegNum;
	% loop over each TR-Pair: 1 fewer pair than number of TRs
	for TRP=1:(TR_n-1)
		% print TR pair iter
		TRP
		% Compute decomposition.
		% pull out adjacent frames
		u = of(N, faces_l, vx_l, fl.TRs{TRP}, fl.TRs{TRP+1}, h, alpha, s);
		% throw u into struct
		us.vf_left{TRPC}=u;
		% now right hemi
		u = of(N, faces_r, vx_r, fr.TRs{TRP}, fr.TRs{TRP+1}, h, alpha, s);
		% tihrow u into struct
		us.vf_right{TRPC}=u;
		% update TR pair counter, which should increase +1 across segments
		TRPC=TRPC+1;
	end
end
% save
save(tsOut,'us')

