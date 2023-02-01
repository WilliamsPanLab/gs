function CalcFC(subcortTS,CortTS_L,CortTS_R,subj,outFP)
%%% ∆∆∆∆
% calculate FC across multiple scales using personalized network boundaries and subcortical boundaries
%%% ∆∆∆∆
% general paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% cifti reading path
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/cifti-matlab/'));

%%%% ∆∆∆ initialize FC as vectors for speed
%%%% subcortical-subcortical
subcortFCsize=32*32;
% remove diagonal
subcortFCsize=subcortFCsize-32;
% remove redundant matrix reflection
subcortFCsize=subcortFCsize/2
% initialize
sub2sub=zeros(subcortFCsize);
% accomp. labels
sub2sub_labs=zeros(subcortFCsize);
%%%% for subcortical-cortical, use full number of networks x 32
subcortCortFCsize=32*464;
% no diagonal, but remove redundant reflection
subcortCortFCsize=subcortCortFCsize/2
sub2cort=zeros(subcortCortFCsize);
% accomp. labels
sub2cort_labs=zeros(subcortCortFCsize);
%%%% for cortical-cortical, use approach from pines 2022
cort2cort=zeros(4495);
cort2cort_labs=zeros(4495);

%%%% ∆∆∆ load in and format time series
dataSub=cifti_read(subcortTS);
dataL=MRIread(CortTS_L).vol;
dataR=MRIread(CortTS_R).vol;
% squeeze to get rid of extra dimensions
TRs_l_g=squeeze(dataL);
TRs_r_g=squeeze(dataR);
% calc time series length
timeSeriesLength=size(TRs_l_g,2);
% Load in SNR/mw mask 
surfML = '/oak/stanford/groups/leanew1/users/apines/data/gp/lh.Mask_SNR.label';
surfMR = '/oak/stanford/groups/leanew1/users/apines/data/gp/rh.Mask_SNR.label';
mwIndVec_l = read_medial_wall_label(surfML);
mwIndVec_r = read_medial_wall_label(surfMR);
% difference of full and invalid is valid vertices
Index_l = setdiff([1:10242], mwIndVec_l);
Index_r = setdiff([1:10242], mwIndVec_r);
% mask and stack left and right time series
ts_l_masked=TRs_l_g(logical(Index_l),:);
ts_r_masked=TRs_r_g(logical(Index_r),:);
ts_both=vertcat(ts_l_masked,ts_r_masked)';

%%%% ∆∆∆ make subcort connectivity matrix right off the bat
subCortCon=corr(dataSub.cdata');
% make upper triangle mask
UTM=triu(true(32),1);
%vectorize subcort FC
sub2sub=subCortCon(UTM);

%%%% ∆∆∆ loop over each K to extract cortical FCs
for K=2:30
	% load in this partition
	partitionFN=['/scratch/users/apines/derivatives/abcd-hcp-pipeline/SingleParcel_1by1/k' num2str(K) '/' subj '/IndividualParcel_Final_sbj1_comp' num2str(K) '_alphaS21_1_alphaL300_vxInfo1_ard0_eta0/final_UV.mat'];
	partition=load(partitionFN);
	% convert to hard parcels
	subj_V=partition.V{1};
	[~, HardParcel]=max(subj_V,[],2);
	% initialize an N by TRs matrix for average TS within parcels over time
	NetworkTimeSeries=zeros(K,timeSeriesLength);
	% loop over each N
	for N=1:K
		% get average TS for this N
		NetworkTimeSeries(N,:)=mean(ts_both(:,HardParcel==N),2);
	end
	% consider adding k index here!!!!!!!
	% cortical subcortical correlation matrix
	sub2cort()=corr(dataSub.cdata',NetworkTimeSeries');
	% Triu
	% cortical cortical correlation matrix
	cort2cort()=corr(NetworkTimeSeries');
	% Triu

% corrmat cortex-cortex

% extract vectorized diagonals

% maintain an equal length vector with k

% maintain an equal length vector with n (within k)

% saveout (or insert out) k label, n label, cor-subcort, cort-cort

% end for each k
end
% saveout subcort-subcort
% saveout cort-subcort
% saveout cort-cort
% save out reference vectors?

% end!
