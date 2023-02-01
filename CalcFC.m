function CalcFC(subcortTS,CortTS_L,CortTS_R,subj,outFP)
%%% ∆∆∆∆
% calculate FC across multiple scales using personalized network boundaries and subcortical boundaries
%%% ∆∆∆∆
% general paths
%%% OMIT FOR MASTER CALL TO REDUCE DILY-DALYING ON GENPATHS
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% cifti reading path
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/cifti-matlab/'));

%%%% ∆∆∆ initialize FC as vectors for speed
%%%% subcortical-subcortical
subcortFCsize=32*32;
% remove diagonal
subcortFCsize=subcortFCsize-32;
% remove redundant matrix reflection
subcortFCsize=subcortFCsize/2;
% initialize
sub2sub=zeros(1,subcortFCsize);
%%%% Cut
% accomp. labels not needed: same 1:32 ordering as OG
% ¬ sub2sub_labs=zeros(subcortFCsize);
%%%% for subcortical-cortical, use full number of networks x 32
% ¬ subcortCortFCsize=32*464;
% no diagonal and no redundancies in row/colnames: stick with OG size
% ¬ sub2cort=zeros(1,subcortCortFCsize);
%%%% Cut
sub2cort=zeros(32,464);
% accomp. labels
sub2cort_labs=zeros(1,464);
%%%% for cortical-cortical, use approach from pines 2022
cort2cort=zeros(1,4495);
cort2cort_labsA=zeros(1,4495);
cort2cort_labsB=zeros(1,4495);

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

%%%% ∆∆∆ Make and index of where each scale should go in network-wise vector
Kind_corres={};
% make an index of which places in network vec align with which scale
Krange=2:30;
for K=Krange
	K_start=((K-1)*(K))/2;
	K_end=(((K-1)*(K))/2)+K-1;
	Kind_corres{K}=K_start:K_end;
end
% and an interator for edges
edgeIter=1;
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
	% indicate where network-level features begin for this scale with triangular numbers formula
	K_start=((K-1)*(K))/2;
	% same, but for end
	K_end=(((K-1)*(K))/2)+K-1;
	% loop over each N
	for N=1:K
		% get average TS for this N
		NetworkTimeSeries(N,:)=mean(ts_both(:,HardParcel==N),2);
	end
	% cortical subcortical correlation matrix: port in in column representing first network at this scale to last network at this scale Kstart to KEnd
	sub2cort(:,K_start:K_end)=corr(dataSub.cdata',NetworkTimeSeries');
	sub2cort_labs(K_start:K_end)=K_start:K_end;	
	% cortical cortical correlation matrix: coincidentally (?) the number of edges at this scale is equal to K_start (triangular number for this scale)
	lastEdge=edgeIter+K_start-1;
	% mask out upper triangle
	Edges=corr(NetworkTimeSeries');
	UTM=triu(true(K),1);
	% plop into df
	cort2cort(edgeIter:lastEdge)=Edges(UTM);
	% to keep track of which network is node A and which network is node B, plop row and column represented in each vector in parallel with edge features
	NodeListRows=K_start:K_end;
	NodeListCols=K_start:K_end;
	[Rows,Cols] = meshgrid(NodeListRows,NodeListCols);
	cort2cort_labsA(edgeIter:lastEdge)=Rows(UTM);
	cort2cort_labsB(edgeIter:lastEdge)=Cols(UTM);
	% update edge iterator
	edgeIter=lastEdge+1;
%%%
end
% saveout subcort-subcort
csvwrite([outFP '_S-S.csv'],sub2sub)
% saveout cort-subcort
csvwrite([outFP '_S-C.csv'],sub2cort)
% saveout cort-cort
csvwrite([outFP '_C-C.csv'],cort2cort)
% save out reference vectors on one run for reference
csvwrite([outFP '_S-C_labs.csv'],sub2cort_labs)
csvwrite([outFP '_C-C_labsA.csv'],cort2cort_labsA)
csvwrite([outFP '_C-C_labsB.csv'],cort2cort_labsB)
