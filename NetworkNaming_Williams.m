% set K Range for multiscale
Krange=2:30;
% add needed paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));
% williams 6 circuits
LH_circuits=gifti('~/groHP_left.func.gii')
RH_circuits=gifti('~/groHP_right.func.gii')
%%%% COMBINE EM AND MASK
% Load in SNR/mw mask 
surfML = '/oak/stanford/groups/leanew1/users/apines/data/gp/lh.Mask_SNR.label';
surfMR = '/oak/stanford/groups/leanew1/users/apines/data/gp/rh.Mask_SNR.label';
mwIndVec_l = read_medial_wall_label(surfML);
mwIndVec_r = read_medial_wall_label(surfMR);
% difference of full and invalid is valid vertices
Index_l = setdiff([1:10242], mwIndVec_l);
Index_r = setdiff([1:10242], mwIndVec_r);
% mask and stack left and right time series
LH_circuits=LH_circuits.cdata(logical(Index_l),1);
RH_circuits=RH_circuits.cdata(logical(Index_r),1);
W6_Label = [LH_circuits' RH_circuits'];
% circuit names
CircuitNames={'DMN', 'CC', 'Atn', 'Rew', 'Sal', 'Thr', 'Other'};
% create a vector that spans all K's, length should be sum of all Ks (2+3+4...+30)
correspondence_over_scales=zeros((length(Krange)*((min(Krange)+max(Krange))/2)),1);
% Make and index of where each scale should go in transmodality vector (largely just to confirm everything is in its right place, in its right place, in its right place, right place)
Kind_corres={};
% make an index of which places in feature vec align with which scale
for K=Krange
	K_start=((K-1)*(K))/2;
	K_end=(((K-1)*(K))/2)+K-1;
	Kind_corres{K}=K_start:K_end;
end
% initialize output DF
df_corres=cell(4,length(correspondence_over_scales));
% make strings vector for colnames
corres_strings=strings(length(correspondence_over_scales),1);
% loop over each scale
for K=Krange
	K
	% get indices for this K
	Kind=Kind_corres{K};
	for N=1:K
		curindex=Kind(N);
		corres_strings(curindex)=strcat('Corres_w6_scores_scale_', num2str(K), '_net', num2str(N));
	end
	% plop in the scale into ouput df
	df_corres(2,Kind)=deal(num2cell(K));
	% plop in the network at this scale (col name) into output df
	df_corres(1,Kind)=cellstr(corres_strings(Kind));
	% get in parcellation for this scale (meaning start the K-loop above here)
	Parcellation_Path = ['/oak/stanford/groups/leanew1/users/apines/data/init_Vs/' num2str(K) '/init_UV.mat'];
	Parcellation_Mat = load(Parcellation_Path);
	Parcellation_Lab=Parcellation_Mat.V{:}; 
	[~, HardParcel]=max(Parcellation_Lab,[],2);
	% loop over hard parcels to get overlap with circuits
	for N=1:K
		% get vertices where this network is predominant
		Index = find(HardParcel == N);
		% get williams 6 membership on these vertices 
		Label_W6=W6_Label(Index);
		% count vertex membership in each
		number=[];
		for W=1:6
			% number of vertices in N that are in this W
			number(W)=length(find(Label_W6==W));
		end
		number
		% extract maximum represented W6
		[~, Max_Index] = max(number);
		% display network name
		disp([num2str(N) ' Max W6 alignment:  ' CircuitNames{Max_Index}]);
		% input this information into output df
		df_corres(3,Kind(N))=cellstr(CircuitNames{Max_Index});
		% calculate proportion of network falling within W6 circuit
		totverts=sum(number);
		prop=(number(Max_Index))/totverts;
		% along with propotion
		df_corres(4,Kind(N))=num2cell(prop);
		disp([num2str(N) ' Max W6 alignment proportion:  ' num2str(prop)]);
	end
end
% saveout
writetable(cell2table(df_corres),'/oak/stanford/groups/leanew1/users/apines/data/gp/networks_wCorrespondence.csv');

