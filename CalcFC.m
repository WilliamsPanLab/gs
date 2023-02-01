function CalcFC(subcortTS,CortTS_L,CortTS_R,percyParcelParent,outFP)
%%% ∆∆∆∆
% calculate FC across multiple scales using personalized network boundaries and subcortical boundaries
%%% ∆∆∆∆

% cifti reading path
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/cifti-matlab/'));

% initialize (k+32)x(k+32)x29 matrix for multiscale fc

% 32x32 for pure subcort first
% 32 x k for subcort x cort
% k x k following pines 2022 for cort x cort
% consider saving out as 3 vectors?

% load in time series
LTS
RTS

% load in subcort time series
subcort=cifti_read(subcortTS);
% make subcort connectivity matrix right off the bat
% loop over each K
for K=2:30
% load in this partition
% initialize an N by TRs matrix for average TS within parcels over time
% loop over each N
for N=1:K
% get average TS for this N
% (or is a loop neccessary? can you vectorize?)
% end for each N
end
% corrmat N's + subcort

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
