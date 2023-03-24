function concat_TS(subj)
% load in subjs
topleveldir='/scratch/users/apines/abcd_images/derivatives/abcd-hcp-pipeline/sub-*'
direc=dir(topleveldir);
% deep gmless cifti for template
sname=subj;
parentfp=['/scratch/users/apines/abcd_images/fmriresults01/derivatives/abcd-hcp-pipeline/' sname '/ses-baselineYear1Arm1/func/'];
rsfp=[parentfp sname '_ses-baselineYear1Arm1_task-rest_p2mm_masked.dtseries.nii'];
sstfp=[parentfp sname '_ses-baselineYear1Arm1_task-SST_p2mm_masked.dtseries.nii'];
nbackfp=[parentfp sname '_ses-baselineYear1Arm1_task-nback_p2mm_masked.dtseries.nii'];
midfp=[parentfp sname '_ses-baselineYear1Arm1_task-MID_p2mm_masked.dtseries.nii'];

% use missing data directory
missingDir=['/oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/MissingDataReports/' sname];

missingFile=[missingDir '/MissingData.txt'];

if ~exist(sstfp,'file')
	disp('missing sst')
	system(['echo sst_missing >> ' missingFile]); 
end

if ~exist(rsfp,'file')
	disp('missing rs')
	system(['echo rs_missing >> ' missingFile]);
end

if ~exist(nbackfp,'file')
	disp('missing nback')
	system(['echo nback_missing >> ' missingFile]);
end

if ~exist(midfp,'file')
	disp('missing mid')
	system(['echo mid_missing >> ' missingFile]);
end


% isolate the masked grayordinatewise time series
% and stack them onto an init oordinate-wise value
Oords=zeros(91282,1);

if exist(rsfp,'file')
	rs=read_cifti(rsfp);
	rsts=rs.cdata;
	Oords=[Oords rsts];
	% initialize output cifti format on this too
	concat=rs;
end

if exist(sstfp,'file')
	sst=read_cifti(sstfp);
	sstts=sst.cdata;
	Oords=[Oords sstts];
	% init on this in case no rs
	concat=sst;
end

if exist(nbackfp,'file')
	nback=read_cifti(nbackfp);
	nbackts=nback.cdata;
	Oords=[Oords nbackts];
	% init on this in case no rs or sst
	concat=nback;
end

if exist(midfp,'file')
	mid=read_cifti(midfp);
	midfpts=mid.cdata;
	Oords=[Oords midfpts];
	% init on this in case it is the only scan
	concat=mid;
end

% remove initialization column of oords
Oords(:,1)=[];

concat.cdata=Oords;
% alter diminfo in cifit to align with new dims
csize=size(Oords);
concat.diminfo{2}.length=csize(2);
fp=[parentfp sname '_ses-baselineYear1Arm1_p2mm_masked_concat.dtseries.nii'];
write_cifti(concat,fp);
end
