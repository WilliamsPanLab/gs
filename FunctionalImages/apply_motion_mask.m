function apply_motion_mask(subj)
% Initialize qc features to saveout
FDvec=zeros(1,4);
FDstringVec={'restFD','MIDFD','SSTFD','nbackFD'};
Missvec=zeros(1,7);
MissstringVec={'restM','MIDM','SSTM','nbackM','CTM','MMM','SCVM'};
RemTRs=0;
Remstring={'RemainingTRs'};
% initiailize empty vector for average length
TRvecNum=[];
% for each "task"
tasks=["rest","MID","SST","nback"];
for t=1:4
	task=tasks(t);
	sname=subj;
	fpParent=['/scratch/users/apines/abcd_images/fmriresults01/derivatives/abcd-hcp-pipeline/' sname '/ses-baselineYear1Arm1/func/'];
	fp=strjoin([fpParent sname '_ses-baselineYear1Arm1_task-' task '_bold_desc-filtered_timeseries.dtseries.nii'],'');
	% not flagging missing tasks in this script, added this conditional to reflect that
	if exist(fp,'file')
	ts_cif=read_cifti(fp);
	ts=ts_cif.cdata;
	% load in mask
	masfp=strjoin([fpParent sname '_ses-baselineYear1Arm1_task-' task '_desc-filtered_motion_mask.mat'],'');
	if exist(masfp,'file')
	mask=load(masfp);
	% get to FD_thresh of .2 mm, corresponds to threshold 21
	maskp2mm=mask.motion_data{1,21}.frame_removal;
	TRwise_mask=logical(maskp2mm);
	% length of mask corresponds to number of TRs
	% 1 indicates flagged for FD over selected threshold, reverse 'em so 0 indicates removal
	TRwise_mask=~TRwise_mask;
	% remove TRs with corresp. flag
	masked_trs=ts(:,TRwise_mask);
	% reconfig cifti metadata to reflect new number of TRs
	newciftiSize=size(masked_trs);
	newTRnum=newciftiSize(2);
	%s-2 because we start at 3
	ts_cif.diminfo{2}.length=newTRnum;
	% overwrite TRs for new file
	ts_cif.cdata=masked_trs;
	% add to remaining TRs
	RemTRs=RemTRs+newTRnum;
	% add remaining FD
	FDvec(t)=mask.motion_data{1,21}.remaining_frame_mean_FD;
	% set output filepath
	ofp=strjoin([fpParent sname '_ses-baselineYear1Arm1_task-' task '_p2mm_masked.dtseries.nii'],'');
	% There is no reason this should be a requried step
	ofp=convertStringsToChars(ofp);
	% write out motion masked cifti
	write_cifti(ts_cif,ofp);
	% write out format_string for use in optical flow segmentation later
	% note the meaning of this mask: https://github.com/DCAN-Labs/dcan_bold_processing/blob/eda67d1c1eed2b31e6a14781d7ab075561314dd4/matlab_code/framewise_displacement/format_generator.m
	% initial 7 are 1s, subsequent x's indicate include the following "v" volumes (or exlcude the following from the motion mask) conversely, +'s indicate exclude the following v volumes (or inc. in motmask)
	% numbers are length of sequence, which is culmatively additive. NOT what is starting TR w/r/t global sequence
	% it also includes between-run frames, which aligns with what we want to exclude
	writelines(mask.motion_data{1,21}.format_string,strjoin([fpParent 'MotionSegments_' task '.txt'],''))
	else
	missingDir=['/oak/stanford/groups/leanew1/users/apines/scripts/abcdImages/MissingDataReports/' sname];
	mkdir(missingDir);
	missingFile=[missingDir '/MissingData.txt'];
	system(['echo motionMask_missing >> ' missingFile]);
	% add to missing vec
	Missvec(t)=1;
	end
	else
	end
end
% check for missing structural data to keep QC feats + missingness in same writes
CTFP=['/scratch/users/apines/abcd_images/imagingcollection01/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/anat/' subj '_ses-baselineYear1Arm1_space-fsLR32k_thickness.dscalar.nii'];
if ~exist(CTFP,'file')
	Missvec(5)=1;
else
end
MMFP=['/scratch/users/apines/abcd_images/imagingcollection01/derivatives/abcd-hcp-pipeline/' subj '/ses-baselineYear1Arm1/anat/' subj '_ses-baselineYear1Arm1_space-fsLR32k_desc-smoothed_myelinmap.dscalar.nii'];
if ~exist(MMFP,'file')
	Missvec(6)=1;
else
end
SCVFP=['/scratch/users/apines/abcd_images/imagingcollection01/derivatives/freesurfer-5.3.0-HCP/' subj '/ses-baselineYear1Arm1/stats/aseg.stats'];
if ~exist(SCVFP,'file')
	Missvec(7)=1;
else
end
% save out QC features
qcDir=['/oak/stanford/groups/leanew1/users/apines/data/gp/QC_Feats/' sname];
mkdir(qcDir)
ValuesVec=[FDvec Missvec RemTRs];
StringsVec=[FDstringVec MissstringVec Remstring];
T=table(ValuesVec','RowNames',StringsVec);
writetable(T,[qcDir '/QC_Feats.csv'],'WriteRowNames',true)

