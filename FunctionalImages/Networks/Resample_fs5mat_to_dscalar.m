function Resample_fs5mat_to_dscalar(subj)

% load in init consensus map
initConsensus=load('/oak/stanford/groups/leanew1/users/apines/data/init_Vs/18/init_UV.mat');
Loading_Mat=load(initConsensus);
V=Loading_Mat.V;
% extract from struct
V=V{:};

% save out as dscalar (but get fs5 map?)
HP=cifti_read('/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii')
HP.cdata(1:59412)=sbj_AtlasLabel_NoMedialWall;
outputfile=['/cbica/projects/abcdfnets/results/SingleParcel_1by1/' subj '/' subj '_Parcel.dscalar.nii'];
cifti_write(HP,outputfile)

% loop over every k and save out sep. func.gii map
for i = 1:17;
	% load file
	fp=[unmasked_GC_folder 'Group_AtlasLoading_Network_' num2str(i) '.dscalar.nii'];
	NetLoads=cifti_read(fp);
	% extract LH
	lhLoad=NetLoads.cdata(1:10242);
	% initialize gifti
	LH_gif=gifti;
	LH_gif.cdata=lhLoad;
	V_lh_File = [unmasked_GC_folder 'Group_lh_Network_' num2str(i) '.func.gii'];
	% save
	save(LH_gif, V_lh_File);
	% extract RH
	rhLoad=NetLoads.cdata(10243:20484);
	% initialize gifti
	RH_gif=gifti;
	RH_gif.cdata=rhLoad;
	V_rh_File = [unmasked_GC_folder 'Group_rh_Network_' num2str(i) '.func.gii'];
	% save
	save(RH_gif, V_rh_File);
	i
end

