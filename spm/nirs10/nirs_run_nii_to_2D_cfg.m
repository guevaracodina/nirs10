function convert_nii_to_2Dtopo = nirs_run_nii_to_2D_cfg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Utilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

convert_files         = cfg_files; %Select NIRS.mat for this subject
convert_files.name    = 'Select files to convert'; % The displayed name
convert_files.tag     = 'convert_files';       %file names
convert_files.filter  = 'image';
convert_files.ufilter = '.*';
convert_files.num     = [0 Inf];     % Number of inputs required
convert_files.help    = {'Select files to convert.'}';

convert_folders         = cfg_files; %Select NIRS.mat for this subject
convert_folders.name    = 'Select folders to convert'; % The displayed name
convert_folders.tag     = 'convert_folders';       %file names
convert_folders.filter = 'dir';
convert_folders.ufilter = '.*';
convert_folders.num     = [0 Inf];     % Number of inputs required
convert_folders.help    = {'Select folders of files to convert.'
    'Attempt will be made to convert all the content of specified folders.'}';

% Executable Branch
convert_nii_to_2Dtopo      = cfg_exbranch;
convert_nii_to_2Dtopo.name = 'Convert 2D nifti files to 2D topo files';
convert_nii_to_2Dtopo.tag  = 'convert_nii_to_2Dtopo';
convert_nii_to_2Dtopo.val  = {convert_files convert_folders};
convert_nii_to_2Dtopo.prog = @nirs_run_nii_to_2D;
convert_nii_to_2Dtopo.vout = @nirs_cfg_vout_nii_to_2D;
convert_nii_to_2Dtopo.help = {'Convert 2D topographich nifti (.nii) files '
    'to 2D topographic files (NIRS_SPM format).'}';

%make NIRS.mat available as a dependency -- not coded up here
function vout = nirs_cfg_vout_nii_to_2D(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
