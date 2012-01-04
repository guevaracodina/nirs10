function checkreconstruct1 = nirs_run_checkreconstruct_cfg
NIRSmat         = cfg_files; %Select NIRS.mat for this subject
NIRSmat.name    = 'NIRS.mat'; % The displayed name
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [1 Inf];     % Number of inputs required
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vox_tempcours         = cfg_entry;
vox_tempcours.name    = 'Voxel';
vox_tempcours.tag     = 'vox_tempcours';
vox_tempcours.strtype = 'r';
vox_tempcours.num     = [0 3];
vox_tempcours.help    = {'Enter position otherwise click on the image to choose to voxel.'}';

roi_tempcours         = cfg_files;
roi_tempcours.name    = 'ROI';
roi_tempcours.tag     = 'roi_tempcours';
roi_tempcours.filter  = 'image';
roi_tempcours.ufilter = '.nii';
roi_tempcours.num     = [1 1];
roi_tempcours.help    = {'Mask with activation.'};

tempcours         = cfg_choice;
tempcours.name    = 'Temporal course';
tempcours.tag     = 'tempcours';
tempcours.values  = {vox_tempcours roi_tempcours};
tempcours.val     = {vox_tempcours};
tempcours.help    = {'Choose voxel or ROI.'};

t_stats      = cfg_branch;
t_stats.name = 'T stats of HbO and HbR in the whole volume';
t_stats.tag  = 't_stats';
t_stats.help = {''};

to_do         = cfg_choice;
to_do.tag     = 'to_do';
to_do.name    = 'Checking...';
to_do.values  = {t_stats tempcours};
to_do.val     = {tempcours};
to_do.help    = {''
    '-- Tikhonov is the simplest method of regularization'
    '-- Extended Tikhonov uses Li et al. model'}';

outreconstruct_Hb         = cfg_files;
outreconstruct_Hb.name    = 'Select Hb output files';
outreconstruct_Hb.tag     = 'outreconstruct_Hb';
outreconstruct_Hb.ufilter = {'.nii'};
outreconstruct_Hb.num     = [1 Inf];
outreconstruct_Hb.help    = {'.'};

% Executable Branch
checkreconstruct1      = cfg_exbranch;
checkreconstruct1.name = 'Check reconstruction';
checkreconstruct1.tag  = 'checkreconstruct1';
checkreconstruct1.val  = {NIRSmat outreconstruct_Hb to_do};
checkreconstruct1.prog = @nirs_run_checkreconstruct;
checkreconstruct1.vout = @nirs_cfg_vout_checkreconstruct;
checkreconstruct1.help = {'Check reconstruction.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_checkreconstruct(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
