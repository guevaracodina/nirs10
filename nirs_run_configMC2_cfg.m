function configMC1 = nirs_run_configMC2_cfg


NIRSmat         = cfg_files; %Select NIRS.mat for this subject
NIRSmat.name    = 'NIRS.mat'; % The displayed name
NIRSmat.tag     = 'NIRSmat';       %file names
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';
NIRSmat.num     = [1 Inf];     % Number of inputs required
NIRSmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed

DelPreviousData      = cfg_menu;
DelPreviousData.tag  = 'DelPreviousData';
DelPreviousData.name = 'Delete Previous data file';
DelPreviousData.labels = {'True','False'};
DelPreviousData.values = {1,0};
DelPreviousData.val  = {0};
DelPreviousData.help = {'Delete the previous data file.'}';

CreateNIRSCopy_false         = cfg_branch;
CreateNIRSCopy_false.tag     = 'CreateNIRSCopy_false';
CreateNIRSCopy_false.name    = 'Do not copy NIRS structure';
CreateNIRSCopy_false.help    = {'Do not copy NIRS structure.'
    'This will write over the previous NIRS.mat'}';

NewNIRSdir         = cfg_entry;
NewNIRSdir.name    = 'Directory for NIRS.mat';
NewNIRSdir.tag     = 'NewNIRSdir';
NewNIRSdir.strtype = 's';
NewNIRSdir.val{1}    = 'NewDir';
NewNIRSdir.num     = [1 Inf];
NewNIRSdir.help    = {'Directory for NIRS.mat.'}';

CreateNIRSCopy         = cfg_branch;
CreateNIRSCopy.tag     = 'CreateNIRSCopy';
CreateNIRSCopy.name    = 'Create new directory and copy NIRS structure';
CreateNIRSCopy.val     = {NewNIRSdir};
CreateNIRSCopy.help    = {'Create new directory and copy NIRS structure there.'}';

%Common to most modules: for creating a new directory and copying NIRS.mat
NewDirCopyNIRS           = cfg_choice;
NewDirCopyNIRS.name      = 'Create new directory and copy NIRS.mat';
NewDirCopyNIRS.tag       = 'NewDirCopyNIRS';
NewDirCopyNIRS.values    = {CreateNIRSCopy_false CreateNIRSCopy};
NewDirCopyNIRS.val       = {CreateNIRSCopy_false};
NewDirCopyNIRS.help      = {'Choose whether to overwrite the NIRS.mat structure'
    'or to create a new directory'
    'and copy the NIRS.mat structure there'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configure input files for Monte Carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

latest_mcim         = cfg_entry;
latest_mcim.tag     = 'latest_mcim';
latest_mcim.name    = 'Latest ROI';
latest_mcim.val     = {'latestROI'};
latest_mcim.help    = {'Latest ROI generated selected.'};

mcim_in         = cfg_files;
mcim_in.name    = 'MC segmented volume';
mcim_in.tag     = 'mcim_in';
mcim_in.filter = 'image';
mcim_in.ufilter = '.nii';
mcim_in.num     = [1 1];
mcim_in.help    = {'Select MC segmented volume for this subject.'};

mcim_cfg           = cfg_choice;
mcim_cfg.name      = 'Image';
mcim_cfg.tag       = 'mcim_cfg';
mcim_cfg.values    = {latest_mcim mcim_in};
mcim_cfg.val       = {latest_mcim};
mcim_cfg.help      = {'Choose latest ROI simulated or specify the segmented image.'};

MC_CUDAchoice    = cfg_menu;
MC_CUDAchoice.name   = 'Configuration file type';
MC_CUDAchoice.tag    = 'MC_CUDAchoice';
MC_CUDAchoice.labels = {'MCX: Monte Carlo Extreme','tMCimg','Both'};
MC_CUDAchoice.values = {1,2,3};
MC_CUDAchoice.def    = @(val)nirs_get_defaults('configMC1.MC_CUDAchoice', val{:});
MC_CUDAchoice.help   = {'Choose type of configuration files to generate.'};

no_pve      = cfg_branch;
no_pve.name = 'No pertubation';
no_pve.tag  = 'no_pve';
no_pve.help = {'No a priori data about perturbation. PVE won''t be calculated.'};

% calc_pve         = cfg_files;
% calc_pve.name    = 'Mask for PVE';
% calc_pve.tag     = 'calc_pve';
% calc_pve.filter  = 'image';
% calc_pve.ufilter = '.nii';
% calc_pve.num     = [1 1];
% calc_pve.help    = {'BOLD, ASL or any anatomical mask (from create mask module : not coded yet).'};
%
pve_bold      = cfg_branch;
pve_bold.name = 'BOLD mask';
pve_bold.tag  = 'pve_bold';
pve_bold.help = {'PVE will be calculated with respect to last BOLD mask.'};

pve_asl      = cfg_branch;
pve_asl.name = 'ASL mask';
pve_asl.tag  = 'pve_asl';
pve_asl.help = {'PVE will be calculated with respect to last ASL mask.'};

pve_anat         = cfg_files;
pve_anat.name    = 'User-chosen mask';
pve_anat.tag     = 'pve_anat';
pve_anat.filter  = 'image';
pve_anat.ufilter = '.nii';
pve_anat.num     = [1 1];
pve_anat.help    = {['PVE will be calculated with respect to a user-chosen mask.' ...
    ' The image can contain any value and have any resolution different from ' ...
    'that of the simulation volume. It will be thresholded above 0 and transformed' ...
    ' to the space of the simulation medium.']};

pve_cfg           = cfg_choice;
pve_cfg.name      = 'Perturbation Mask';
pve_cfg.tag       = 'pve_cfg';
pve_cfg.values    = {no_pve pve_bold pve_asl pve_anat};
pve_cfg.val       = {pve_bold};
pve_cfg.help      = {'A priori data about spatial location of a perturbation. It will be used fro reconstructions and PVF in perturbation mask (BOLD or other) will be computed.'};

dpf_cfg        = cfg_menu;
dpf_cfg.name   = 'Evaluate DPF';
dpf_cfg.tag    = 'dpf_cfg';
dpf_cfg.labels = {'Yes','No'};
dpf_cfg.values = {1,0};
dpf_cfg.val{1} = 1;
dpf_cfg.help   = {'Evaluate DPF thanks to MonteCarlo simulation.'};

% est ce qu'il y a pas un pb du au fait qu'il attend un directory ???
MC_configdir         = cfg_entry;
MC_configdir.tag     = 'MC_configdir';
MC_configdir.name    = 'Monte Carlo configuration files directory';
MC_configdir.strtype = 's';
MC_configdir.num     = [1 Inf];
MC_configdir.def     = @(val)nirs_get_defaults('configMC1.MC_configdir', val{:});
MC_configdir.help    = {'Directory to put Monte Carlo configuration files.'
    'NO Longer USED'}';

MC_nam         = cfg_entry;
MC_nam.tag     = 'MC_nam';
MC_nam.name    = 'Monte Carlo simulation name';
MC_nam.strtype = 's';
MC_nam.num     = [1 Inf];
MC_nam.val{1}  = 'sim';
MC_nam.help    = {'Name of Monte Carlo simulation.'
    'If a simulation has already been run, '
    'the current date will be automatically added to the name specified'}';

mulitt      = cfg_branch;
mulitt.name = 'Litterature coefficients';
mulitt.tag  = 'mulitt';
mulitt.help = {'No a priori data about perturbation. PVE won''t be calculated.'};

muTRS         = cfg_files;
muTRS.name    = 'File containing subjects mua and mus';
muTRS.tag     = 'muTRS';
muTRS.ufilter = '.mat';
muTRS.num     = [1 1];
muTRS.help    = {'Values from TRS for exemple. DESCRIBE FILE'};

mu_subj           = cfg_choice;
mu_subj.name      = 'Use subject particular parameters';
mu_subj.tag       = 'mu_subj';
mu_subj.values    = {mulitt muTRS};
mu_subj.val       = {mulitt};
mu_subj.help      = {'Parameters measured thanks to TRS.'};


%--------------------------------------------------------------------------
nphotons         = cfg_entry;
nphotons.name    = 'Number of photons'; % The displayed name
nphotons.tag     = 'nphotons';       %file names
nphotons.strtype = 'r';
nphotons.num     = [1 1];     % Number of inputs required
nphotons.def = @(val)nirs_get_defaults('configMC1.nphotons', val{:});
nphotons.help    = {'Input number of photons (not currently used).'};

seed         = cfg_entry;
seed.name    = 'Random seed';
seed.tag     = 'seed';
seed.strtype = 'r';
seed.num     = [1 1];
seed.def = @(val)nirs_get_defaults('configMC1.seed', val{:});
seed.help    = {'Input random seed.'};

modulationFreq         = cfg_entry;
modulationFreq.name    = 'Modulation Frequency'; % The displayed name
modulationFreq.tag     = 'modulationFreq';       %file names
modulationFreq.strtype = 'r';
modulationFreq.num     = [1 1];     % Number of inputs required
modulationFreq.def = @(val)nirs_get_defaults('configMC1.modulationFreq', val{:});
modulationFreq.help    = {'Modulation Frequency; leave at 0 for CW (continuous wave operation).'};

deltaT         = cfg_entry;
deltaT.name    = 'deltaT'; % The displayed name
deltaT.tag     = 'deltaT';       %file names
deltaT.strtype = 'r';
deltaT.num     = [1 1];     % Number of inputs required
deltaT.def = @(val)nirs_get_defaults('configMC1.deltaT', val{:});
deltaT.help    = {'deltaT: interval between time gates in seconds.'};

numTimeGates         = cfg_entry;
numTimeGates.name    = 'numTimeGates'; % The displayed name
numTimeGates.tag     = 'numTimeGates';       %file names
numTimeGates.strtype = 'r';
numTimeGates.num     = [1 1];     % Number of inputs required
numTimeGates.def = @(val)nirs_get_defaults('configMC1.numTimeGates', val{:});
numTimeGates.help    = {'Number of time gates; total duration is number of time gates times deltaT.'};

radiis         = cfg_entry;
radiis.name    = 'Source Radii'; % The displayed name
radiis.tag     = 'radiis';       %file names
radiis.strtype = 'r';
radiis.num     = [1 1];     % Number of inputs required
radiis.def = @(val)nirs_get_defaults('configMC1.radiis', val{:});
radiis.help    = {'Input radius of the sources.'};

radiid         = cfg_entry;
radiid.name    = 'Detector radii'; % The displayed name
radiid.tag     = 'radiid';       %file names
radiid.strtype = 'r';
radiid.num     = [1 1];     % Number of inputs required
radiid.def = @(val)nirs_get_defaults('configMC1.radiid', val{:});
radiid.help    = {'Input radius of the detectors.'};

voxelSize         = cfg_entry;
voxelSize.name    = 'Voxel Size'; % The displayed name
voxelSize.tag     = 'voxelSize';       %file names
voxelSize.strtype = 'r';
voxelSize.num     = [1 1];     % Number of inputs required
voxelSize.def = @(val)nirs_get_defaults('configMC1.voxelSize', val{:});
voxelSize.help    = {'Input voxel Size.'};

perturbationPpties_l1 = cfg_entry;
perturbationPpties_l1.name    = 'Perturbation first wavelength'; % The displayed name
perturbationPpties_l1.tag     = 'perturbationPpties_l1';       %file names
perturbationPpties_l1.strtype = 'r';
perturbationPpties_l1.num     = [1 4];     % Number of inputs required
perturbationPpties_l1.def = @(val)nirs_get_defaults('configMC1.perturbationPpties_l1', val{:});
perturbationPpties_l1.help    = {['Action on grey matter only: ',...
    'Perturbation properties Delta(\mu_a,\mu_s, g, n) for first wavelength (default = 830 nm).']};

perturbationPpties_l2 = cfg_entry;
perturbationPpties_l2.name    = 'Perturbation second wavelength'; % The displayed name
perturbationPpties_l2.tag     = 'perturbationPpties_l2';       %file names
perturbationPpties_l2.strtype = 'r';
perturbationPpties_l2.num     = [1 4];     % Number of inputs required
perturbationPpties_l2.def = @(val)nirs_get_defaults('configMC1.perturbationPpties_l2', val{:});
perturbationPpties_l2.help    = {['Action on grey matter only: ',...
    'Perturbation properties Delta(\mu_a,\mu_s, g, n) for first wavelength (default = 830 nm).']};

MC_parameters      = cfg_branch;
MC_parameters.tag  = 'MC_parameters';
MC_parameters.name = 'Parameters';
MC_parameters.val  = {nphotons seed modulationFreq deltaT numTimeGates radiis radiid voxelSize ...
    perturbationPpties_l1 perturbationPpties_l2};
MC_parameters.help = {'Parameters'};

% Executable Branch
configMC1      = cfg_exbranch;
configMC1.name = 'Configure Monte Carlo inputs';
configMC1.tag  = 'configMC1';
configMC1.val  = {NIRSmat MC_nam mcim_cfg MC_CUDAchoice dpf_cfg pve_cfg MC_configdir mu_subj MC_parameters};
configMC1.prog = @nirs_run_configMC2;
configMC1.vout = @nirs_cfg_vout_configMC;
configMC1.help = {'Generate configuration input files for Monte Carlo simulation.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_configMC(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});