function liom_HDM = nirs_run_liom_HDM_cfg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hemodynamic modeling - updated Feb. 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nameHDM         = cfg_entry;
nameHDM.name    = 'HDM Name';
nameHDM.tag     = 'nameHDM';
nameHDM.strtype = 's';
nameHDM.num     = [0 Inf];
nameHDM.val     = {''};
nameHDM.help    = {'Enter name for this HDM calculation.'}';

% % % Data % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
which_subjects_bold         = cfg_files; %
which_subjects_bold.name    = 'Select BOLD SPM folders for each subject'; 
which_subjects_bold.tag     = 'which_subjects_bold';       
which_subjects_bold.filter = 'dir';
which_subjects_bold.ufilter = '.*';
which_subjects_bold.num     = [1 Inf];    
which_subjects_bold.help    = {'Select folders for each subject containing'
    'first level GLM SPM analysis (SPM.mat) for BOLD data.'}';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
which_subjects_flow         = cfg_files; %
which_subjects_flow.name    = 'Select flow SPM folders for each subject'; 
which_subjects_flow.tag     = 'which_subjects_flow';       
which_subjects_flow.filter = 'dir';
which_subjects_flow.ufilter = '.*';
which_subjects_flow.num     = [1 Inf];    
which_subjects_flow.help    = {'Select folders for each subject containing'
    'first level GLM SPM analysis (SPM.mat) for flow data.'}';

% % % Modality % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xSPM_BOLD           = cfg_branch;
xSPM_BOLD.name      = 'BOLD only';
xSPM_BOLD.tag       = 'xSPM_BOLD';
xSPM_BOLD.val       = {which_subjects_bold};
xSPM_BOLD.help      = {''};

xSPM_ASL           = cfg_branch;
xSPM_ASL.name      = 'Flow only';
xSPM_ASL.tag       = 'xSPM_ASL';
xSPM_ASL.val       = {which_subjects_flow};
xSPM_ASL.help      = {''};

xSPM_BOLD_ASL           = cfg_branch;
xSPM_BOLD_ASL.name      = 'BOLD and flow from stat maps';
xSPM_BOLD_ASL.tag       = 'xSPM_BOLD_ASL';
xSPM_BOLD_ASL.val       = {which_subjects_bold which_subjects_flow}; 
xSPM_BOLD_ASL.help      = {'Not coded yet'};

xSPM_BOLD_ASL_V2           = cfg_branch;
xSPM_BOLD_ASL_V2.name      = 'BOLD and flow from ASL data';
xSPM_BOLD_ASL_V2.tag       = 'xSPM_BOLD_ASL_V2';
xSPM_BOLD_ASL_V2.val       = {which_subjects_bold which_subjects_flow};
xSPM_BOLD_ASL_V2.help      = {''}';

xSPM_Modalities           = cfg_choice;
xSPM_Modalities.name      = 'Data';
xSPM_Modalities.tag       = 'xSPM_Modalities';
xSPM_Modalities.values    = {xSPM_BOLD xSPM_BOLD_ASL xSPM_BOLD_ASL_V2 xSPM_ASL};
xSPM_Modalities.val       = {xSPM_BOLD};
xSPM_Modalities.help      = {'Choose data type: BOLD, BOLD & flow, flow only'};

Model_Choice = nirs_dfg_Model_Choice;

% % % Data otions % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
which_session = nirs_dfg_which_session;
which_condition = nirs_dfg_which_condition;

echo_time         = cfg_entry;
echo_time.tag     = 'echo_time';
echo_time.name    = 'Echo time (s)';
echo_time.help    = {'For BOLD, enter echo time in seconds. For BOLD+ASL, enter echo time of BOLD only'};
echo_time.strtype = 'e';
echo_time.num     = [1 1];
echo_time.val     = {0.030};

% % % ROIs % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genericROI = nirs_dfg_3D_roi;

% % % Options % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dp_start         = cfg_entry;
dp_start.tag     = 'dp_start';
dp_start.name    = 'Data points to remove (start)';
dp_start.help    = {'Enter number of data points to remove at the beginning of each file (to remove artefacts due to filtering)'};
dp_start.strtype = 'e';
dp_start.num     = [1 1];
dp_start.val     = {2};

dp_end         = cfg_entry;
dp_end.tag     = 'dp_end';
dp_end.name    = 'Data points to remove (end)';
dp_end.help    = {'Enter number of data points to remove at the end of each file (to remove artefacts due to filtering)'};
dp_end.strtype = 'e';
dp_end.num     = [1 1];
dp_end.val     = {2};

removeWhitening           = cfg_menu;
removeWhitening.name      = 'Remove Whitening Filter';
removeWhitening.tag       = 'removeWhitening';
removeWhitening.labels    = {'Yes' 'No'};
removeWhitening.values    = {1,0};
removeWhitening.val       = {0};
removeWhitening.help      = {'Remove SPM whitening filter on BOLD data prior to extracting VOIs'}';

priorFile         = cfg_files; 
priorFile.name    = 'Select spm_hdm_priors file to set covariance'; % The displayed name
priorFile.tag     = 'priorFile';       %file names
priorFile.filter = 'm';
priorFile.ufilter = '^spm_hdm_priors.*\.m$';
[tmpdir tmpfil tmpext] = fileparts(which('spm'));
theFile = fullfile(tmpdir,'spm_hdm_priors.m');
priorFile.num     = [0 1];     % Number of inputs required
priorFile.val     = {''};% BUG?????? {theFile};
priorFile.help    = {'Select spm_hdm_priors.m file used to set prior covariances on parameters.' ...
    'If omitted, the default spm values will be used (spm_hdm_priors.m).'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EM parameters
EM_parameters = nirs_dfg_hdm_EM;
[generate_figures save_figures] = nirs_dfg_generate_figures;
%Simulation options
simuOn = nirs_dfg_hdm_simu_options(1);

% Executable Branch
liom_HDM      = cfg_exbranch;
liom_HDM.name = 'LIOM Hemodynamic Modelling';
liom_HDM.tag  = 'liom_HDM';
liom_HDM.val  = {nameHDM xSPM_Modalities Model_Choice ...
     echo_time which_condition which_session  genericROI ...
    generate_figures save_figures dp_start dp_end ...
    priorFile removeWhitening EM_parameters ...
    simuOn};
liom_HDM.prog = @nirs_run_liom_HDM;
liom_HDM.vout = @nirs_cfg_vout_liom_HDM;
liom_HDM.help = {'NIRS_SPM Hemodynamic Modeling.'};

function vout = nirs_cfg_vout_liom_HDM(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});