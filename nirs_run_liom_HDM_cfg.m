function liom_HDM = nirs_run_liom_HDM_cfg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hemodynamic modeling - new version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xSPM_spmmat         = cfg_files; %Select NIRS.mat for this subject
xSPM_spmmat.name    = 'Select xSPM.mat (for BOLD)'; % The displayed name
xSPM_spmmat.tag     = 'xSPM_spmmat';       %file names
xSPM_spmmat.filter = 'mat';
xSPM_spmmat.ufilter = '^xSPM';
xSPM_spmmat.num     = [0 1];     % Number of inputs required
xSPM_spmmat.help    = {'Select xSPM.mat of BOLD estimation for this subject or group.'};

xSPM_spmmat_ASL         = cfg_files;
xSPM_spmmat_ASL.name    = 'Select xSPM.mat (for ASL)'; % The displayed name
xSPM_spmmat_ASL.tag     = 'xSPM_spmmat_ASL';       %file names
xSPM_spmmat_ASL.filter = 'mat';
xSPM_spmmat_ASL.ufilter = '^xSPM';
xSPM_spmmat_ASL.num     = [0 1];     % Number of inputs required
xSPM_spmmat_ASL.help    = {'Select xSPM.mat of ASL estimation for this subject or group.'};

% VOI_spmmat         = cfg_files; %Select NIRS.mat for this subject
% VOI_spmmat.name    = 'Select VOI_....mat (for BOLD)'; % The displayed name
% VOI_spmmat.tag     = 'VOI_spmmat';       %file names
% VOI_spmmat.filter = 'mat';
% VOI_spmmat.ufilter = '^VOI';
% VOI_spmmat.num     = [0 1];     % Number of inputs required
% VOI_spmmat.help    = {'Select VOI_....mat of BOLD estimation for this subject (or group).'};
%
% VOI_spmmat_ASL         = cfg_files;
% VOI_spmmat_ASL.name    = 'Select VOI_....mat (for ASL)'; % The displayed name
% VOI_spmmat_ASL.tag     = 'VOI_spmmat_ASL';       %file names
% VOI_spmmat_ASL.filter = 'mat';
% VOI_spmmat_ASL.ufilter = '^VOI';
% VOI_spmmat_ASL.num     = [0 1];     % Number of inputs required
% VOI_spmmat_ASL.help    = {'Select VOI_....mat of ASL estimation for this subject (or group).'};

which_subjects_ASL         = cfg_files; %Select NIRS.mat for this subject
which_subjects_ASL.name    = 'Select SPM folders for each subject'; % The displayed name
which_subjects_ASL.tag     = 'which_subjects_ASL';       %file names
which_subjects_ASL.filter = 'dir';
which_subjects_ASL.ufilter = '.*';
which_subjects_ASL.num     = [1 Inf];     % Number of inputs required
which_subjects_ASL.help    = {'Select folders for each subject containing'
    'first level GLM SPM analysis for ASL (UR1 folders).'}';


spmmat         = cfg_files; %Select NIRS.mat for this subject
spmmat.name    = 'Select SPM.mat (for BOLD)'; % The displayed name
spmmat.tag     = 'spmmat';       %file names
spmmat.filter = 'mat';
spmmat.ufilter = '^SPM.mat$';
spmmat.num     = [1 1];     % Number of inputs required
spmmat.help    = {'Select SPM.mat of BOLD estimation for this subject (or group).'};

spmmat_ASL         = cfg_files;
spmmat_ASL.name    = 'Select SPM.mat (for ASL)'; % The displayed name
spmmat_ASL.tag     = 'spmmat_ASL';       %file names
spmmat_ASL.filter = 'mat';
spmmat_ASL.ufilter = '^SPM.mat$';
spmmat_ASL.num     = [1 1];     % Number of inputs required
spmmat_ASL.help    = {'Select SPM.mat of ASL estimation for this subject (or group).'};


xSPM_BOLD           = cfg_branch;
xSPM_BOLD.name      = 'BOLD only';
xSPM_BOLD.tag       = 'xSPM_BOLD';
xSPM_BOLD.val       = {spmmat xSPM_spmmat}; % VOI_spmmat};
xSPM_BOLD.help      = {''};

xSPM_ASL           = cfg_branch;
xSPM_ASL.name      = 'ASL only';
xSPM_ASL.tag       = 'xSPM_ASL';
xSPM_ASL.val       = {spmmat_ASL xSPM_spmmat_ASL}; % VOI_spmmat_ASL};
xSPM_ASL.help      = {''};

xSPM_BOLD_ASL           = cfg_branch;
xSPM_BOLD_ASL.name      = 'BOLD and ASL from stat maps';
xSPM_BOLD_ASL.tag       = 'xSPM_BOLD_ASL';
xSPM_BOLD_ASL.val       = {spmmat spmmat_ASL xSPM_spmmat xSPM_spmmat_ASL which_subjects_ASL};
xSPM_BOLD_ASL.help      = {'Not coded yet'};

xSPM_BOLD_ASL_V2           = cfg_branch;
xSPM_BOLD_ASL_V2.name      = 'BOLD and ASL from flow and BOLD data';
xSPM_BOLD_ASL_V2.tag       = 'xSPM_BOLD_ASL_V2';
xSPM_BOLD_ASL_V2.val       = {spmmat xSPM_spmmat}; % flow_subj1 bold_subj1};
xSPM_BOLD_ASL_V2.help      = {'Location and selection of Flow and BOLD '
    'functional data is hard-coded for Michele''s Project'}';

xSPM_Modalities           = cfg_choice;
xSPM_Modalities.name      = 'Modalities: BOLD, BOLD + ASL estimation, ASL only';
xSPM_Modalities.tag       = 'xSPM_Modalities';
xSPM_Modalities.values    = {xSPM_BOLD xSPM_BOLD_ASL xSPM_BOLD_ASL_V2 xSPM_ASL};
xSPM_Modalities.val       = {xSPM_BOLD};
xSPM_Modalities.help      = {'Choose data type: BOLD, BOLD + ASL, ASL only'};

which_session     = cfg_entry;
which_session.name    = 'Which session?';
which_session.tag     = 'which_session';
which_session.strtype = 'r';
which_session.val     = {1};
which_session.num     = [0 1];
which_session.help    = {'Enter session to analyze'}';

which_subjects         = cfg_files; %Select NIRS.mat for this subject
which_subjects.name    = 'Select SPM folders for each subject'; % The displayed name
which_subjects.tag     = 'which_subjects';       %file names
which_subjects.filter = 'dir';
which_subjects.ufilter = '.*';
which_subjects.num     = [1 Inf];     % Number of inputs required
which_subjects.help    = {'Select folders for each subject containing'
    'first level GLM SPM analysis (UR0 or UR3 folders).'}';

nameROI         = cfg_entry;
nameROI.name    = 'ROI name';
nameROI.tag     = 'nameROI';
nameROI.strtype = 's';
nameROI.num     = [0 Inf];
nameROI.val     = {''};
nameROI.help    = {'Enter name for ROI. If left blank, ROIs will be enumerated.'}';

radiusROI         = cfg_entry;
radiusROI.name    = 'ROI radius value';
radiusROI.tag     = 'radiusROI';
radiusROI.strtype = 'r';
radiusROI.num     = [1 1];
radiusROI.val     = {5};
radiusROI.help    = {'Radius value in mm'}';

coordinateROI         = cfg_entry;
coordinateROI.name    = 'ROI coordinates';
coordinateROI.tag     = 'coordinateROI';
coordinateROI.strtype = 'r';
coordinateROI.num     = [1 3];
%coordinateROI.val     = {};
coordinateROI.help    = {'Enter MNI coordinates [x y z] in mm for center of ROI'}';

whichROI         = cfg_branch;
whichROI.tag     = 'whichROI';
whichROI.name    = 'Define ROI';
whichROI.val     = {nameROI coordinateROI radiusROI};
whichROI.help    = {'Define ROI'}';

genericROI         = cfg_repeat;
genericROI.tag     = 'genericROI';
genericROI.name    = 'Define ROIs';
genericROI.help    = {'Define here the ROIs to be analyzed'}';
genericROI.values  = {whichROI};
genericROI.num     = [1 Inf];

HDMdisplay           = cfg_menu;
HDMdisplay.name      = 'Display individual results';
HDMdisplay.tag       = 'HDMdisplay';
HDMdisplay.labels    = {'Yes' 'No'};
HDMdisplay.values    = {1,0};
HDMdisplay.val       = {0};
HDMdisplay.help      = {'Display output for each subject'}';

save_figures           = cfg_menu;
save_figures.name      = 'Save figures';
save_figures.tag       = 'save_figures';
save_figures.labels    = {'Yes' 'No'};
save_figures.values    = {1,0};
save_figures.val       = {1};
save_figures.help      = {'Save figures for each subject'}';

nameHDM         = cfg_entry;
nameHDM.name    = 'HDM Name';
nameHDM.tag     = 'nameHDM';
nameHDM.strtype = 's';
nameHDM.num     = [0 Inf];
nameHDM.val     = {''};
nameHDM.help    = {'Enter name for this HDM calculation.'}';

% StimuliSign     = cfg_entry;
% StimuliSign.name    = 'Stimuli sign';
% StimuliSign.tag     = 'Stimuli';
% StimuliSign.strtype = 'r';
% StimuliSign.val     = {1};
% StimuliSign.num     = [0 Inf];
% StimuliSign.help    = {'Enter sign of response to each stimulus, as a vector.'}';

echo_time         = cfg_entry;
echo_time.tag     = 'echo_time';
echo_time.name    = 'Echo time (s)';
echo_time.help    = {'For BOLD, enter echo time in seconds. For BOLD+ASL, enter echo time of BOLD only'};
echo_time.strtype = 'e';
echo_time.num     = [1 1];
echo_time.val     = {0.030};

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
removeWhitening.val       = {1};
removeWhitening.help      = {'Remove SPM whitening filter on BOLD data prior to extracting VOIs'}';

restscans         = cfg_files;
restscans.tag     = 'restscans';
restscans.name    = 'Resting state scans';
restscans.help    = {'Select resting state scans on which the simulated data will be added.'
    'The code was developed to work specifically on the 4D.nii scans in folder'
    '11-ep2d_bold_ax_4x4x4_TR_1010 of Michele subject 28'
    'Assumptions are that these scans are  compatible with the SPM and xSPM structures'
    '(same number number of scans in particular)'}';
restscans.filter = 'image';
restscans.ufilter = '.*';
restscans.val     = {''};
restscans.num     = [0 Inf];

restscans_BOLD         = cfg_files;
restscans_BOLD.tag     = 'restscans_BOLD';
restscans_BOLD.name    = 'Resting state scans for BOLD';
restscans_BOLD.help    = {'ONLY for simulating BOLD+ASL'}';
restscans_BOLD.filter = 'image';
restscans_BOLD.ufilter = '.*';
restscans_BOLD.val     = {''};
restscans_BOLD.num     = [0 Inf];

restscans_ASL         = cfg_files;
restscans_ASL.tag     = 'restscans_ASL';
restscans_ASL.name    = 'Resting state scans for ASL';
restscans_ASL.help    = {'ONLY for simulating BOLD+ASL'}';
restscans_ASL.filter = 'image';
restscans_ASL.ufilter = '.*';
restscans_ASL.val     = {''};
restscans_ASL.num     = [0 Inf];

simuA         = cfg_entry;
simuA.tag     = 'simuA';
simuA.name    = 'Signal Amplitude';
simuA.help    = {'Enter signal amplitude, as a percentage of the BOLD signal (e.g. enter 1 for a 1% amplitude)'};
simuA.strtype = 'e';
simuA.num     = [1 1];
simuA.val     = {1};

simuS         = cfg_entry;
simuS.tag     = 'simuS';
simuS.name    = 'Stimuli to simulate';
simuS.help    = {'Enter array of stimuli types to simulated.'
    'Enter 0 to include all stimuli types.'}';
simuS.strtype = 'e';
simuS.num     = [1 Inf];
simuS.val     = {0};

simuP         = cfg_entry;
simuP.tag     = 'simuP';
simuP.name    = 'Parameters to randomize';
simuP.help    = {'Enter array of parameters to be sampled.'
    'Enter 0 to randomize all parameters.'}';
simuP.strtype = 'e';
simuP.num     = [1 Inf];
simuP.val     = {1};

simuR         = cfg_entry;
simuR.tag     = 'simuR';
simuR.name    = 'Parameter range';
simuR.help    = {'For each parameter specified, enter range to be sampled as a percentage of the prior value.'
    'Parameters will be sampled randomly uniformly between 0.75 and 1.25 times the prior value.'
    'If only one number is specified, it will be applied to all parameters to be sampled.'}';
simuR.strtype = 'e';
simuR.num     = [1 Inf];
simuR.val     = {25};

simuPrior         = cfg_entry;
simuPrior.tag     = 'simuPrior';
simuPrior.name    = 'Prior values of parameters';
simuPrior.help    = {'Enter array of prior values of the previously specified parameters to be sampled.'
    'If nothing is entered, the default prior values will be used.'}';
simuPrior.strtype = 'e';
simuPrior.num     = [0 Inf];
simuPrior.val     = {''};

simuIt         = cfg_entry;
simuIt.tag     = 'simuIt';
simuIt.name    = 'Number of simulations';
simuIt.help    = {'Enter number of simulations'};
simuIt.strtype = 'e';
simuIt.num     = [1 1];
simuIt.val     = {1};

simuNoise           = cfg_menu;
simuNoise.name      = 'Include baseline noise';
simuNoise.tag       = 'simuNoise';
simuNoise.labels    = {'Yes' 'No'};
simuNoise.values    = {1,0};
simuNoise.val       = {1};
simuNoise.help      = {'Include noise background scans; if No, the simulated data will be noiseless, i.e. on 0 background.'}';

simuUpsample         = cfg_entry;
simuUpsample.tag     = 'simuUpsample';
simuUpsample.name    = 'Data upsampling factor';
simuUpsample.help    = {'Enter an upsampling factor (max = 16)'};
simuUpsample.strtype = 'e';
simuUpsample.num     = [1 1];
simuUpsample.val     = {1};

simuYes         = cfg_branch;
simuYes.tag     = 'simuYes';
simuYes.name    = 'HDM on simulated data';
simuYes.val     = {simuIt simuA simuS simuP simuPrior simuR simuUpsample simuNoise restscans restscans_BOLD restscans_ASL};
simuYes.help    = {'Perform HDM on real data'}';

simuNo         = cfg_branch;
simuNo.tag     = 'simuNo';
simuNo.name    = 'HDM on real data';
simuNo.val     = {};
simuNo.help    = {'No simulations; perform HDM on real data'}';

simuOn         = cfg_choice;
simuOn.tag     = 'simuOn';
simuOn.name    = 'Perform simulations';
simuOn.values  = {simuNo simuYes};
simuOn.val     = {simuNo};
simuOn.help    = {'Perform simulations'}';


Model_Choice           = cfg_menu;
Model_Choice.name      = 'Choice of Model';
Model_Choice.tag       = 'Model_Choice';
Model_Choice.labels    = {'Buxton-Friston' 'Zheng-Mayhew' 'Boas-Huppert'};
Model_Choice.values    = {0,1,2};
Model_Choice.val       = {0};
%Model_Choice.def  = @(val)nirs_get_defaults('readNIRS.boxy1.save_bin_dot', val{:});
Model_Choice.help      = {'Choose hemodynamic model: Buxton-Friston, '
    'Zheng-Mayhew, or 1-Compartment Boas-Huppert Model'}';

Stimuli     = cfg_entry;
Stimuli.name    = 'Stimuli identification numbers';
Stimuli.tag     = 'Stimuli';
Stimuli.strtype = 'r';
Stimuli.val     = {1};
Stimuli.num     = [0 Inf];
Stimuli.help    = {'Enter stimuli numbers to include as a Matlab row '
    'vector (get from the design matrix associated with the data file)'}';

% Executable Branch
liom_HDM      = cfg_exbranch;
liom_HDM.name = 'LIOM Hemodynamic Modelling';
liom_HDM.tag  = 'liom_HDM';
liom_HDM.val  = {xSPM_Modalities Model_Choice Stimuli which_session ...
    which_subjects genericROI HDMdisplay save_figures nameHDM echo_time ...
    dp_start dp_end removeWhitening simuOn};
liom_HDM.prog = @nirs_run_liom_HDM;
liom_HDM.vout = @nirs_cfg_vout_liom_HDM;
liom_HDM.help = {'NIRS_SPM Hemodynamic Modeling.'};

function vout = nirs_cfg_vout_liom_HDM(job)
vout = cfg_dep;
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});