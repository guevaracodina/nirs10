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

% % % Modality % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xSPM_BOLD           = cfg_branch;
xSPM_BOLD.name      = 'BOLD only';
xSPM_BOLD.tag       = 'xSPM_BOLD';
xSPM_BOLD.val       = {};
xSPM_BOLD.help      = {''};

xSPM_ASL           = cfg_branch;
xSPM_ASL.name      = 'Flow only';
xSPM_ASL.tag       = 'xSPM_ASL';
xSPM_ASL.val       = {};
xSPM_ASL.help      = {''};

xSPM_BOLD_ASL           = cfg_branch;
xSPM_BOLD_ASL.name      = 'BOLD and flow from stat maps';
xSPM_BOLD_ASL.tag       = 'xSPM_BOLD_ASL';
xSPM_BOLD_ASL.val       = {}; 
xSPM_BOLD_ASL.help      = {'Not coded yet'};

xSPM_BOLD_ASL_V2           = cfg_branch;
xSPM_BOLD_ASL_V2.name      = 'BOLD and flow from ASL data';
xSPM_BOLD_ASL_V2.tag       = 'xSPM_BOLD_ASL_V2';
xSPM_BOLD_ASL_V2.val       = {};
xSPM_BOLD_ASL_V2.help      = {'!!! Location and selection of flow and BOLD '
    'functional data is hard-coded for Michele''s project !!!'}';

xSPM_Modalities           = cfg_choice;
xSPM_Modalities.name      = 'Modalities';
xSPM_Modalities.tag       = 'xSPM_Modalities';
xSPM_Modalities.values    = {xSPM_BOLD xSPM_BOLD_ASL xSPM_BOLD_ASL_V2 xSPM_ASL};
xSPM_Modalities.val       = {xSPM_BOLD};
xSPM_Modalities.help      = {'Choose data type: BOLD, BOLD + flow, flow only'};


Model_Choice           = cfg_menu;
Model_Choice.name      = 'Model choice';
Model_Choice.tag       = 'Model_Choice';
Model_Choice.labels    = {'Buxton-Friston' 'Zheng-Mayhew' 'Boas-Huppert'};
Model_Choice.values    = {0,1,2};
Model_Choice.val       = {0};
%Model_Choice.def  = @(val)nirs_get_defaults('readNIRS.boxy1.save_bin_dot', val{:});
Model_Choice.help      = {'Choose hemodynamic model: Buxton-Friston, '
    'Zheng-Mayhew, or 1-Compartment Boas-Huppert Model'}';



% % % Data % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

which_subjects         = cfg_files; %
which_subjects.name    = 'Select SPM folders for each subject'; 
which_subjects.tag     = 'which_subjects';       
which_subjects.filter = 'dir';
which_subjects.ufilter = '.*';
which_subjects.num     = [1 Inf];    
which_subjects.help    = {'Select folders for each subject containing'
    'first level GLM SPM analysis (UR0 or UR3 folders).'}';

which_session     = cfg_entry;
which_session.name    = 'Which session(s)?';
which_session.tag     = 'which_session';
which_session.strtype = 'r';
which_session.val     = {1};
which_session.num     = [0 1];
which_session.help    = {'Enter session numbers (based on those included in SPM.mat)'}';

which_condition     = cfg_entry;
which_condition.name    = 'Which condition(s)?';
which_condition.tag     = 'which_condition';
which_condition.strtype = 'r';
which_condition.val     = {1};
which_condition.num     = [0 Inf];
which_condition.help    = {'Enter stimuli numbers to include as a Matlab row '
    'vector (get from the design matrix associated with the data file)'}';

echo_time         = cfg_entry;
echo_time.tag     = 'echo_time';
echo_time.name    = 'Echo time (s)';
echo_time.help    = {'For BOLD, enter echo time in seconds. For BOLD+ASL, enter echo time of BOLD only'};
echo_time.strtype = 'e';
echo_time.num     = [1 1];
echo_time.val     = {0.030};


% % % ROIs % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% % % Options % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_figures           = cfg_menu;
save_figures.name      = 'Save figures';
save_figures.tag       = 'save_figures';
save_figures.labels    = {'Yes' 'No'};
save_figures.values    = {1,0};
save_figures.val       = {1};
save_figures.help      = {'Save figures.'}';

generate_figures      = cfg_menu;
generate_figures.tag  = 'generate_figures';
generate_figures.name = 'Show figures';
generate_figures.labels = {'Yes','No'};
generate_figures.values = {1,0};
generate_figures.val  = {0};
generate_figures.help = {'Show figures. When selecting this option, the figures will stay opened after the code has completed.'}';

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


% % % Simulations % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simuA         = cfg_entry;
simuA.tag     = 'simuA';
simuA.name    = '% Signal amplitude';
simuA.help    = {'Enter signal amplitude, as a percentage of the signal (e.g. enter 1 for a 1% amplitude)'};
simuA.strtype = 'e';
simuA.num     = [1 1];
simuA.val     = {100};

simuS         = cfg_entry;
simuS.tag     = 'simuS';
simuS.name    = 'which_condition to simulate';
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

simuIt         = cfg_entry;
simuIt.tag     = 'simuIt';
simuIt.name    = 'Number of simulations';
simuIt.help    = {'Enter number of simulations'};
simuIt.strtype = 'e';
simuIt.num     = [1 1];
simuIt.val     = {1};

noiseNo             = cfg_branch;
noiseNo.tag         = 'noiseNo';
noiseNo.name        = 'No noise in simulated data';
noiseNo.val         = {}; 
noiseNo.help        = {'Simulated data will be noiseless (0 baseline).'};


% % % Experimental noise % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

noiseYes             = cfg_branch;
noiseYes.tag         = 'noiseYes';
noiseYes.name        = 'Add noise to simulated data';
noiseYes.val         = {restscans restscans_BOLD restscans_ASL}; 
noiseYes.help        = {['Noise will be added to simulated data before inversion.'...
    ' This noise is added in the form of resting state data measured experimentally.' ...
    ' Specify those experimental scans below.']};

simuNoise           = cfg_choice;
simuNoise.name      = 'Include baseline noise';
simuNoise.tag       = 'simuNoise';
simuNoise.values    = {noiseYes noiseNo};
simuNoise.val       = {noiseYes};
simuNoise.help      = {'Include noise background scans; if No, the simulated data will be noiseless, i.e. on 0 background.'}';

simuUpsample         = cfg_entry;
simuUpsample.tag     = 'simuUpsample';
simuUpsample.name    = 'Data upsampling factor';
simuUpsample.help    = {'Enter an upsampling factor (max = 16)'};
simuUpsample.strtype = 'e';
simuUpsample.num     = [1 1];
simuUpsample.val     = {1};

simuInterp         = cfg_entry;
simuInterp.tag     = 'simuInterp';
simuInterp.name    = 'Data interpolation factor';
simuInterp.help    = {'Enter an interpolation factor for simulated data. Data will first be '...
    ' simulated at 16 time bins per TR, then decimated to TR, or between 1-1/16 of TR '...
    ' if an upsampling factor is specified. Then, noise will be added if specified. ' ...
    ' Only then, data will be interpolated by this factor. Large interpolation factors will increase '...
    ' computation time significantly.' };
simuInterp.strtype = 'e';
simuInterp.num     = [1 1];
simuInterp.val     = {1};

simuR         = cfg_entry;
simuR.tag     = 'simuR';
simuR.name    = 'Parameter range';
simuR.help    = {'For each parameter specified, enter range to be sampled as a percentage of the prior value.'
    'Parameters will be sampled randomly uniformly between (100-range) and (100+range) % of the prior value.'
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


distr_uniform           = cfg_branch;
distr_uniform.tag       = 'distr_uniform';
distr_uniform.name      = 'Uniform';
distr_uniform.val       = {simuPrior simuR};
distr_uniform.help      = {'Simulated priors drawn pseudo-randomly from a uniform distribution.'};

simuMean1         = cfg_entry;
simuMean1.tag     = 'simuMean1';
simuMean1.name    = 'Mean for Gaussian #1';
simuMean1.help    = {'Enter array of prior values of the previously specified parameters to be sampled.'
    'If nothing is entered, default prior values will be used.'};
simuMean1.strtype = 'e';
simuMean1.num     = [0 Inf];
simuMean1.val     = {''};

% SimuMean2         = cfg_entry;
% SimuMean2.tag     = 'SimuMean2';
% SimuMean2.name    = 'Mean for Gaussian #2';
% SimuMean2.help    = {'Enter array of prior values of the previously specified parameters to be sampled.'
%     'If nothing is entered, default prior values will be used.'};
% SimuMean2.strtype = 'e';
% SimuMean2.num     = [0 Inf];
% SimuMean2.val     = {''};

simuMean21         = cfg_entry;
simuMean21.tag     = 'simuMean21';
simuMean21.name    = '% Difference between Gaussian means';
simuMean21.help    = {['Enter array, for all parameters, of % values that ' ...
    ' specify the difference between the means of the 2 Gaussians distributions.' ...
    ' For example, if you specify "25" for one parameter, the first half of parameters are drawn ' ...
    'from a Gaussian distribution with mean M specified above, and the second half ' ...
    'from a Gaussian with 1.25*M.' ...
    'If nothing is entered, all parameters will be drawn from 1 Gaussian distribution.']};
simuMean21.strtype = 'e';
simuMean21.num     = [0 Inf];
simuMean21.val     = {''};

simuR1         = cfg_entry;
simuR1.tag     = 'simuR1';
simuR1.name    = 'Parameter range #1';
simuR1.help    = {'For each parameter specified, enter range to be sampled as a percentage of the prior value.'
    'Parameters will be sampled from a gaussian distribution with sigma = this range.'
    'If only one number is specified, it will be applied to all parameters to be sampled.'}';
simuR1.strtype = 'e';
simuR1.num     = [1 Inf];
simuR1.val     = {25};

simuR2         = cfg_entry;
simuR2.tag     = 'simuR2';
simuR2.name    = 'Parameter range #2';
simuR2.help    = {'For each parameter specified, enter range to be sampled as a percentage of the prior value.'
    'Parameters will be sampled from a gaussian distribution with sigma = this range.'
    'If only one number is specified, it will be applied to all parameters to be sampled.'}';
simuR2.strtype = 'e';
simuR2.num     = [1 Inf];
simuR2.val     = {25};

distr_bimodal           = cfg_branch;
distr_bimodal.tag       = 'distr_bimodal';
distr_bimodal.name      = '2 Gaussians';
distr_bimodal.val       = {simuMean1 simuMean21 simuR1 simuR2};
distr_bimodal.help      = {['Half of simulated priors drawn pseudo-randomly '...
    'each from a gaussian distribution.']};

simuParamDistr         = cfg_choice;
simuParamDistr.tag     = 'simuParamDistr';
simuParamDistr.name    = 'Parameter distribution';
simuParamDistr.values  = {distr_uniform distr_bimodal};
simuParamDistr.val     = {distr_uniform};
simuParamDistr.help  = {['Type of distribution from which are pseudo-randomly' ...
    'drawn the parameters used to simulate data (priors).']};

simuYes         = cfg_branch;
simuYes.tag     = 'simuYes';
simuYes.name    = 'HDM on simulated data';
simuYes.val     = {simuIt simuA simuS simuP simuParamDistr simuUpsample simuInterp simuNoise};
simuYes.help    = {'Perform HDM on simulated data'}';

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EM parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Niterations         = cfg_entry; 
Niterations.name    = 'Maximum number of EM iterations';
Niterations.tag     = 'Niterations';       
Niterations.strtype = 'r';
Niterations.num     = [1 1];     
Niterations.val     = {128};
Niterations.help    = {'Maximum number of EM iterations. 128 is the basic number.'
    'Increase to 512 or more if necessary.'}';

dFcriterion         = cfg_entry; 
dFcriterion.name    = 'Convergence criterion';
dFcriterion.tag     = 'dFcriterion';       
dFcriterion.strtype = 'r';
dFcriterion.num     = [1 1];     
dFcriterion.val     = {1};
dFcriterion.help    = {'Convergence criterion on changes of free energy F.'
    'Changes in F less than this value are required for convergence in E and in M steps.'
    '1e-2 is the SPM8 default value.'}';

LogAscentRate         = cfg_entry; 
LogAscentRate.name    = 'Initial log ascent rate';
LogAscentRate.tag     = 'LogAscentRate';       
LogAscentRate.strtype = 'r';
LogAscentRate.num     = [1 1];     
LogAscentRate.val     = {-2};
LogAscentRate.help    = {'Initial log ascent rate: control initial rate of movement in parameter space'}';

MaxLogAscentRate         = cfg_entry; 
MaxLogAscentRate.name    = 'Maximal log ascent rate';
MaxLogAscentRate.tag     = 'MaxLogAscentRate';       
MaxLogAscentRate.strtype = 'r';
MaxLogAscentRate.num     = [1 1];     
MaxLogAscentRate.val     = {4};
MaxLogAscentRate.help    = {'Maximal absolute value of log ascent rate: control minimal/maximal rate of movement in parameter space'}';

spm_integrator      = cfg_menu;
spm_integrator.tag  = 'spm_integrator';
spm_integrator.name = 'Choose ODE integrator';
spm_integrator.labels = {'spm_int','spm_int_ode','spm_int_J'};
spm_integrator.values = {'spm_int','spm_int_ode','spm_int_J'};
spm_integrator.val  = {'spm_int'};
spm_integrator.help = {'Choose integrator to use for ordinary differential equations.'
    'spm_int is fastest, spm_int_ode most precise, spm_int_J in between'}';

Mstep_iterations         = cfg_entry; 
Mstep_iterations.name    = 'Maximum number of M step iterations';
Mstep_iterations.tag     = 'Mstep_iterations';       
Mstep_iterations.strtype = 'r';
Mstep_iterations.num     = [1 1];     
Mstep_iterations.val     = {8};
Mstep_iterations.help    = {'Maximum number of M-step iterations. 8 is the standard choice.'}';

EM_parameters         = cfg_branch;
EM_parameters.tag     = 'EM_parameters';
EM_parameters.name    = 'Parameters of EM algorithm';
EM_parameters.val     = {Niterations spm_integrator dFcriterion LogAscentRate MaxLogAscentRate Mstep_iterations}; 
EM_parameters.help    = {'Parameters of Expectation Maximization algorithm.'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Executable Branch
liom_HDM      = cfg_exbranch;
liom_HDM.name = 'LIOM Hemodynamic Modelling';
liom_HDM.tag  = 'liom_HDM';
liom_HDM.val  = {nameHDM xSPM_Modalities Model_Choice which_subjects ...
     echo_time which_condition which_session  genericROI ...
    generate_figures save_figures  ...
    dp_start dp_end ...
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