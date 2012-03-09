function wls_bglm_estimate = nirs_run_liom_GLM_estimate_cfg
[NIRSmat redo1 NIRSmatCopyChoice] = get_common_NIRSmat(1,'GLM');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIOM General Linear Model Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NIRS_SPM_which_GLM      = cfg_menu;
NIRS_SPM_which_GLM.tag  = 'NIRS_SPM_which_GLM';
NIRS_SPM_which_GLM.name = 'Which GLM to estimate';
NIRS_SPM_which_GLM.labels = {'First','All', 'Last'};
NIRS_SPM_which_GLM.values = {1,2,3};
NIRS_SPM_which_GLM.val = {1};
NIRS_SPM_which_GLM.help = {'Choose which GLM (if more than one available) to estimate.'};

% Executable Branch
wls_bglm_estimate      = cfg_exbranch;
wls_bglm_estimate.name = 'LIOM GLM Estimation';
wls_bglm_estimate.tag  = 'wls_bglm_estimate';
wls_bglm_estimate.val  = {NIRSmat redo1 NIRSmatCopyChoice NIRS_SPM_which_GLM};
wls_bglm_estimate.prog = @nirs_run_liom_GLM_estimate;
wls_bglm_estimate.vout = @nirs_cfg_vout_liom_GLM_estimate;
wls_bglm_estimate.help = {'LIOM GLM Estimation: WLS (wavelet least square)'
    'and Bayesian GLM.'}';

function vout = nirs_cfg_vout_liom_GLM_estimate(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'NIRS.mat';
vout.src_output = substruct('.','NIRSmat');
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});