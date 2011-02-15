function nirs_LIOM_specify_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Sets the defaults for NIRS coregistration

global nirs10

nirs10.model_specify.wls_bglm_specify.LiomDeleteLarge = 0;
nirs10.model_specify.wls_bglm_specify.channel_pca = 0;
nirs10.model_specify.wls_bglm_specify.lpf_butter.lpf_butter_On.lpf_butter_freq = 0.667; %1/1.5s
nirs10.model_specify.wls_bglm_specify.hpf_butter.hpf_butter_On.hpf_butter_freq = 0.01;
nirs10.model_specify.wls_bglm_specify.volt = {2};
nirs10.model_specify.wls_bglm_specify.derivs = {[0 0]};
nirs10.model_specify.wls_bglm_specify.wls_or_bglm.NIRS_SPM.nirs_noise = 0;
nirs10.model_specify.wls_bglm_specify.time_res = {10};
nirs10.model_specify.wls_bglm_specify.units = {1};
nirs10.model_specify.wls_bglm_specify.GLM_include_cardiac = 1;
nirs10.model_specify.wls_bglm_specify.GLM_include_Mayer = 0;
nirs10.model_specify.wls_bglm_specify.wls_or_bglm.WLS.WLS_J0 = 7;
nirs10.model_specify.wls_bglm_specify.wls_or_bglm.WLS.WLS_threshold_drift = 0.1;
nirs10.model_specify.wls_bglm_specify.wls_or_bglm.WLS.WLS_L0 = 0;
nirs10.model_specify.wls_bglm_specify.wls_or_bglm.BGLM.BGLM_fmax = 0.15;
nirs10.model_specify.wls_bglm_specify.wls_or_bglm.BGLM.BGLM_degre = 2;
nirs10.model_specify.wls_bglm_specify.wls_or_bglm.BGLM.BGLM_threshold_drift = 0.1;
nirs10.model_specify.wls_bglm_specify.GenerateHbT = 1;
nirs10.model_specify.wls_bglm_specify.flag_window = 0;