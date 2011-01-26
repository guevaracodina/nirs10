function nirs_LIOM_specify_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Sets the defaults for NIRS coregistration

global nirs10

nirs10.model_specify.wls_bglm_specify.LiomDeleteLarge = 0;
%nirs10.model_specify.wls_bglm_specify.wls_or_bglm = 3;
nirs10.model_specify.wls_bglm_specify.channel_pca = 0;
nirs10.model_specify.wls_bglm_specify.lpf_butter.lpf_butter_On.lpf_butter_freq = 1/1.5;
nirs10.model_specify.wls_bglm_specify.hpf_butter.hpf_butter_On.hpf_butter_freq = 0.01;
nirs10.model_specify.wls_bglm_specify.volt = {2};
nirs10.model_specify.wls_bglm_specify.derivs = {[0 0]};
nirs10.model_specify.wls_bglm_specify.wls_or_bglm.NIRS_SPM.nirs_noise = 0;
nirs10.model_specify.wls_bglm_specify.time_res = {1};
nirs10.model_specify.wls_bglm_specify.units = {1};
nirs10.model_specify.wls_bglm_specify.time_res = {1};
%nirs10.model_specify.wls_bglm_specify.NIRS_SPM.hpf_wavelet_iter = {4};
%nirs10.model_specify.wls_bglm_specify.NIRS_SPM.hpf_dct_cutoff = {128};
%nirs10.model_specify.wls_bglm_specify.nirs_hpf = {hpf_none};
%nirs10.model_specify.wls_bglm_specify.NIRS_SPM.fwhm1 = {1.5};
%nirs10.model_specify.wls_bglm_specify.nirs_lpf = {lpf_none};
nirs10.model_specify.wls_bglm_specify.GLM_include_cardiac = 1;
nirs10.model_specify.wls_bglm_specify.GLM_include_Mayer = 0;
nirs10.model_specify.wls_bglm_specify.wls_or_bglm.WLS.WLS_J0 = 7;
nirs10.model_specify.wls_bglm_specify.wls_or_bglm.WLS.WLS_threshold_drift = 0.2;
nirs10.model_specify.wls_bglm_specify.wls_or_bglm.WLS.WLS_L0 = 0;
nirs10.model_specify.wls_bglm_specify.wls_or_bglm.BGLM.BGLM_fmax = 0.15;
nirs10.model_specify.wls_bglm_specify.wls_or_bglm.BGLM.BGLM_degre = 2;
nirs10.model_specify.wls_bglm_specify.wls_or_bglm.BGLM.BGLM_threshold_drift = 0.1;
