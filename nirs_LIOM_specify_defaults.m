function nirs_LIOM_specify_defaults
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Sets the defaults for NIRS coregistration

global nirs10

nirs10.model_specify.wls_bglm_specify.LiomDeleteLarge = 2;
nirs10.model_specify.wls_bglm_specify.wls_or_bglm = 3;
nirs10.model_specify.wls_bglm_specify.channel_pca = 0;
nirs10.model_specify.wls_bglm_specify.lpf_butter = 0;
nirs10.model_specify.wls_bglm_specify.hpf_butter = 1;
nirs10.model_specify.wls_bglm_specify.volt = {2};
nirs10.model_specify.wls_bglm_specify.derivs = {[0 0]};
nirs10.model_specify.wls_bglm_specify.nirs_noise = {0};
nirs10.model_specify.wls_bglm_specify.time_res = {1};
nirs10.model_specify.wls_bglm_specify.units = {1};
nirs10.model_specify.wls_bglm_specify.time_res = {1};
nirs10.model_specify.wls_bglm_specify.hpf_wavelet_iter = {4};
nirs10.model_specify.wls_bglm_specify.hpf_dct_cutoff = {128};
%nirs10.model_specify.wls_bglm_specify.nirs_hpf = {hpf_none};
nirs10.model_specify.wls_bglm_specify.fwhm1 = {1.5};
%nirs10.model_specify.wls_bglm_specify.nirs_lpf = {lpf_none};