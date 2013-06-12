function Y = nirs_SCKS_filter(Y,LPF,HPF)
if LPF.lpf_gauss_On
    K = get_K(1:size(Y.y,1),LPF.fwhm1,Y.dt);
    Y.y = nirs_filter_HPF_LPF_WMDL(K,Y.y);   
end
Y.X0 = [];
switch HPF.hpf_filter
    case 1
        k       = size(Y.y,1);
        n       = fix(2*(k*Y.dt)*HPF.cosine_freq + 1);
        X0      = spm_dctmtx(k,n);
        Y.X0 = X0(:,2:end);
        if HPF.band_pass_filter
            n2      = fix(2*(k*Y.dt)*HPF.cosine_freq_band(1) + 1);
            X2      = spm_dctmtx(k,n2);
            n3      = fix(2*(k*Y.dt)*HPF.cosine_freq_band(2) + 1);
            X3      = spm_dctmtx(k,n3);
            X = X3(:,(size(X2,2)+1):end);
            Y.X0 = [Y.X0 X];
        end
    case 2
        Y.y = ButterHPF(1/Y.dt,HPF.hpf_butter_freq,HPF.hpf_butter_order,Y.y);
end