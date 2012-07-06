function M = nirs_get_hdm_data(M,Y,c1,s1)
Y0 = Y;
Y = Y.y;
%filtering and baseline correction
HPF = M.O.HPF;
LPF = M.O.LPF;
if HPF.hpf_butter_On
    Y = ButterHPF(1/M.dt,HPF.hpf_butter_freq,HPF.hpf_butter_order,Y);
end
if M.IC.include_HbR
    %Y = private_baseline_correction(M,Y,1,c1,s1);
    Y(:,1) = Y(:,1)/M.O.baseline_HbR; %to get percent change
end
if M.IC.include_HbT || M.IC.include_HbO
    %Y = private_baseline_correction(M,Y,2,c1,s1);
    Y(:,end) =Y(:,end)/(M.O.baseline_HbO+M.O.baseline_HbR); %to get percent change
end

if LPF.lpf_gauss_On
    K = get_K(1:size(Y,1),LPF.fwhm1,M.dt);
    for i=1:size(Y,2)
        y = Y(:,i)';
        y = nirs_filter_HPF_LPF_WMDL(K,y')';
        Y(:,i) = y;
    end
end
M.Y = Y0;
M.Y.y = Y;

function y = private_baseline_correction(M,y,modality,c1,s1)
if ~isempty(M.O.baseline_correction)
    switch M.O.baseline_choice
        case 0
        case 1
            pctle = M.O.baseline_correction{modality}(s1,c1);
            y_offset = prctile(y,pctle);
            y = y-y_offset;
        case 2
            y = y-M.O.baseline_correction{modality}(s1,c1);
    end
end