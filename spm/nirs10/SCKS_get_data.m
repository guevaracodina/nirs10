function M = SCKS_get_data(M,Y)
%number of modalities
l = size(Y,2);
%data length
ns = size(Y,1);
HPF = M.HPF;
LPF = M.LPF;
M.HbRnorm = M.O.baseline_HbR;
M.HbOnorm = M.O.baseline_HbO;
M.HbTnorm = M.HbOnorm + M.HbRnorm;
    if HPF.hpf_butter_On
        ty = ButterHPF(1/M.dt,HPF.hpf_butter_freq,HPF.hpf_butter_order,ty);
    end    
    ty = private_baseline_correction(M,ty,1,r1,s1);
    Y(:,cl) = ty/M.HbRnorm; %to get percent change
end
if includeHbT
    cl = cl+1;
    ty1 = ROI{r1}{s1,cHbR};
    ty2 = ROI{r1}{s1,cHbO};
    if HPF.hpf_butter_On
        ty1 = ButterHPF(1/M.dt,HPF.hpf_butter_freq,HPF.hpf_butter_order,ty1);
        ty2 = ButterHPF(1/M.dt,HPF.hpf_butter_freq,HPF.hpf_butter_order,ty2);
    end
    ty = ty1+ty2;
    ty = private_baseline_correction(M,ty,2,r1,s1);
    Y(:,cl) =ty/M.HbTnorm; %to get percent change
end
if includeFlow
    cl = cl+1;
    ty = ROI{r1}{s1,cFlow};
    if HPF.hpf_butter_On
        ty = ButterHPF(1/M.dt,HPF.hpf_butter_freq,HPF.hpf_butter_order,ty);
    end
    ty = private_baseline_correction(M,ty,3,r1,s1);
    Y(:,cl) = ty;
end
if LPF.lpf_gauss_On
    K = get_K(1:size(Y,1),LPF.fwhm1,M.dt);
    Y = ioi_filter_HPF_LPF_WMDL(K,Y);
    for i=1:size(Y,2)
        y = Y(:,i)';
        y = ioi_filter_HPF_LPF_WMDL(K,y')';
        Y(:,i) = y;
    end
end
% if M.S.simuOn && ~M.S.simuNoise
%     Y = zeros(size(Y)); %null background
% end

M.Y.y = Y;
M.Y.dt = M.dt;

function y = private_baseline_correction(M,y,modality,r1,s1)
if ~isempty(M.O.baseline_correction)
    if size(M.O.baseline_correction{modality},1) > 1 || size(M.O.baseline_correction{modality},2) > 1
        r2 = find(r1==M.O.selected_ROIs);
        s2 = find(s1==M.O.selected_sessions);
    else
        r2 = 1;
        s2 = 1;
    end
    switch M.O.baseline_choice
        case 0
        case 1
            pctle = M.O.baseline_correction{modality}(s2,r2);
            y_offset = prctile(y,pctle);
            y = y-y_offset;
        case 2
            y = y-M.O.baseline_correction{modality}(s2,r2);
    end
end