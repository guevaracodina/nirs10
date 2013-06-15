function CINE = prepare_cine_core_call(Z,W,NIRS,CINE)
NC = NIRS.Cf.H.C.N;
fs = NIRS.Cf.dev.fs;
%Information on how to use onsets
OnsetChoice = Z.onsetInfo.OnsetChoice;
onset_delay = Z.onsetInfo.onset_delay;
onset_duration = Z.onsetInfo.onset_duration;
time_resolution = Z.onsetInfo.time_resolution;
baseline_offset = Z.onsetInfo.baseline_offset;
baseline_duration = Z.onsetInfo.baseline_duration;
if length(onset_delay) == 1
    onset_delays = onset_delay:time_resolution:(onset_delay+onset_duration);
else
    onset_delays = onset_delay;
end
LOD = length(onset_delays);

%Filters
%HPF - Butterworth infinite impulse response filter
hpf_butter = Z.CineFilters.hpf_butter;
if isfield(hpf_butter,'hpf_butter_On')
    hpf_butter_freq = hpf_butter.hpf_butter_On.hpf_butter_freq;
    hpf_butter_order = hpf_butter.hpf_butter_On.hpf_butter_order;
    HPFbutter = 1;
else
    if isfield(hpf_butter,'remove_linear')
        HPFbutter = 2;
        hpf_butter_freq = 0; %not used
        hpf_butter_order = 3; %not used
    else
        HPFbutter = 0; %no high pass filter
        hpf_butter_freq = 0; %not used
        hpf_butter_order = 3; %not used
    end
end

%LPF - filter from NIRS_SPM
nirs_lpf = Z.CineFilters.nirs_lpf;
if isfield(nirs_lpf,'lpf_gauss')
    FWHM = nirs_lpf.lpf_gauss.fwhm1;
    LPF = 'gaussian';
else
    if isfield(nirs_lpf,'lpf_hrf')
        LPF = 'hrf';
    else
        if isfield(nirs_lpf,'lpf_none')
            LPF = 'none';
        else
            disp('Unrecognized low pass filter');
        end
    end
end

%Build K structure for low pass filter
K.HParam.type = 'none';
if isempty(strfind(LPF, 'hrf')) == 0 % hrf smoothing
    K.LParam.type = 'hrf';
elseif isempty(strfind(LPF, 'gaussian')) == 0 % Gaussian smoothing
    K.LParam.FWHM = FWHM;
    K.LParam.type = 'Gaussian';
else
    K.LParam.type = 'none';
end

%Z.Avg = 0; %This sets whether we are looking at averaged rather than GLM data
%Objective is to specify W.beta, W.res, W.var and W.corr_beta
%then pass that to contrast_core
try
    %Sessions are processed individually
    Nsess = length(NIRS.Dt.fir.pp(end).p);
    for f1=1:Nsess
        if any(Z.sessions == 0) || any(f1 == Z.sessions)
            %Load NIRS data
            Y = fopen_NIR(NIRS.Dt.fir.pp(end).p{f1},NC);
            Y = Y';
            %apply filters
            %HPF
            switch HPFbutter
                case 1
                    KY = ButterHPF(fs,hpf_butter_freq,hpf_butter_order,Y);
                case 2
                    nS = size(Y,1);
                    mX = linspace(0,round(nS/fs),nS);
                    mX = [mX' ones(nS,1)];
                    pmX = pinv(mX);
                    KY = Y - mX * (pmX * Y);
            end
            %LPF
            svec = 1:size(Y,1);
            K = struct('HParam',K.HParam,...
                'row',svec ,'RT',1/fs,'LParam',K.LParam);
            K = spm_filter_HPF_LPF_WMDL(K);
            K.row = svec; %1:length(svec);
            KY = spm_filter_HPF_LPF_WMDL(K,KY);
            
            %Get onset information
            U = NIRS.Dt.fir.Sess(f1).U;
            
            %loop over onset types
            for u1 = OnsetChoice
                %loop over onsets
                ons = U(OnsetChoice).ons;
                nons = length(ons);
                for k0 = 1:nons
                    beta = zeros(LOD,NC);%stimulus types by channels
                    covbeta = zeros(LOD,NC);
                    %Baseline
                    bo = round(fs*(ons(k0)-baseline_offset));
                    bidx = round((bo-fs*baseline_duration)):bo;
                    if bidx(1) > 1 %else skip onset
                        B0 = std(KY(bidx,:),0,1).^2; %Baseline variance
                        M0 = mean(KY(bidx,:),1); %Baseline mean
                        for t0=1:LOD
                            idx1 = round(fs*(ons(k0)+onset_delays(t0)));
                            idx = idx1:round(idx1+fs*onset_duration-1);
                            if idx(end) <= size(KY,1);
                                beta(t0,:) = mean(KY(idx,:),1)-M0;
                                covbeta(t0,:) = B0;
                                %t = SPM.xX.beta./SPM.xX.covbeta.^(1/2);
                            else
                                LODmax = t0-1;
                                beta = beta(1:LODmax,:);
                                covbeta = covbeta(1:LODmax,:);
                                disp(['Only ' int2str(LODmax) ' time points can be included for onset ' int2str(k0) ' of onset type ' int2str(u0)]);
                            end
                        end
                    end
                    W.beta = beta;
                    W.covbeta = covbeta;
                    %W.t = t;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    CINE = cine_core(Z,W,CINE,f1,u1,k0);                   
                end
            end
        end
    end
catch exception
    disp(exception.identifier)
    disp(exception.stack(1));
    disp('Problem in prepare_cine_core_call');
end