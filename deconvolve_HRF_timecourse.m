function [estiHRF dm] = deconvolve_HRF_timecourse(Y, onsets, fs, time_before_spk, time_after_spk, EvtofInterest, EvtofInterest_name, chromophore, decov)

Nc = size(Y,1);
switch chromophore
    case 1 %HbO
        YC = Y(1:(Nc/3),:);
    case 2 %HbR
        YC = Y((Nc/3+1):(2*Nc/3),:);
    case 3 %HbT
        YC = Y((2*Nc/3+1):Nc,:);
end

ons_sz = [];
ons_spk = [];
for i0 = 1 : length(onsets)
    name_ons = onsets(i0).name{1};
    if ~isempty(EvtofInterest_name) && strcmpi(name_ons, EvtofInterest_name)
        spk_ons = onsets(i0).ons;
        EvtofInterest_name = name_ons;
    elseif  strcmpi(name_ons, 'sz')
        ons_sz = [ons_sz i0];
    elseif (~isempty(strfind(name_ons, 'spk'))) || (~isempty(strfind(name_ons, 'Spk'))) || ~isempty(strfind(name_ons, 'SPK'))
        if i0 ~= EvtofInterest
            ons_spk = [ons_spk i0];
        end
    end
end
if EvtofInterest ~= 0
    spk_ons = onsets(EvtofInterest).ons;
    EvtofInterest_name = onsets(EvtofInterest).name;
end
disp(['Total number of ' int2str(length(spk_ons)) ' Event of interest (' EvtofInterest_name ') in the session']);
disp('-----------------------------------------------------------------');

PointsBefore = round(abs(time_before_spk)*fs);
PointsAfter = round(abs(time_after_spk)*fs);
total_num = PointsBefore+PointsAfter+1;
total_length = size(YC, 2);
poly_order = decov.poly_order;
X = zeros(total_length - total_num + 1,total_num + 1 + poly_order); % size of the design matrix

%Define X
if ~isempty(spk_ons)
    % Stimulus term
    spk_ons_idx = round(spk_ons*fs); %+ 33; % Seconds -> Index in time sequence, 32-bin offset
    f0 = zeros(1,size(YC,2));
    f0(spk_ons_idx) = 1; %Stimulus sequence
    idx = 1 : (total_length-total_num+1);
    for i0 = 0 : (total_num-1)
        X(:, total_num-i0) = f0(idx+i0)';
    end
    
    % Constant
    X(:,total_num + 1) = 1;
    
    % Polynormial term
    poly = (PointsAfter + 1) : (total_length - PointsBefore);
    for i0 = 1 : poly_order
        X(:, total_num + 1 + i0) = (poly .^ i0)';
    end
end

%Define Z
Z = YC(:, (PointsAfter + 1) : (total_length - PointsBefore));

%Estabish GLM
pXpX = pinv(X'*X);
%pXpX = spm_sp('pinvxpx',X);
beta = pXpX * X' * Z'; %OLS
%Residual
res = Z' - X*pXpX*X'*Z'; %e=(I-X*pinv(X'X)*X')Z'
ResSS = sum(res.^2); %Residual sum of squares = Z*R*Z' = Z*(I-H)*Z'

if decov.kl == 1 %if ols return
    estiHRF = beta(1:total_num,:)';
    dm.X = X;
    dm.Z = Z;
    dm.beta = beta;
    dm.ResSS = ResSS;
elseif decov.kl == 2 %elseif ridge regression
    k_mtd = decov.k_val;
    rv_file = decov.rv_file;
    if k_mtd == 1 %derive k value with Hoerl's single iterative method
        %Calculate k
        XTX = X'*X;
        p = rank(XTX);
        pre_coloring = 1;
        if pre_coloring
            %calculate standard deviation of beta
            try
                load(rv_file);
            catch
                HParam.type = 'none';
                LParam.type = 'hrf';
                row = 1 : size(Z,2);
                K = struct( 'HParam', HParam,...
                    'RT', 1/fs ,...
                    'row', row,...
                    'LParam', LParam);
                K = spm_filter_HPF_LPF_WMDL(K);
                S = K.KL; %Tempoal autocorrelation kernel
                pKX = pXpX * X'; %pKX = pinv(X'X)*X'
                %calculate effective DOF
                Var1 = S * S'; %Var1 = V;
                Var2 = Var1 - X * (pKX * Var1); %R = I - X*pinv(X'X)X'; Var2 = RV;
                trRV = sum(diag(Var2));
            end
            ervar = ResSS/trRV;%Usual estimator of variance of error of betas in a least-square scheme
        end
        BSS = beta' * beta;
        k = (p*ervar)./diag(BSS)';
    end
    %Establish ridge regression, channel-wise
    beta_r = zeros(size(beta));
    ResSS_r = zeros(size(ResSS));
    for k0 = 1 : length(k)
        ki = k(k0);
        XpX = XTX + ki.*eye(size(XTX,1));
        pXpX_r = pinv(XpX);
        beta_r(:,k0) = pXpX_r * X' *Z(k0,:)';
        bi = beta_r(:,k0) - beta(:,k0);
        ResSS_r(k0) = ResSS(k0) + bi'*XTX*bi;
    end
    estiHRF = beta_r(1:total_num,:)';
    dm.X = X;
    dm.Z = Z;
    dm.beta = beta_r;
    dm.ResSS = ResSS_r;    
end
