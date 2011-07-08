%Buxton-Friston VS Volterra expansion - forward model
%Load stims
%load('W:\epiNIRSb\epi127SD\dataSPM\Std\StatV0.004ord5\SPM.mat');
%Pos_Volt_Contrast_in_X = 3;
load('W:\epiNIRSb\epi133EA\dataSPM\Std\StatV0.004ord5\SPM.mat');
Pos_Volt_Contrast_in_X = 2;
Z = []; Z.b = []; Z.r = [];
%Boolean for figures
gen_fig = 0;
Slen = 900; %session length, in seconds
%loop over sessions
for s = 1:length(SPM.Sess) %session
    U = SPM.Sess(s).U(1);
    U.u = U.u(33:end); 
    ns = size(U.pst,2);
    Y.y = zeros(ns,1);
    Y.dt = U.dt;
    % priors (3 modes of hemodynamic variation)
    %--------------------------------------------------------------------------
    [pE,pC] = spm_hdm_priors(1,3);
    pE(7) = 1; %efficacy
    % model
    %--------------------------------------------------------------------------
    M.f     = 'spm_fx_hdm';
    M.g     = 'spm_gx_hdm';
    M.x     = [0 0 0 0]'; 
    M.pE    = pE;    
    M.pC    = pC;
    M.m     = 1;
    M.n     = 4;
    M.l     = 1;
    M.N     = ns;
    M.dt    = U.dt;
    M.IS    = 'spm_int';
    M.ns    = ns; 
    y   = Y.y;
    IS  = inline([M.IS '(P,M,U)'],'P','M','U');
    nr   = length(spm_vec(y))/ns;       % number of samples and responses

    % prior moments
    %--------------------------------------------------------------------------
    pE    = M.pE;
    pC    = M.pC;

    % dimension reduction of parameter space
    %--------------------------------------------------------------------------
    V     = spm_svd(pC,exp(-32));
    np    = size(V,2);                    % number of parameters (effective)
    ip    = (1:np)';

    p     = [V'*(spm_vec(M.P) - spm_vec(M.pE))];
    Ep    = spm_unvec(spm_vec(pE) + V*p(ip),pE);

    %result of forward model
    f = spm_int(Ep,M,U);
    %GLM using 1st and 2nd Volterra
    %f = X b
    X = SPM.xXn{s}.X(:,[1 Pos_Volt_Contrast_in_X]);
    K = spm_sp('Set', X);
    K.X = X;
    pX = spm_sp('x-', K);
    b = pX * f;  
    if gen_fig
        lin = linspace(1,Slen,round(Slen/M.dt));
        figure; plot(lin,f,'k'); hold on; plot(lin, X*b,'b'); hold on; 
        b(2)/b(1)
        figure; plot(lin,f,'k'); hold on; plot(lin, X*[b(1); -b(1)],'b'); hold on; 
    end
    %Store results
    Z.f{s} = f;
    Z.b = [Z.b b];
    Z.r = [Z.r b(2)/b(1)];
    Z.X{s} = X;
    Z.pX{s} = pX;
end

