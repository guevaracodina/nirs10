function [SCKS]  = nirs_SCKS2(SCKS,fig)
% FORMAT SCKS  = spm_SCKS2(SCKS,fig)
% NONAUGMENTED state vector
%_________________________________________________________________________
% Square-root Cubature Kalman Filters [2] & Square-root Rauch-Tang-Striebel
% Smoother (SCKF-SCKS [1]).
%==========================================================================
% This function performs joint estimation of the states, input and parameters
% of the model that is described as a stochastic continuous-discrete
% state-space in terms of nonlinear blind deconvolution. The state equations
% must have the form of ordinary differential equations, where the
% discretization is performed through local-linearization scheme [3].
% Additionally, the parameter noise covariance is estimated online via
% stochastic Robbins-Monro approximation method [4], and the measurement noise
% covariance is estimated online as well by using combination of varitional
% Bayesian (VB) approach with nonlinear filter/smoother [5].
%__________________________________________________________________________
%
% SCKS.M  - model structure (based on DEM [6] in SPM8 toolbox)
% SCKS.Y  - response variable, output or data
%    fig  - 1 = display estimates, 0 = do not display
%__________________________________________________________________________
%
% generative model:
%--------------------------------------------------------------------------
%   M(1).f  = dx/dt = f(x,v,P)    {inline function, string or m-file}
%   M(1).g  = y(t)  = g(x,v,P)    {inline function, string or m-file}
%
%   M(1).xP = state error covariance matrix
%   M(1).uP = input error variance
%   M(1).wP = parameter error covariance matrix
%
%   M(1).pE = prior expectation of p model-parameters
%   M(1).pC = prior covariances of p model-parameters
%   M(1).pP = prior process covariance of p model-parameters   % pas
%   trouvé dans le programme
%   M(1).ip = paramter indices
%   M(1).cb = constrain on parameters [lower, upper];
%
%   M(1).Q  = precision components on observation noise
%   M(1).V  = fixed precision (input noise)
%   M(1).W  = precision on state noise (approximated by annealing)
%
%   M(i).m  = number of inputs v(i + 1);
%   M(1).n  = number of states x(i);
%   M(1).l  = number of output v(i);
%
%   M(1).Qf      = form of measarument noise cov estimate:
%                  'auto'(=default),'min','mean'
%   M(1).E.nN    = number of SCKF-SCKS algorithm iterations
%   M(1).E.Itol  = tolerance value for SCKF-SCKS convergence
%   M(1).E.nD    = number of integration step between observations
%   M(1).VB.N    = number of VB algorithm iterations
%   M(1).VB.Itol = tolerance value for VB convergence
%   M(1).VB.l    = VB scaling factor;
%
%   conditional moments of model-states - q(u)
%--------------------------------------------------------------------------
%   qU.x{1}  = Conditional expectation of hidden states (backward estimate)
%   qU.x{2}  = Conditional expectation of hidden states (forward estimate)
%   qU.v{1}  = Conditional expectation of input (backward estimate)
%   qU.v{2}  = Conditional expectation of input (forward estimate)
%   qU.r{1}  = Conditional prediction of response
%   qU.z{1}  = Conditional prediction error
%   qU.S{1}  = Conditional covariance: cov(x) (states - backward estimate)
%   qU.S{2}  = Conditional covariance: cov(x) (states - forward estimate)
%   qU.C{1}  = Conditional covariance: cov(u) (input - backward estimate)
%   qU.C{2}  = Conditional covariance: cov(u) (input - forward estimate)
%
% conditional moments of model-parameters - q(p)
%--------------------------------------------------------------------------
%   qP.P    = Conditional expectation
%   qP.C    = Conditional covariance
%
%      F    = log-likelihood
%__________________________________________________________________________
% Copyright (c) Brno University of Technology (2011)
% Martin Havlicek 20-03-2011
%
% References:
% [1] Havlicek et al. (2011) Modeling neuronal responses in fMRI using
%     cubature Kalman filter. Neuroimage, In Press.
% [2] Arasaratnam, I., Haykin, S. (2009) Cubature Kalman Filters. IEEE
%     Transactions on Automatic Control 54, 1254-1269.
% [3] Jimenez, J.C. (2002) A simple algebraic expression to evaluate the
%     local linearization schemes for stochastic differential equations*
%     1. Applied Mathematics Letters 15, 775-780.
% [4] Van der Merwe, R., 2004. Sigma-point Kalman filters for probabilistic
%     inference in dynamic state-space models. Ph.D.thesis, Oregon Graduate
%     Institute of Science and Technology.
% [5] Sarkka, S., Hartikainen, J. (2011?) Extension of VB-AKF to Estimation
%     of Full Covariance and Non-Linear Systems. In Press.
% [6] Friston, K.J., et al. (2008) DEM: a variational treatment of dynamic
%     systems. Neuroimage 41, 849-885.
%__________________________________________________________________________

% check model specification
%--------------------------------------------------------------------------

nanInd=0;
SCKS.SCKS_OK = 1;
%PS = SCKS.PS; %Parameter structure
M  = SCKS.M;
%M  = ioi_SCKS_DEM_M_set(M);
M  = nirs_SCKS_DEM_M_set(M);
% get integration step dt:
dt= M(1).E.dt;    % default 1
nD = M(1).E.nD;
dt = dt/nD;        % integration step

% INITIALISATION:
% =========================================================================
% interpolate observation according to integration step
%--------------------------------------------------------------------------
y     = SCKS.Y.y;            % observations
if size(y,1)>size(y,2)     % check the dimensions
    y = y';
end
% intrepolate if dt < 1:
y    = interp1(y',[1:(1/nD):size(y,2)])';
if size(y,1)>size(y,2)     % check dimensions again
    y = y';
end
[N T]  = size(y);          % number of time points

% initial condition:
%--------------------------------------------------------------------------
x     = M(1).x;            % states
u     = M(2).v;            % input
%pE    = spm_vec(M(1).PS.pE);  % all model parameter
pE    = spm_vec(M(1).pE);  % all model parameter
ip    = M(1).ip;           % parameter indices to be estimated
theta = pE(ip);            % selected parameters

try cb  = M(1).cb;               catch, cb = []; end; % parameter constrains
try tE  = spm_vec(SCKS.pP.P{1}); catch, tE = []; end; % true paramters for display (if available)

% covariances (square-roots)
%--------------------------------------------------------------------------
sR      = cell(1,T);
[sR{:}] = deal(sparse(real(spm_sqrtm(inv(M(1).V)))*sqrt(dt)));   % observation noise variance %should be *sqrt(dt)**sqrt(N)

sQ      = sparse(real(spm_sqrtm(inv(M(1).W)))*sqrt(dt));         % hidden state noise variance
if ~isempty(M(2).v)
    sV    = sparse(real(spm_sqrtm(inv(M(2).V)))*sqrt(dt));       % input noise variance
else
    sV    = [];
end

% process error covariances (square-roots)
%--------------------------------------------------------------------------
Sx = sparse(real(spm_sqrtm(M(1).xP))*sqrt(dt));
if ~isempty(M(2).v)
    Su = sparse(real(spm_sqrtm(M(1).uP))*sqrt(dt));
else
    Su = [];
end
if ~isempty(ip)
    Sw = sparse(real(spm_sqrtm(M(1).wP(ip,ip)))*sqrt(dt));
    sW = sparse(real(spm_sqrtm(M(1).pC(ip,ip)))*sqrt(dt));  % parameter noise variance
    dv = diag(sW);
else
    Sw = [];
    sW = [];
end

% number of states, inputs and parameters:
%--------------------------------------------------------------------------
nx = size(Sx,1);             % number of states
nu = size(sV,1);             % number of states
nw = size(sW,1);             % number of paramters
no = size(sR{1},1);          % number of observations

% concatenate state vector and square-root error covariance:
%--------------------------------------------------------------------------
xc      = [x(:); u(:); theta(:)];
xx      = zeros(nx+nu+nw,T);
xx(:,1) = xc;
Sc      = cell(1,T);
[Sc{:}] = deal(sparse(nx+nu+nw,nx+nu+nw));
Sc{1}   = blkdiag(Sx,Su,Sw);

% get vector indices for components of concatenated state vector
xmask = [ones(1,nx),ones(1,nu)*2,ones(1,nw)*3,ones(1,no)*4];
xind  = find(xmask==1);
uind  = find(xmask==2);
wind  = find(xmask==3);
clear xmask;

% Precalculate cubature points:
%--------------------------------------------------------------------------
n          = nx + nu + nw;       % total state vector dimension
nPts       = 2*n;                         % number of cubature points
CubPtArray = sqrt(n)*[eye(n) -eye(n)];    % cubature points array

% augment paramter matrix by number of cubature points:
pE = pE(:,ones(1,nPts));

% setting for VB: observation noise estimation:
if ~isempty(M(1).Q)
    try, iter0  = M(1).VB.N;   catch, iter0 = 3;              end
    try, lambda = M(1).VB.l;   catch, lambda = 1-exp(-2);     end
    try, Ir     = M(1).VB.Ir;  catch, Ir = ones(size(sR{1})); end
    beta = eye(no);
    alpha = 1;
    %      alpha =1/1e-4;%test
    [sR{:}] = deal(sqrt(beta/alpha));
    
    iter  = iter0;
    MSE0  = zeros(no,1);
    RR0   = zeros(no,T-1);
    VBrun = [];
    dl    =  (1-lambda);
else
    iter0 = 1;
    iter  = iter0;
    RR    = repmat(diag(sR{1}),1,T-1);
    VBrun = [];
end

% prepare matrix template for integration by Local linearization scheme:
%--------------------------------------------------------------------------
EXPm = repmat({[ones(nx),2*ones(nx,1);zeros(1,nx+1)]},1,nPts);
EXPm = sparse(blkdiag(EXPm{:}));
xt   = repmat([zeros(1,nx) 1],1,nPts)';

OnesNpts = ones(1,nPts);
xPred    = zeros(n,nPts);
yPred    = zeros(N,nPts);
Xf       = cell(1,T-1);
[Xf{:}]  = deal(sparse(nx+nu+nw,nPts));
x1f      = zeros(nx+nu+nw,T-1);

% Initialize display:
%--------------------------------------------------------------------------
try, M(1).nograph; catch, M(1).nograph = 0; end
if ~M(1).nograph
    f1 = spm_figure('Create','Graphics','SCKF-SCKS estimates');
    movegui(f1,'northeast'); drawnow;
end

% ==================================================================
% Iteration scheme:
% ==================================================================
% get maximum number of iterations and tolerance:
try,  ItolVB = M(1).VB.Itol;  catch,  ItolVB = 1e-4;      end
try,  Itol   = M(1).E.Itol;   catch,  Itol   = 1e-3;      end
try,  RUN    = M(1).E.nN;     catch,  RUN    = 60;        end
try,  ap     = M(1).E.RM;     catch,  ap     = [1e3 1e6]; end   % Robins-Monro approximation paramters

MLdiff0  = 0;
MLdiff0  =1e-4; %test
mloglik0 = 0;
ML       = [];
VBrun    = RUN;
EXEC     = 0;
t0       = tic;
% =========================================================================
% Iteration loop (until convergence)
% =========================================================================

for run = 1:RUN
    t1 = tic;
    mloglik = -log(2.*pi).*(T/dt);
    % ==================================================================
    %   Forward pass:
    % ==================================================================
    for t = 2:T
        %ww=0.9995; %default No big changes
        sQ = diag(diag((1/sqrt(M(1).A)-1)*Sc{t-1}(xind,xind)));
        if ~isempty(sW), sW = (1/sqrt(M(1).Ap)-1)*sW; end
        Xi =  xc(:,OnesNpts) + Sc{t-1}*CubPtArray;
        %------------------------------------------------------------------
        % PREDICTION STEP:
        %------------------------------------------------------------------
        xPred(uind,:) = Xi(uind,:);%
        xPred(wind,:) = Xi(wind,:);%
        
        % parameter constrain:
        if ~isempty(cb) && ~isempty(ip)
            xPred(wind,:) = min(cb(:,2*OnesNpts),xPred(wind,:)); % upper constrain
            xPred(wind,:) = max(cb(:,1*OnesNpts),xPred(wind,:)); % lower constrain
        end
        
        pE(ip,:)      = xPred(wind,:);
        %PS.pE = pE;
        % propagation of cubature points through nonlinear function:
        %------------------------------------------------------------------
        %f             = M(1).f(Xi(xind,:),xPred(uind,:),PS);
        f             = M(1).f(Xi(xind,:),xPred(uind,:),pE,M(1));
        % integration by local-linearization scheme:
        %------------------------------------------------------------------
        %dfdx          = spm_diff_all(M(1).f,Xi(xind,:),xPred(uind,:),PS,1);
        dfdx          = spm_diff_all(M(1).f,Xi(xind,:),xPred(uind,:),pE,M(1),1);
        dx            = expmall(dfdx,f,dt,EXPm)*xt;
        xPred(xind,:) = Xi(xind,:) + reshape(dx(~xt),nx,nPts);
        % mean prediction:
        %------------------------------------------------------------------
        x1            = sum(xPred,2)/nPts;
        X             = (xPred-x1(:,OnesNpts))/sqrt(nPts);
        Xf{t-1}       = X;  % store for the backwards run (then no need to propaget through the nonlinear fcn)
        x1f(:,t-1)    = x1; % store for the backwards run
        [~,S]         = qr([X blkdiag(sQ,sV,sW)]',0);
        S             = S';
        
        % in the case of VB observation noise estimation:
        %---------------------------------------------------------
        if ~isempty(M(1).Q) && iter~=1
            if isfield(M(1),'alpha') && M(1).alpha==1  %  % this is former version from SPM_SCKS2
                alpha     = lambda.*alpha + 1;
            elseif  isfield(M(1),'alpha') && M(1).alpha==2 % this is what I think is ok
                alpha     = lambda.*alpha + 1/N;
            else    %Default :this is modification from Havlicek and his comment in SPM_SCKS2new
                alpha     = lambda.*alpha + dt/N; % this my modification for this particular case
                % otherwise it should be + 1
            end
            beta      = lambda.*beta;
            beta0     = beta;
        end
        
        Xi            = x1(:,OnesNpts) + S*CubPtArray;
        X             = (Xi-x1(:,OnesNpts))/sqrt(nPts);
        pE(ip,:)      = Xi(wind,:);
        %PS.pE = pE;
        % propagate cubature points through observation function:
        %yPred = M(1).g(Xi(xind,:),Xi(uind,:),PS);
        yPred = M(1).g(Xi(xind,:),Xi(uind,:),pE,M(1));
        y1    = sum(yPred,2)/nPts;
        Y     = (yPred-y1(:,OnesNpts))/sqrt(nPts);
        
        Pxy    = X*Y';                % cross covariance
        resid  = y(:,t) - y1;         % innovations
        %------------------------------------------------------------------
        % UPDATE STEP:
        %------------------------------------------------------------------
        % VB estimation of sR (iteratively)
        for it = 1:iter
            % VB part - update of square-root measurement noise cov
            if ~isempty(M(1).Q) && iter~=1
                Rtype = 'diag';
                switch(Rtype)
                    case('full')
                        sR{t-1} = spm_sqrtm(beta./alpha);
                    case('diag')
                        sR{t-1} = diag(sqrt(diag(beta./alpha)));
                    case('mean')
                        sR{t-1} = mean(sqrt(diag(beta./alpha)))*eye(no);
                end
                sR{t-1}= sR{t-1}.*Ir;
            end
            
            
            [~,Sy] = qr([Y sR{t-1}]',0);
            Sy     = Sy';
            
            
            K   = (Pxy/Sy')/Sy;        % Kalman gain
            if isnan(K)
                nanInd=nanInd+1;
                if nanInd==100
                    SCKS.SCKS_OK = 0;
                    return
                end
            end
            
            % state (input,parameter) estimates:
            xc = x1 + K*(resid);
            
            % estimate of process error covarinace:
            [~,S]   = qr([(X - K*Y) K*sR{t-1}]',0);
            S       = S';
            
            % check parameter constrain:
            if ~isempty(cb) && ~isempty(ip)
                xc(wind) = min(cb(:,2),xc(wind)); % upper constrain
                xc(wind) = max(cb(:,1),xc(wind)); % lower constrain
            end
            
            % VB part:
            if ~isempty(M(1).Q) && iter~=1
                Xi       = xc(:,OnesNpts) + S*CubPtArray;
                pE(ip,:) = Xi(wind,:);
                %PS.pE = pE;
                %yPred(:) = M(1).g(Xi(xind,:),Xi(uind,:),PS); % no additive noise here!
                yPred(:) = M(1).g(Xi(xind,:),Xi(uind,:),pE,M(1)); % no additive noise here!
                D        = (y(:,t*OnesNpts)-yPred)/sqrt(nPts);
                beta     = beta0 + D*D';
            end
            
        end
        
        Sc{t}   = S;
        xx(:,t) = xc;
        
        % Maximum log-Likelihood
        %------------------------------------------------------------------
        mloglik = mloglik - log(det(Sy*Sy')) - resid'/(Sy*Sy')*resid;
        
        % Robins-Monro stochastic approximation for of parameters noise cov
        %------------------------------------------------------------------
        if ~isempty(ip)
            subKG = K(wind,:);
            dv    = sqrt((1-1/ap(1))*(dv.^2) + 1/ap(1)*diag(subKG*(subKG*resid*resid')'));
            sW    = diag(dv);
            ap(1) = min(ap(1)+dt,ap(2));
        end
        if ~isempty(M(1).Q) && iter~=1
            RR(:,t-1) = diag(sR{t-1});
        end
        
        
    end
    xxf = xx;
    Sf  = Sc;
    if ~isempty(M(1).Q) && iter~=1
        lambda = min(lambda + dl*(1-exp(-run)), 1);
    end
    
    
    %----------------------------------------------------------------------
    % END of forward pass
    % ---------------------------------------------------------------------
    
    % ==================================================================
    %   Backward pass:
    % ==================================================================
    for t = T-1:-1:1
        
        % Square-root Cubature Rauch-Tung-Striebel smoother
        %------------------------------------------------------------------
        % evaluate cubature points:
        
        Xi =  xx(:,t*OnesNpts) + Sc{t}*CubPtArray;
        
        x01      = xx(:,t);
        X01      = (Xi -  x01(:,OnesNpts))/sqrt(nPts);
        [~,S]    = qr([Xf{t} blkdiag(sQ,sV,sW)]',0);
        S        = S';
        
        Pxy      = X01*Xf{t}';      % cross covariance
        K        = (Pxy/S')/S;  % Kalman gain
        
        % smoothed estimate of the states (input,parameters)
        % and process error covariance:
        xx(:,t)  = xx(:,t) + K*(xx(:,t+1) - x1f(:,t));
        [~,S]    = qr([X01 - K*Xf{t}, K*blkdiag(sQ,sV,sW),K*Sc{t+1}]',0);
        S        = S';
        Sc{t}    = S;
        
        % check parameter constrain:
        if ~isempty(cb) && ~isempty(ip)
            xx(wind,t) = min(cb(:,2),xx(wind,t)); % upper constrain
            xx(wind,t) = max(cb(:,1),xx(wind,t));
        end
        
    end
    xxb = xx;
    Sb  = Sc;
    %----------------------------------------------------------------------
    % END of backward pass
    %----------------------------------------------------------------------
    
    str{1} = sprintf('SCKS: %i (1:%i)',run,iter);
    
    % iteration condition for measurement noise estimate:
    % iterate until stabilization of sR estimate
    %------------------------------------------------------------------
    if ~isempty(M(1).Q) && iter0(run)~=1
        MSE     = mean((RR -(RR0)).^2,2);
        RR0     = RR;
        MSEdiff = abs(MSE - MSE0);
        MSE0    = MSE;
        iter0(run+1) = iter0(run);
        if MSEdiff<ItolVB  % (till it gets stable)
            switch(lower(M(1).Qf))
                case('all')
                    % take all
                case('mean')
                    sR      = cell(1,T);
                    [sR{:}] = deal(diag(mean(RR,2)));
                case('mean-all')
                    sR      = cell(1,T);
                    [sR{:}] = deal(eye(no)*(mean(mean(RR,2))));
                case('min')
                    RRs     = sort(RR,2,'descend');
                    sR      = cell(1,T);
                    [sR{:}] = deal(diag(mean(RRs(:,round(T*0.90):end),2)));
                case('auto')
                    dlim    = min(RR,[],2);
                    ulim    = max(RR,[],2);
                    if all(ulim./dlim<4)
                        % take all
                    else
                        RRs     = sort(RR,2,'descend');
                        sR      = cell(1,T);
                        [sR{:}] = deal(diag(mean(RRs(:,round(T*0.90):end),2)));
                    end
            end
            %
            iter0(run+1) = 1;
            iter     = iter0(run+1);
            mloglik0 = 0;
            VBrun    = run;
        end
    else
        iter0(run+1) = 1;
    end
    %----------------------------------------------------------------------
    
    % log-likelihood difference:
    MLdiff(run) = (mloglik-mloglik0);
    ML(run)     = mloglik;
    
    timed  = toc(t1);
    str{2} = sprintf('F:%.4e',ML(end));
    str{3} = sprintf('dF:%.4e',MLdiff(end));
    str{4} = sprintf('(%.2e sec)',timed);
    fprintf('%-16s%-16s%-16s%-16s\n',str{:})
    
    
    %First old stopping condition from Havlicek:
% %     if RUN>1 && (~isempty(ip) || ~isempty(M(1).Q))
% %         if run==2
% %             MLdiff0 = MLdiff(run);
% %         elseif run>2
% %             if MLdiff0<MLdiff(run),
% %                 MLdiff0 = MLdiff(run);
% %             end
% %         end
% %         if (((MLdiff(run)/MLdiff0)<Itol || run==RUN) || (isempty(ip)&&MLdiff(run)<Itol)) && iter0 == 1,
% %       
    %Second old stopping condition from Havlicek
% %     if (MLdiff(run)<0 && run>1 && ((isempty(M(1).Q) && sum(iter0(run:end))==2) || ...
% %                                    (~isempty(M(1).Q) && sum(iter0(run-1:end))==3)))
 
                               
    if (MLdiff(run)<0 && run>1 && abs(MLdiff(run)) < 1e-1) %&& ((isempty(M(1).Q) && sum(iter0(run:end))==2) || ...
        %   (~isempty(M(1).Q) && sum(iter0(run-1:end))==3)))
        EXEC  = 1;
    else
        
        XXf = xxf;
        XXb = xxb;
        SSf = Sf;
        SSb = Sb;
    end    
        % plot estimates:
        %----------------------------------------------------------------------
        if ~M(1).nograph
            doplotting(SCKS,M,XXf,XXb,SSf,SSb,ML,f1,T,wind,ip,run,RR,VBrun);
        end
    %end
    %----------------------------------------------------------------------
    % stopping condition:
    if RUN>1 %&& (~isempty(ip) || ~isempty(M(1).Q))
        if run==2
            MLdiff0 = MLdiff(run);
        elseif run>2
            if MLdiff0<MLdiff(run),
                MLdiff0 = MLdiff(run);
            end
        end
%         if ( (abs(MLdiff(run)/MLdiff0)<Itol || run==RUN) || ...
%                 (isempty(ip) && MLdiff(run)<Itol) || EXEC ) && run>1,
%             
        if EXEC && run> 1
            timed = toc(t0);
            if (run~=RUN || EXEC)
                fprintf('Converged (in %2.2e sec)\n',timed);
            else
                fprintf('Reached the maximum of iterations (in %2.2e sec)\n',timed);
            end
            if ~isempty(ip)
                fprintf('Estimated parameters: %4.2f\n',mean(full(XXb(wind,:)),2));
            end
            pE(ip,1) = mean(XXb(wind,:),2);
            %PS.pE = pE;
            %yy       = M(1).g(XXb(xind,:),XXb(uind,:),PS);
            yy       = M(1).g(XXb(xind,:),XXb(uind,:),pE,M(1));
            res      = y - yy;
            
            try SCKS = rmfield(SCKS,'qU'); end
            try SCKS = rmfield(SCKS,'qP'); end
            try SCKS = rmfield(SCKS,'qH'); end
            % save results into structure:
            SCKS.qU.x{2} = XXf(xind,1:nD:end);
            SCKS.qU.x{1} = XXb(xind,1:nD:end);
            SCKS.qU.v{2} = XXf(uind,1:nD:end);
            SCKS.qU.v{1} = XXb(uind,1:nD:end);
            SCKS.qU.r{1} = yy;
            SCKS.qU.z{1} = res;
            if ~isempty(ip)
                SCKS.qP.P{1} = XXb(wind,1:nD:end);
                SCKS.qP.P{2} = XXf(wind,1:nD:end);
            end
            SCKS.F       = ML;
            
            for i = 1:nD:T
                j = 1 + (i - 1)/nD;
                SCKS.qU.S{1}(:,j) = diag(SSb{i}(xind,xind));
                SCKS.qU.S{2}(:,j) = diag(SSf{i}(xind,xind));
                SCKS.qU.C{1}(:,j) = diag(SSb{i}(uind,uind));
                SCKS.qU.C{2}(:,j) = diag(SSf{i}(uind,uind));
                if ~isempty(ip)
                    SCKS.qP.C{1}(:,j) = diag(SSb{i}(wind,wind));
                    SCKS.qP.C{2}(:,j) = diag(SSf{i}(wind,wind));
                end
            end
            SCKS.run=run;
            return
        end
        
        mloglik0 = mloglik;
        xc       = [xx([xind,uind],1); mean(xx(wind,:),2)];
        xx(:,1)  = xc;
        Sc{1}    = Sc{end};
        
    else
        pE(ip,1) = mean(XXb(wind,:),2);
        %PS.pE = pE;
        %yy       = M(1).g(XXb(xind,:),XXb(uind,:),PS);
        yy       = M(1).g(XXb(xind,:),XXb(uind,:),pE);
        res      = y - yy;
        
        try SCKS = rmfield(SCKS,'qU'); end
        try SCKS = rmfield(SCKS,'qP'); end
        try SCKS = rmfield(SCKS,'qH'); end
        % save results into structure:
        SCKS.qU.x{2} = XXf(xind,1:nD:end);
        SCKS.qU.x{1} = XXb(xind,1:nD:end);
        SCKS.qU.v{2} = XXf(uind,1:nD:end);
        SCKS.qU.v{1} = XXb(uind,1:nD:end);
        SCKS.qU.r{1} = yy;
        SCKS.qU.z{1} = res;
        SCKS.qP.P{1} = XXb(wind,1:nD:end);
        SCKS.qP.P{2} = XXf(wind,1:nD:end);
        SCKS.F       = ML;
        
        for i = 1:nD:T
            j = 1 + (i - 1)/nD;
            SCKS.qU.S{1}(:,j) = diag(SSb{i}(xind,xind));
            SCKS.qU.S{2}(:,j) = diag(SSf{i}(xind,xind));
            SCKS.qU.C{1}(:,j) = diag(SSb{i}(uind,uind));
            SCKS.qU.C{2}(:,j) = diag(SSf{i}(uind,uind));
            if ~isempty(ip)
                SCKS.qP.C{1}(:,j) = diag(SSb{i}(wind,wind));
                SCKS.qP.C{2}(:,j) = diag(SSf{i}(wind,wind));
            end
            
        end
        SCKS.run=run;
        return
    end
end

%==========================================================================
%==========================================================================


%--------------------------------------------------------------------------
% Plot estimates at each iteration:
%--------------------------------------------------------------------------
function doplotting(SCKS,M,xxf,xxb,Sf,Sb,ML,f1,T,wind,ip,run,RR,VBrun)
if SCKS.DO.save_figures || SCKS.DO.generate_figures
    figure(f1);
    set(f1,'Renderer','painter');
    clf(f1);
    for p=1:2
        subplot(3,3,[1:3]+3*(p-1)),
        hax = gca;
        si    = spm_invNcdf(1 - 0.05);
        s     = [];
        if p == 1,
            xxfig = xxf;
            Sfig  = Sf;
            tit   = 'SCKF - forward pass';
        else
            xxfig = xxb;
            Sfig  = Sb;
            tit   = 'SCKS - backward pass';
        end
        for i = 1:T
            s = [s abs(diag(Sfig{i}))];
        end
        
        % conditional covariances
        %------------------------------------------------------------------
        j           = [1:size(xxfig(:,:),1)];
        ss          = si*s(j,:);
        [ill indss] = sort(full(mean(ss,2)),'descend');
        
        pf = plot(1:T,xxfig,'linewidth',1.5);
        set(hax,'xlim',[1,T],'nextplot','add')
        box(hax,'on')
        for ic = 1:size(xxfig,1)
            col0 = get(pf(indss(ic)),'color');
            col  = (ones(1,3)-col0)*0.65 + col0;
            fill([(1:T) fliplr(1:T)],[(xxfig(indss(ic),:) + ss(indss(ic),:)) fliplr((xxfig(indss(ic),:) - ss(indss(ic),:)))],...
                'r',...
                'FaceColor',col,...
                'EdgeColor',col);
            hold on;
            COL{ic} = col0;
        end
        for ic = 1:size(xxfig,1)
            plot(xxfig(indss(ic),:),'color',COL{ic},'linewidth',0.75);
        end
        title(tit);
        grid(hax,'on');
        axis(hax,'tight');
    end
    
    subplot(3,3,7)
    h = plot([1:length(ML)],ML);
    if ~isempty(M(1).Q)
        AYlim = get(gca,'Ylim');
        if VBrun>=run
            bkg = ones(1,run)*max(AYlim);
        else
            bkg = [ones(1,VBrun)*max(AYlim),ones(1,abs(VBrun-run))*min(AYlim)];
        end
        a = area(bkg,min(AYlim));
        axis([1 run+1 AYlim(1) AYlim(2)]); hold on;
        set(a(1),'FaceColor',ones(1,3)*0.6,'EdgeColor',ones(1,3)*0.6)
        h = plot([1:length(ML)],ML);
        axis([1 run+1 AYlim(1) AYlim(2)]); hold on;
        if VBrun>=run
            text(run/2+1,mean(AYlim),'VB-SCKS','HorizontalAlignment','center','VerticalAlignment','top');
        else
            text(VBrun/2+1,mean(AYlim),'VB-SCKS','HorizontalAlignment','center','VerticalAlignment','top');
            text(VBrun+(abs(VBrun-run)/2)+1,mean(AYlim),'SCKS','HorizontalAlignment','center','VerticalAlignment','top');
        end
        set(gca,'Layer','top')
        hold off
    end
    set(h,'color','k','Marker','o','Markersize',4,'MarkerFaceColor','k','linewidth',1);
    title('Log-Likelihood');
    % disp(mean(xxb(wind,:),2)')
    if ~isempty(ip)
        subplot(3,3,8)
        b1 = bar(ip,mean(xxb(wind,:),2)','FaceColor','k');
        set(b1,'BarWidth',0.5); hold on;
        title('Parameters');
    end
    
    if ~isempty(M(1).Q)
        subplot(3,3,9)
        if VBrun>=run
            plot(RR'); hold on;
        else
            plot(RR'); hold on;
            switch(lower(M(1).Qf))
                case('all')
                    sR  = RR';
                case('mean')
                    sR  = repmat(mean(RR,2),1,T)';
                case('min')
                    RRs = sort(RR,2,'descend');
                    sR  = repmat(mean(RRs(:,round(T*0.90):end),2),1,T)';
                case('auto')
                    dlim    = min(RR,[],2);
                    ulim    = max(RR,[],2);
                    if all(ulim./dlim<4)
                        sR  = RR';
                    else
                        RRs  = sort(RR,2,'descend');
                        sR  = repmat(mean(RRs(:,round(T*0.90):end),2),1,T)';
                    end
            end
            plot([1:T-1],sR,'r','linewidth',2);
        end
        axis(gca,'tight');
        title('Std(R)');
        hold off;
    end
    drawnow;
end
if SCKS.DO.save_figures
    nirs_SCKS_print_figures(SCKS,f1,['EstIt' gen_num_str(run,2)]);
end
% commentaire
% la nouvelle version pratiquement juste une différence lorsuqe qU n'est
% pas vide