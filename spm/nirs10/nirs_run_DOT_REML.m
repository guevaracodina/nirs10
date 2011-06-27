function [lambda,Beta,Stats]=nirs_run_DOT_REML(Y,X,Beta_prior,Qn,Qp,maxIter)
% Restricted Maximum Likelihood code.  This program was modified from the
% SPM function spm_reml to be optimized for the dimensions of optical data.
% For reference see
%       Friston KJ. Statistical parametric mapping : the analysis of functional
%       brain images. London: Academic; 2007.
%

if(~exist('maxIter'))
    maxIter=35;  %Max # of iterations of REML code
end

if(size(X,1)<size(X,2))
   % If X is not full rank, then we can make this MUCH faster
   %  noting that the hyper-parameters of the reduced problem are the same as
   %  the original forward model.  Note- this was not done in the paper, but is
   %  a worthwhile extension in future work.  This trick was discovered after
   %  we finished the paper in prep of this demo.

    %  Y = U*S*V'*Beta
    %  Y = U*(S*V'*Beta) --> Y=U*S*Beta2;
    %  cov(Beta2) = S*V'*Q*V*S';
    [U,S,V]=svd(full(X),'econ');

    for idx=1:length(Qp)
        Qp2{idx}=S*V'*Qp{idx}*V*S';
    end
    Beta_prior=S*V'*Beta_prior;
    [lambda,Beta,Stats]=nirs_run_DOT_REML(Y,U,Beta_prior,Qn,Qp2,maxIter);

    %lambda is right, but the Stats are not directly related to the ones we want.  So we
    %recompute. 
    Cn=1E-14*speye(size(Qn{1},1));
    for idx=1:length(Qn)
        Cn=Cn+Qn{idx}*exp(lambda(idx));
    end
    Cp=1E-14*speye(size(Qp{1},1));
    for idx2=1:length(Qp)
        Cp=Cp+Qp{idx2}*exp(lambda(idx+idx2));
    end
    
    foo=Cp*X';
    iCe=inv(Cn+X*foo);
    C_beta_y=Cp-foo*iCe*foo';%%%OUT OF MEMORY
    Xt_iCe=X'*inv(Cn);
    %The rest will be calculated the same as for the "else" term so we will
    %wait until after the If/End

else

    %%Else, run the normal model

    %% Set up heirarchical model
    Y = [Y; sparse(size(X,2),size(Y,2))];
    X = [X; speye(size(X,2))];

    %Set up the extended covariance model by concatinating the measurement
    %and parameter noise terms
    Q=cell(length(Qn)+length(Qp),1);
    for idx=1:length(Qn)
        Q{idx}=blkdiag(Qn{idx},sparse(size(Qp{1},1),size(Qp{1},2))); % Build block diagonal matrix from Qn & Qp matrices
    end
    for idx2=1:length(Qp)
        Q{idx+idx2}=blkdiag(sparse(size(Qn{1},1),size(Qn{1},2)),Qp{idx2});
    end


    lambda = ones(length(Q),1);  %Initial guess of lambda. 
    tol=1E-4;  %REML goes till tolorence (or max iter) 
    
    t=256;
    dF    = Inf;  %Initial decent 
    cnt=0;  %This is a bookkeeping param for display purposes

    %% This function was adapted from the SPM function spm_reml.m
    %  which is part of the SPM toolkit (http://www.fil.ion.ucl.ac.uk/spm)
    %

    for iter=1:maxIter
        Ce = 1E-14*speye(size(Q{1},1));  %Make sure it stays in numerical precision

        % E-step:

        %Bound lambda to avoid numerircal prec. issues
        lambda=max(lambda,-16);
        lambda=min(lambda,16);

        for i = 1:length(Q)
            Ce = Ce + Q{i}*exp(lambda(i));
        end

        iCe = inv(Ce);
        Xt_iCe = X' * iCe;
        Xt_iCe_X = Xt_iCe * X;
        C_beta_y = inv(Xt_iCe_X);  %Estimate of covariance of beta given the measurements

        
        % M-step:
        P = iCe - (iCe*X)*C_beta_y*Xt_iCe;
        PY=P*Y;
        for i=1:size(Q,1)
            PQ_i{i}=P*Q{i};
        end

        for i = 1:size(Q,1)
            PQ = PQ_i{i};
            PQt=PQ';
            g(i,1) = -0.5*trace(PQ)*exp(lambda(i)) + 0.5*PY'*Q{i}*PY*exp(lambda(i));
            for j = i:size(Q,1)
                PQj = PQ_i{j};
                H(i,j) = -0.5*sum(sum(PQt.*PQj))*exp(lambda(i)+lambda(j));
                H(j,i)=H(i,j);
            end
        end

        %Now update the lambda.  dLambda = -inv(H)*g
        I=eye(size(H,1));
        dL = (expm(H*t) - I)*inv(H)*g;
        lambda = lambda + dL;

        df    = g'*dL;
        if df > dF - exp(-8), t = max(2,t/2); end %retune the regularization if req.
        dF    = df;

        for c=1:cnt, fprintf('\b'); end
        cnt=fprintf('%-5s: %i %5s%e','  ReML Iteration',iter,'...',full(dF));

        if abs(dF) < tol, break; end

    end
end

%Now, put the final pieces together
Beta= C_beta_y * Xt_iCe * Y;


Stats.tstat.beta=Beta;
Stats.tstat.Cov_beta=C_beta_y;
Stats.tstat.dfe=size(X,2);
Stats.tstat.t=Stats.tstat.beta./sqrt(diag(C_beta_y));
%Stats.tstat.pval=2*tcdf(-abs(Stats.tstat.t),Stats.tstat.dfe);

fprintf('\n');

return