function out = nirs_run_ReMLreconstruct(job)
% Achieve image segmentation after New Segment
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% from Ted Huppert ReML code

% Clément Bonnéry
% 2010-09

NIRS = job.NIRS;
% 
% NIRS.nirsfile = job.nirsfile
% nirsfile = load('NIRS.nirsfile','-mat');

% % % % covariances inter and intra NIRS signals
% % %
% % % % intra (temporal covariance)
% % % %-- sur OD : sur chaque paire
% % % % inter (spatial covariance)
% % % %-- on calcule les covariances entre les OD qui constituent les signaux
% % % %directement enregistres sur les paires
% % %
% % % %-- est ce qu'on devrait calculer la correlation entre HbO et HbR ????????
% % % %(en prenant soin de correler HbO et -HbR disons...)

n_pairs = size(NIRS.nirs_file.d,2)/2;
% SAME wavelength
for n_wl = 1:2
    for i = 1:n_pairs
        for j = 1:n_pairs
            Q{(n_wl-1)*n_pairs+i,(n_wl-1)*n_pairs+j} = xcorr(NIRS.nirs_file.d(:,(n_wl-1)*n_pairs+i),NIRS.nirs_file.d(:,(n_wl-1)*n_pairs+j));
            Q{i,(n_wl-1)*n_pairs+j} = xcorr(NIRS.nirs_file.d(:,i),NIRS.nirs_file.d(:,n_pairs+j));
            Q{(n_wl-1)*n_pairs+i,j} = xcorr(NIRS.nirs_file.d(:,n_pairs+i),NIRS.nirs_file.d(:,j));
        end
    end
end

[C,h,Ph,F,Fa,Fc] = spm_reml_sc(YY,X,Q,NIRS.nirs_file.size(d,1),-32,256,V);
% ReML estimation of covariance components from y*y' - proper components
% FORMAT [C,h,Ph,F,Fa,Fc] = spm_reml_sc(YY,X,Q,N,[hE,hC,V]);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components
% N   - number of samples
%
% hE  - hyperprior expectation in log-space [default = -32]
% hC  - hyperprior covariance  in log-space [default = 256]
% V   - fixed covariance component
%
% C   - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) ReML hyperparameters h
% Ph  - (q x q) conditional precision of log(h)
%
% F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective
%
% Fa  - accuracy
% Fc  - complexity (F = Fa - Fc)
%
% Performs a Fisher-Scoring ascent on F to find MAP variance parameter
% estimates.  NB: uses weakly informative log-normal hyperpriors.
% See also spm_reml for an unconstrained version that allows for negative
% hyperparameters
%
%__________________________________________________________________________
%      spm_reml_sc: positivity constraints on covariance parameters

% on resout à t donné

%get anatomical and functional datas
%-> anatomical datas : 5 layer segmented image



%-> functional datas :
%   -- position of sources and detectors
%   -- HbO and HbT [in SD pairs domain]
%   -- HbO and HbT hemodynamic response ?? (prior temporel...)
%   -- ??
out{1} =1;
return

% function out = nirs_run_reconstruction(job)
% % Achieve image segmentation after New Segment
% %_______________________________________________________________________
% % Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% % from Ted Huppert ReML code
%
% % Clément Bonnéry
% % 2010-09
%
%
% if(~exist('K'))
%     K=35;  %Max # of iterations of REML code
% end
%
% % on resout à t donné
%
% %get anatomical and functional datas
% %-> anatomical datas : 5 layer segmented image
% %-> functional datas :
% %   -- position of sources and detectors
% %   -- HbO and HbT [in SD pairs domain]
% %   -- HbO and HbT hemodynamic response ?? (prior temporel...)
% %   -- ??
%
% %set the model (write down explicitly the system)
% % From SPM book p 283______________________________________________________
%
% % setting of the model : augmentation to embody priors in error covariance
% X = [X,X;speye(size(X,1))];
% y = [y;zeros(size(X,2));eta_beta];
% C_beta = sparse(i,j,C_beta,m,n);
% C_epsilon = C_beta;
%
% lambda = ones(length(Q),1);  %Initial guess of lambda.
% tol = 1E-4; %REML goes till tolerance (or max iter)
% t = 256;
% dF = Inf;   %Initial decent
% cnt = 0;    %This is a bookkeeping param for display purposes
%
%
% % building Q matrixes______________________________________________________
%
%
% % until convergence EM algorithm___________________________________________
% for iter=1:K
%
%     % E-Step :
%     %______________________________________________________________________
%     %   -- C_epsilon
%     for i = 1:length(Q)
%         C_epsilon = C_epsilon + Q{i}*lambda(i);
%     end
%     %   -- C_beta|y
%     iC_epsilon = inv(C_epsilon);
%     Xt_iC_epsilon = X' * iC_epsilon;
%     Xt_iC_epsilon_X = Xt_iC_epsilon * X;
%     C_beta_y = inv(Xt_iC_epsilon_X);  %Estimate of covariance of beta given the measurements
%     % -- n_beta_y : jamais utilise...........
%
%     % M-Step :
%     %______________________________________________________________________
%     %   -- P
%     P = iC_epsilon - (iC_epsilon*X)*C_beta_y*Xt_iC_epsilon;
%     %   -- g_i and H_ij
%     PY=P*Y;
%     for i=1:size(Q,1)
%         PQ_i{i}=P*Q{i};
%     end
%
%     for i = 1:size(Q,1)
%         PQ = PQ_i{i};
%         PQt=PQ';
%         g(i,1) = -0.5*trace(PQ)*exp(lambda(i)) + 0.5*PY'*Q{i}*PY*exp(lambda(i));
%         for j = i:size(Q,1)
%             PQj = PQ_i{j};
%             H(i,j) = -0.5*sum(sum(PQt.*PQj))*exp(lambda(i)+lambda(j));
%             H(j,i)=H(i,j);
%         end
%     end
%
%     %Now update the lambda.  dLambda = -inv(H)*g
%     %     I=eye(size(H,1));
%     dlambda = -inv(H)*g;%(expm(H*t) - I)*inv(H)*g;
%     lambda = lambda + dlambda;
%
%     df    = g'*dlambda;
%     if df > dF - exp(-4), t = max(2,t/2); end %retune the regularization if req.????????????????
%     dF    = df;
%
%     for c=1:cnt, fprintf('\b'); end
%     cnt=fprintf('%-5s: %i %5s%e','  ReML Iteration',iter,'...',full(dF));
%
%     if dF < tol, break; end
%
% end
%
% %Now, put the final pieces together
% Beta= C_beta_y * Xt_iCe * Y;
%
% fprintf('\n');
% out{1} =1;
% return