function TOPO = prepare_constrast_core_call(Z,W,SPM,TOPO)
xXn = SPM.xXn;
NC = SPM.xY.Cf;
Z.Avg = 0; %This sets whether we are looking at averaged rather than GLM data
%Objective is to specify W.beta, W.res, W.var and W.corr_beta
%then pass that to contrast_core
try
    %loop over 2 options
    if Z.GroupMultiSession
        %Sessions are processed as one
        Nsess = 1;
    else
        %Sessions are processed individually
        Nsess = length(xXn);
    end
    for f1=1:Nsess
        if any(Z.sessions == 0) || any(f1 == Z.sessions)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %W.beta, W.res, W.var and W.corr_beta
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if Z.GroupMultiSession
                W.beta = [];
                W.res = [];
                W.var = [];
                W.var0 = [];
                W.corr_beta = [];
                %pile-up all the betas from all the sessions and all channels
                for f2 = 1:length(xXn)
                    if any(Z.sessions == 0) || any(f2 == Z.sessions)
                    W.beta = [W.beta; xXn{f2}.beta];
                    if isfield(xXn{f2},'res')
                        tmp_res = fopen_NIR(xXn{f2}.res,NC);
                        W.res = [W.res tmp_res];
                    else
                        Z.Avg = 1;
                        if f2 == 1
                            W.covbeta = xXn{f2}.covbeta;
                        else
                            W.covbeta = [W.covbeta;  xXn{f2}.covbeta];
                        end
                    end
                    %How to calculate those? for now, just take the average
                    if f2 == 1
                        tmp_var = xXn{f2}.ResSS./xXn{f1}.trRV;
                        tmp_varch = xXn{f2}.ResSSch./xXn{f1}.trRV;
                        tmp_corr = xXn{f2}.Bcov;
                    else
                        tmp_var = tmp_var + xXn{f2}.ResSS./xXn{f1}.trRV;
                        tmp_varch = tmp_varch + xXn{f2}.ResSSch./xXn{f1}.trRV;
                        tmp_corr = blkdiag(tmp_corr, xXn{f2}.Bcov);
                    end
                    end
                end
                if ~Z.Avg
                    W.var = tmp_var/length(xXn);
                    W.varch = tmp_varch/length(xXn);
                    W.corr_beta = tmp_corr;
                end
            else
                W.beta = xXn{f1}.beta;
                if isfield(xXn{f1},'res')
                    W.res = fopen_NIR(xXn{f1}.res,NC);
                else
                    Z.Avg = 1;
                    W.covbeta = xXn{f1}.covbeta;
                end
                try
                    if ~Z.Avg
                        %for NIRS_SPM method
                        W.var = xXn{f1}.ResSS./xXn{f1}.trRV;
                        W.varch = xXn{f1}.ResSSch./xXn{f1}.trRV; %channel by channel
                        %covariance of beta estimates
                        W.corr_beta = xXn{f1}.Bcov;
                    end
                catch exception
                    try
                        %for WLS and BGLM methods
                        W.corr_beta = xXn{f1}.Bvar;
                    catch exception2
                        disp(exception2.identifier);
                        disp(exception2.stack(1));
                        disp(exception.identifier);
                        disp(exception.stack(1));
                        disp('Problem with Bcov or ResSS or trRV if using NIRS_SPM, or Bvar if using WLS/BGLM');
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if Z.GroupMultiSession
                TOPO = constrasts_core(Z,W,TOPO,TOPO.xX,0);
            else
                TOPO = constrasts_core(Z,W,TOPO,xXn{f1},f1);
            end
        end
    end
catch exception
    disp(exception.identifier)
    disp(exception.stack(1));
    disp('Problem in prepare_contrast_call');
end