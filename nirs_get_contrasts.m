function TOPO = nirs_get_contrasts(SPM0,Z,TF,TOPO)
%Organization. Either:
%1- the user has specified contrasts (TF not empty) or
%2- contrasts need to be automatically generated (TF empty), for
%each session (sessions cannot be grouped)

%Make sure all the outputs are assigned:
TOPO.xCon = [];
%SPM and xCon are for contrasts at the subject level
%SS_SPM and SSxCon are arrays for each of the sessions
try
    SPM = SPM0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~Z.automated_contrasts && Z.GroupMultiSession
        %1- get spm contrasts, which are put into SPM.xCon, using structure TF

        %Construct the full design matrix over all sessions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp1 = []; tmp2 = []; tmp3 = 0; tmp4 = 0; tmp5 = 0; tmp6 = 0; tmp7 = []; tmp8 = [];
        for s1=1:length(SPM.xXn)
            if Z.sessions == 0 || any(s1 == Z.sessions)
                %filtered design matrix
                tmp1 = blkdiag(tmp1,SPM.xXn{s1}.xKXs.X);
                %design matrix
                tmp2 = blkdiag(tmp2,SPM.xXn{s1}.X);
                %residual sum of squares of full model
                tmp3 = tmp3 + SPM.xXn{s1}.ResSS;
                tmp8 = tmp8 + SPM.xXn{s1}.ResSSch;
                tmp5 = tmp5 + SPM.xXn{s1}.trRV;
                tmp6 = tmp6 + SPM.xXn{s1}.trRVRV;
                tmp7 = blkdiag(tmp7,SPM.xXn{s1}.Bcov);
            end
        end
        SPM.xX.xKXs = tmp1;
        SPM.xX.X = tmp2;
        SPM.xX.ResSS = tmp3;
        SPM.xX.ResSSch = tmp8;
        %Not clear from Satterthwaite approximation whether
        %to use tmp5/tmp6 or tmp4 for # of degrees of freedom
        %for multiple sessions
        SPM.xX.trRV = tmp5;
        SPM.xX.trRVRV = tmp6;
        SPM.xX.erdf = tmp5^2/tmp6;
        SPM.xX.erdf_check = tmp4;
        SPM.xX.corr_beta = tmp7;
        SPM.xX.var = SPM.xX.ResSS/SPM.xX.trRV;
        SPM.xCon = [];
        %Generate the SPM-type contrasts over all sessions
        disp('Subject-level contrasts:');
        SPM = nirs_spm_run_con(TF,SPM);
        TOPO.xX = SPM.xX;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else %if Z.automated_contrasts || ~Z.GroupMultiSession
        TF0 = TF;
        %loop over selected sessions
        C_done = 0;
        for s1=1:length(SPM.xXn)
            if Z.sessions == 0 || any(s1 == Z.sessions)
                if Z.automated_contrasts
                    TF = TF0;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %2- Generate automated contrasts if required
                    %number of regressors for one session
                    nr = size(SPM.xXn{s1}.X,2);
                    sC.contrastT = {}; sC.contrastF = {};
                    sC.contrastT_name = {}; sC.contrastF_name = {};
                    if Z.NonlinearEpilepsyOn
                        sC = nirs_nonlinear_epilepsy_contrasts(SPM,nr);
                    else
                        for i0=1:nr-1
                            sC.contrastT{i0} = [zeros(1,i0-1) 1 zeros(1,nr-i0)];
                            try
                                sC.contrastT_name{i0} = ['_' validate_name(SPM.Sess(1).U(i0).name{1})];
                            catch %exception
                                %disp(exception.identifier);
                                %disp('Using default names for contrasts');
                                sC.contrastT_name{i0} = ['C' int2str(i0)];
                            end
                        end
                    end
                    for j1=1:length(sC.contrastT)
                        TF.consess{j1}.tcon.name = sC.contrastT_name{j1};
                        TF.consess{j1}.tcon.convec = sC.contrastT{j1};
                        TF.consess{j1}.tcon.sessrep = 'none';
                    end
                    for j1=1:length(sC.contrastF)
                        TF.consess{j1+length(contrastT)}.fcon.name = sC.contrastF_name{j1};
                        TF.consess{j1+length(contrastT)}.fcon.convec = sC.contrastF{j1};
                        TF.consess{j1+length(contrastT)}.fcon.sessrep = 'none';
                    end
                end
                if ~C_done
                    SPM = SPM0;
                    SPM.xCon = [];
                    SPM.xX = SPM.xXn{s1};
                    %contrasts get added to previously generated contrasts stored
                    %in SPM structure
                    disp('Subject-level contrasts:');
                    SPM = nirs_spm_run_con(TF,SPM);
                    C_done = 1;
                end
                %Need the single-session contrasts to cope with the
                %situation where the contrasts are different from session
                %to session -- this happens a lot in epilepsy, but also in
                %cognitive studies (Maude study)
                SPM1 = SPM0;
                SPM1.xCon = [];
                SPM1.xX = SPM1.xXn{s1};
                disp(['Contrasts for session ' int2str(s1) ':']);
                SPM1 = nirs_spm_run_con(TF,SPM1);
                TOPO.SSxCon{s1} = SPM1.xCon;
            end
        end
    end
    TOPO.xCon = SPM.xCon;
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Could not generate contrasts');
end