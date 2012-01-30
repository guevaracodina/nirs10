function sC = nirs_nonlinear_epilepsy_contrasts(SPM,nr)
contrastT = [];
contrastF = [];
contrastT_name = [];
contrastF_name = [];
try
    %number of confounds -- no longer used?
    ncf = NIRSconfounds.NumChConfoundsActual;
catch
    ncf = 0;
end
%Need to check if GLM was run with or without derivs
switch SPM.xBF.name
    case {'hrf','Gamma functions'}        
        %canonical HRF - no derivatives - no F contrasts
        %1st Volterra kernel
        contrastT{1} = [1 zeros(1,nr-1)];
        contrastT_name{1} = 'T1';
        %2nd Volterra kernel
        if SPM.job.volt > 1
            if nr > 5+ncf %Careful, this might not be the correct number
                %if there are more confounding regressors
                %assume 2 stimuli - only take the first one
                contrastT{2} = [0 0 1 zeros(1,nr-3)];
                contrastT_name{2} = 'T2';
            else
                %assume only 1 stimulus
                contrastT{2} = [0 1 zeros(1,nr-2)];
                contrastT_name{2} = 'T2';
            end
        end
        if SPM.job.volt == 3
            if nr > 5+ncf
                %assume 2 stimuli - only take the first one
                contrastT{3} = [0 0 0 0 0 0 1 zeros(1,nr-7)];
                contrastT_name{3} = 'T3';
            else
                %assume only 1 stimulus
                contrastT{3} = [0 0 1 zeros(1,nr-3)];
                contrastT_name{3} = 'T3';
            end
        end
    case 'hrf (with time derivative)'
        %Automated T contrasts
        %1st Volterra kernel
        contrastT{1} = [1 zeros(1,nr-1)];
        contrastT_name{1} = 'T1';
        %2nd Volterra kernel
        if SPM.job.volt > 1
            if nr > 10+ncf %Careful, this might not be the correct number
                %if there are more confounding regressors
                %assume 2 stimuli - only take the first one
                contrastT{2} = [0 0 0 0 1 zeros(1,nr-7)];
                contrastT_name{2} = 'T2';
            else
                %assume only 1 stimulus
                contrastT{2} = [0 0 1 zeros(1,nr-4)];
                contrastT_name{2} = 'T2';
            end
        end
        %Automated F contrasts - careful, F contrasts need to be a cell
        contrastF{1} = {[eye(2) zeros(2,nr-2)]};
        contrastF_name{1} = 'F1';
        if SPM.job.volt > 1
            %temp_mat = [ zeros(2,1) eye(2) zeros(2,3); zeros(1,5) 1];
            mat_Volt = [eye(3)]; %[eye(3) zeros(3,6); zeros(3,3) temp_mat];
            if nr > 10+ncf
                contrastF{2} = {[zeros(3,4) mat_Volt zeros(3,nr-7)]};
                contrastF_name{2} = 'F2';
            else
                contrastF{2} = {[zeros(3,2) mat_Volt zeros(3,nr-5)]};
                contrastF_name{2} = 'F2';
            end
        end
        
        %Automated F contrasts
    case 'hrf (with time and dispersion derivatives)'
        %Automated T contrasts
        %1st Volterra kernel
        contrastT{1} = [1 zeros(1,nr-1)];
        contrastT_name{1} = 'T1';
        %2nd Volterra kernel
        if SPM.job.volt > 1
            if nr > 12+ncf %Careful, this might not be the correct number
                %if there are more confounding regressors
                %assume 2 stimuli - only take the first one
                contrastT{2} = [0 0 0 0 0 0 1 zeros(1,nr-7)];
                contrastT_name{2} = 'T2';
            else
                %assume only 1 stimulus
                contrastT{2} = [0 0 0 1 zeros(1,nr-4)];
                contrastT_name{2} = 'T2';
            end
        end
        %Automated F contrasts - careful, F contrasts need to be a cell
        contrastF{1} = {[eye(3) zeros(3,nr-3)]};
        contrastF_name{1} = 'F1';
        if SPM.job.volt > 1
            %temp_mat = [ zeros(2,1) eye(2) zeros(2,3); zeros(1,5) 1];
            mat_Volt = eye(6); %[eye(3) zeros(3,6); zeros(3,3) temp_mat];
            if nr > 12+ncf
                contrastF{2} = {[zeros(6,6) mat_Volt zeros(6,nr-12)]};
                contrastF_name{2} = 'F2';
            else
                contrastF{2} = {[zeros(6,3) mat_Volt zeros(6,nr-9)]};
                contrastF_name{2} = 'F2';
            end
        end
end
sC.contrastT = contrastT;
sC.contrastF = contrastF;
sC.contrastT_name = contrastT_name;
sC.contrastF_name = contrastF_name;