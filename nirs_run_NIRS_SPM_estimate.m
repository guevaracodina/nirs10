function out = nirs_run_NIRS_SPM_estimate(job)
%Load NIRS.mat information
%clear NIRS
%load(job.NIRSmat{1,1});
Dmx_files = job.Dmx_files;
%Loop over files
for f1=1:size(Dmx_files,1)
    fname_SPM = Dmx_files{f1};
    clear SPM_nirs
        % NIRS-SPM batch script for 'Estimation' routine, which estimates the GLM
        % parameters and temporal correlation. 
        % fname_SPM : the name of file which results from model specifiction
        % e.g.,'...\NIRS_SPM_v3_1\Sample_data\categorical_indiv_HbO\SPM_indiv_HbO.mat';
        % fname_nirs : the name of the nirs file
        % e.g.,'...\NIRS_SPM_v3_1\Sample_data\converted_NIRS.mat';
        % example usage
        % >> fname_SPM =
        % 'C:\NIRS_SPM_v3_1\Sample_data\categorical_indiv_HbO\SPM_indiv_HbO.mat';
        % >> fname_nirs = 'C:\NIRS_SPM_v3_1\Sample_data\converted_NIRS.mat';
        % >> [SPM_nirs] = estimation_batch(fname_SPM, fname_nirs);

%         try
%             [pathn, name, ext, versn] = fileparts(fname_SPM);
%             pathn = [pathn filesep];
%         catch
%             index = find(fname_SPM == filesep);
%             pathn = fname_SPM(1:pathn(end));
%         end
    %load(fname_nirs);
    load(fname_SPM);
    %load HbO/HbR concentrations
    load(SPM_nirs.nirs.fname);
    switch SPM_nirs.nirs.step
        case 'estimation'
            disp(['The process of estimation has been already done for ' fname_SPM ' Do model specification step again.']);
            %SPM_nirs = [];
            %break;
            
        case 'specification'
            disp(['Begin model parameter estimation for ' fname_SPM]);
            switch SPM_nirs.nirs.Hb
                case 'HbO'
                    Y = nirs_data.oxyData;
                case 'HbR'
                    Y = nirs_data.dxyData;
                case 'HbT'
                    tf = isfield(nirs_data, 'tHbData');
                    if tf == 0
                        Y = nirs_data.oxyData + nirs_data.dxyData;
                    elseif tf == 1
                        Y = nirs_data.tHbData;
                    end
            end
            if isfield(SPM_nirs.xVi, 'V') == 1
                SPM_nirs = rmfield(SPM_nirs, 'xVi');
                [SPM_nirs] = precoloring(SPM_nirs, Y);
            elseif isfield(SPM_nirs.xVi, 'V') == 0
                [SPM_nirs] = prewhitening(SPM_nirs, Y, pathn);
            end
            disp('Completed.');
            save(fname_SPM, 'SPM_nirs');
        otherwise
    end      
end        
%save(job.NIRSmat{1,1},'NIRS');
%out.NIRSmat{1} = fullfile(NIRS.subj_path,'NIRS.mat');
out = [];
end

