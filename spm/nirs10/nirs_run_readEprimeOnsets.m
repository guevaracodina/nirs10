function out = nirs_run_readEprimeOnsets(job)
%
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Mahnoush Amiri
% 2011-07

% Loop over subjects
for iSubj=1:size(job.NIRSmat,1)
    
    % Load NIRS.mat
    %     try
    [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{iSubj,1},job.NIRSmatCopyChoice,job.force_redo);
    job.NIRSmat{iSubj,1} = newNIRSlocation;
    
    % Columns to be read from excelEprime
    excelp = job.excelEprime{:};
    
    %
    %             NIRS.Dt.fir.stax.p{1} = staxp;
    [sDtp,name,ext] = fileparts(job.NIRSmat{iSubj,1});
    
    %             try % try to get columns automatically
    %                 jobH.subj.sDtp = sDtp;
    %                 jobH.subj.helmet.staxp = staxp;
    %                 nirs_criugm_getHelmet_autosave(jobH);
    %             catch % launch GUI for the user to select the right points
    jobH.subj.sDtp = sDtp;
    jobH.subj.excelp = excelp;
    nirs_readEprimeOnsets(jobH); % get informative columns from excelEprime
    fig=gcf;
    waitfor(fig,'BeingDeleted','On');
    
    %             end
    
    % Data from selected columns is saved in a .mat structure which
    % is deleted immediately after updating NIRS
    load(fullfile(sDtp, 'NIRS_eprime.mat'));
    NIRS.Dt.aux.eprime = eprime;
    NIRS.Dt.aux.eprime.p = excelp;
    delete(fullfile(sDtp, 'NIRS_eprime.mat'));
    
    % Get Sess.mat which contains "duration", "name"
    % and "onsets" for all sessions of a subject
    jobS.NIRSmat = job.NIRSmat{iSubj,1};
    jobS.eprime = eprime;
    jobS.excelp = excelp;
    Sess = onsets_eprime(jobS);
    
    NIRS.Dt.fir.Sess = Sess;
    
    
    save(job.NIRSmat{iSubj,1},'NIRS');
    
    %     catch exception
    %         disp(exception{1});
    %         disp(['readEprime failed for the ' int2str(iSubj) 'th subject.']);
    %     end
end
out.NIRSmat = job.NIRSmat; %job.NIRSmat{iSubj};
