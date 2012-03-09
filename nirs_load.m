function [NIRS newNIRSlocation] = nirs_load(NIRSmat,jobNIRSmatCopyChoice,force_redo)
%Load NIRS.mat. If it is to be written in a new location:
%1- if force_redo, then write or possibly overwrite a previous NIRS.mat in that location
%2- if force_redo == 0, then load the file in that location if it exists,
%or write the NIRS.mat there
try
    if isfield(jobNIRSmatCopyChoice,'NIRSmatCopy')
        NewNIRSdir = jobNIRSmatCopyChoice.NIRSmatCopy.NewNIRSdir;
        NewDirCopyNIRS = 1;
    else
        NewDirCopyNIRS = 0;
    end
    if spm_existfile(NIRSmat)
        load(NIRSmat);
        if NewDirCopyNIRS
            [dirN fil1 ext1] =fileparts(NIRSmat);
            dir2 = fullfile(dirN,NewNIRSdir);
            if ~exist(dir2,'dir'), mkdir(dir2); end;
            newNIRSlocation = fullfile(dir2,'NIRS.mat');
            if spm_existfile(newNIRSlocation) && ~force_redo
                load(newNIRSlocation);
            else
                save(newNIRSlocation,'NIRS');
            end
        else
            newNIRSlocation = NIRSmat;
        end
    else
        disp(['Looking for NIRS.mat file that does not exist at location: ' NIRSmat]);
        NIRS = [];
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
end
