function [rendered_MNI run_contrast_OK NIRS] = nirs_load_TopoData(job,NIRS,run_contrast_OK)
%load topographic information (formerly known as preproc_info)
try
    %If present, this loads NIRS_SPM's rendered_MNI structure
    fname_ch = job.TopoData{1}; %if user specified
    if ~isempty(fname_ch)
        load(fname_ch);
        NIRS.Dt.ana.rend = fname_ch; %possibly overwrite previous choice
    else
        if isfield(NIRS.Dt.ana,'rend')
            fname_ch = NIRS.Dt.ana.rend;
            load(fname_ch); %load rendered_MNI
        else
            disp(['No TopoData structure for subject ' int2str(Idx) ...
                '. Rerun coregistration or even first module and make sure that TopoData is generated.']);
            run_contrast_OK = 0;
        end
    end
catch  exception
    disp(exception.identifier);
    disp(exception.stack(1));
    run_contrast_OK = 0;
    rendered_MNI = [];
end