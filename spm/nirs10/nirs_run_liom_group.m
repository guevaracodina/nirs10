function out = nirs_run_liom_group(job)
%try to get factorial design specification
try
    %this is not coded up at all
    out = nirs_spm_run_factorial_design(job);
end
Z = get_contrast_group_common_options(job);
Z.StatStr = 'EC';
%Run simple group level analysis as a one sample t-test
Z.FFX = job.FFX_or_RFX;
Z.p_value = job.contrast_p_value;
Z.LKC = job.StatMethod;
% if ~Z.LKC
%     %at the group level, two ways to get corrected statistics:
%     %1-LKC
%     %2-Bonferroni
%     %currently using output_unc for Bonferroni? That's not clear and that's
%     %confusing
%     Z.output_unc = 1;
%     %Also, if at the contrast module, LKC was run, then some of the older
%     %options will not be available
% end
Z.group_session_to_average = job.group_session_to_average;
Z.group_dir_name = job.group_dir_name;
Z.simple_sum = job.simple_sum;
Z.min_s = 2;
Z.nS = size(job.NIRSmat,1);
nS = Z.nS;
try 
    number_dir_to_remove = job.number_dir_to_remove;
catch
    number_dir_to_remove = 3;
end
%SPM contrasts or automatic contrasts
if isfield(job.ContrastChoice,'user_contrasts')
    if isempty(job.ContrastChoice.user_contrasts.consess)
        disp('No contrast specified. Aborting. Either specify contrasts or use automated-contrasts option');
        return
    else
        Z.automated_contrasts = 0;
        TF.consess = job.ContrastChoice.user_contrasts.consess;
    end
else
    TF.consess = [];
    %automated contrasts if no user-specified contrast
    Z.automated_contrasts = 1;
end
if Z.FFX || nS==1
    %fixed effects: loop over subjects first, as they are treated
    %separately
    nl = nS; %to loop over subjects
    big_TOPO = []; %not used
else
    %RFX - loop over subjects done later
    nl = 1;
    big_TOPO{nS} =[];
    for Idx=1:nS
        %Load NIRS.mat information
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        dir1 = NIRS.SPM{1};
        %load topographic information (formerly known as preproc_info)
        if Idx == 1 %assume same configuration for each subject - could be generalized
            fname_ch = NIRS.Dt.ana.rend;
            load(fname_ch);
        end
        try
            ftopo = NIRS.TOPO;
        catch
            ftopo = fullfile(dir1,'TOPO.mat');
            disp('TOPO not found in NIRS.TOPO');
        end
        TOPO = [];
        load(ftopo);
        %large structure
        big_TOPO{Idx} = TOPO;
    end
    %create a new TOPO at the group level
    TOPO = [];
    TOPO.xCon = big_TOPO{1}.xCon;
end


%Loop over all subjects for FFX but go through only once for RFX
for Idx=1:nl
    %Load NIRS.mat information
    try
        if Z.FFX || nS==1
            NIRS = [];
            load(job.NIRSmat{Idx,1});
            %load topographic information (formerly known as preproc_info)
            fname_ch = NIRS.Dt.ana.rend;
            load(fname_ch);
            ftopo = NIRS.TOPO; 
            [dir1 fil1 ext1] = fileparts(ftopo);
            TOPO = [];
            load(ftopo);
            dir_group = fullfile(dir1, Z.group_dir_name);
        else
            [dir0,dummy,dummy2] = fileparts(job.NIRSmat{1});
            %extract previous directory
            tmp = strfind(dir0,filesep);
            dir_root = dir0(1:tmp(end-number_dir_to_remove));
            dir_group = fullfile(dir_root, Z.group_dir_name);
        end
        if ~exist(dir_group,'dir'), mkdir(dir_group); end
        %store in same directory as first subject
        ftopo = fullfile(dir_group,'TOPO.mat');
        %save a NIRS structure for the group
        newNIRSlocation = fullfile(dir_group,'NIRS.mat');
        NIRS.TOPO = ftopo;
        save(newNIRSlocation,'NIRS');
        job.NIRSmat{nl,1} = newNIRSlocation;
        Z.dir1 = dir_group;
        %contrasts
        TOPO = nirs_get_contrasts_group(Z,TF,TOPO);
        if Z.FFX || nS==1
            fg = 'g';
            TOPOsrc = TOPO;
        else
            fg = 'group';
            TOPOsrc = big_TOPO;
        end
        load Split
        F.split = split;
        F.pathn = Z.dir1;
        %CF: copy figure structure
        CF.GInv = Z.GInv;
        CF.split = split;
        xCon = TOPO.xCon;
        nC = length(xCon);
        CF.nC = nC;
        %Big loop over views
        for v1=1:6
            view_estimated = 0;
            try
                if isfield(TOPOsrc.v{v1},'s1')
                    ns = length(TOPOsrc.v{v1}.s); %number of sessions
                    view_estimated = 1;
                end
            end
            try
                if isfield(TOPOsrc{1}.v{v1},'s1')
                    ns = length(TOPOsrc); %number of subjects
                    view_estimated = 1;
                end
            end
            
            if view_estimated   
                %Structure for passing GLM and interpolation data
                W = [];
                [W.side_hemi W.spec_hemi] = nirs_get_brain_view(v1);
                %View dependent info for figures
                brain = rendered_MNI{v1}.ren;
                if issparse(brain), %does not apply?
                    d = size(brain);
                    B1 = spm_dctmtx(d(1),d(1));
                    B2 = spm_dctmtx(d(2),d(2));
                    brain = B1*brain*B2';
                end;
                msk = brain>1;brain(msk)=1;
                msk = brain<0;brain(msk)=0;
                %brain = brain(end:-1:1,:); %???
                brain = brain * 0.5;
                W.brain = brain;
                W.s1 = size(brain, 1);
                W.s2 = size(brain, 2);
                TOPO.v{v1}.(fg).ns = ns;
                TOPO.v{v1}.(fg).min_s = Z.min_s;
                TOPO.v{v1}.(fg).s1 = W.s1;
                TOPO.v{v1}.(fg).s2 = W.s2;
            
                %Handles for assembled figures
                H = initialize_assembled_figure_handles;
                H = initialize_assembled_figures(Z,H,0,'Group');
                %Loop over chromophores
                for h1=1:3 %including HbT
                    for c1=1:nC %Loop over contrasts
                        %Positive stats
                        [H,TOPO,big_TOPO] = fill_group(H,TOPO,big_TOPO,v1,c1,h1,Z,W,F,CF,xCon,ns,1);
                        %Negative stats
                        [H,TOPO,big_TOPO] = fill_group(H,TOPO,big_TOPO,v1,c1,h1,Z,W,F,CF,xCon,ns,0);
                    end
                end                
                call_save_assembled_figures(Z,W,H,0);
            end %if view_estimated
        end %end for v1
        
        if exist('Sess','var')
            TOPO.Sess = Sess;
            TOPO.Cp = Cp;
        end
        save(ftopo,'TOPO','-v7.3');
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp('Could not do group analysis');
    end
end
out.NIRSmat = job.NIRSmat;