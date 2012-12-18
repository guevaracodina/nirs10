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
number_dir_to_remove = job.number_dir_to_remove;
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
else
    %RFX - loop over subjects done later
    nl = 1;
end

%Loop over all subjects for FFX but go through only once for RFX
for Idx=1:nl
    %Load NIRS.mat information
    try
        if Z.FFX || nS==1
            [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
            job.NIRSmat{Idx,1} = newNIRSlocation;
            %load topographic information (formerly known as preproc_info)
            fname_ch = NIRS.Dt.ana.rend;
            load(fname_ch);
            rendered_MNI0 = rendered_MNI;
            ftopo = NIRS.TOPO;
            [dir1 fil1 ext1] = fileparts(ftopo);
            TOPO = [];
            load(ftopo);
            TOPO.rendered_MNI = rendered_MNI;
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
        if ~(Z.FFX || nS==1)
            job.NIRSmat{nl,1} = newNIRSlocation;
            try
                [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
            end
        end
        job.NIRSmat{Idx,1} = newNIRSlocation;
        
        NIRS.TOPO = ftopo;
        save(newNIRSlocation,'NIRS');
        if ~isfield(NIRS,'flags')
            NIRS.flags = [];
        end
        %Could save time by making these checks before loading all the
        %group data
        if (( (Z.FFX || nS==1) && ~isfield(NIRS.flags,'session_groupOK') || job.force_redo) || ...
                (~(Z.FFX || nS==1) && ~isfield(NIRS.flags,'groupOK') || job.force_redo))
            %disp('Starting group analysis');
            if Z.FFX || nS==1
                big_TOPO = [];
            else
                big_TOPO{nS} = [];
                disp('Wait... loading individual data for group analysis');
                for Idx2=1:nS
                    %Load NIRS.mat information
                    NIRS = [];
                    load(job.NIRSmat{Idx2,1});
                    dir1 = NIRS.SPM{1};
                    %load topographic information (formerly known as preproc_info)
                    fname_ch = NIRS.Dt.ana.rend;
                    %quick fix for Claudine's study:
                    %fname_ch = 'W:\Claudine\SPMDataNT\S003\TopoData.mat';
                    
                    load(fname_ch);
                    if Idx == 1
                        rendered_MNI0 = rendered_MNI;
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
                    big_TOPO{Idx2} = TOPO;
                    big_TOPO{Idx2}.rendered_MNI = rendered_MNI;
                end
                %create a new TOPO at the group level
                TOPO = [];
                TOPO.xCon = big_TOPO{1}.xCon;
            end
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
                        if (Z.FFX || nS==1) && ns == 1
                            view_estimated = 0; %force not running the group if only one subject and only one session
                        else
                            view_estimated = 1;
                        end
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
                    brain = rendered_MNI0{v1}.ren;
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
                    %maximal number of contrasts to group in assembled
                    %figures
                                      
                    Z.nCloop = 4;
                    for c1=1:nC %Loop over contrasts
                        %effective c1 contrast number for the purpose of
                        %assembling figures only
                        Z.c1eff = mod(c1,Z.nCloop);
                        if Z.c1eff == 0
                            Z.c1eff = Z.nCloop;
                        end
                        if Z.c1eff == 1
                            Z.scon = '';
                            Z.sconFig = '';
                            for k0=c1:(c1+Z.nCloop-1)
                                if k0 <= nC
                                    Z.scon = [Z.scon '_' xCon(k0).name];
                                    %nrep = regexprep(xCon(k0).name, '_', ' ');
                                    %Z.sconFig = [Z.sconFig ' ' nrep];       
                                end
                            end
                            %Handles for assembled figures
                            Z.spec_hemi = W.spec_hemi;
                            CF.Z = Z; %horrible, but needed for nirs_copy_figure (near line 17)
                            H = initialize_assembled_figure_handles;
                            H = initialize_assembled_figures(Z,H,0,'Group');
                        end
                        
                        %Loop over chromophores
                        for h1=1:3 %including HbT                        
                            %Positive stats
                            [H,TOPO,big_TOPO] = fill_group(H,TOPO,big_TOPO,v1,c1,h1,Z,W,F,CF,xCon,ns,1);
                            %Negative stats
                            [H,TOPO,big_TOPO] = fill_group(H,TOPO,big_TOPO,v1,c1,h1,Z,W,F,CF,xCon,ns,0);
                        end
                        %save after each group or when reaching the last
                        %contrast
                        if Z.c1eff == Z.nCloop || c1 == nC
                            call_save_assembled_figures(Z,W,H,0);
                        end
                    end
                    
                end %if view_estimated
            end %end for v1
            
            if exist('Sess','var')
                TOPO.Sess = Sess;
                TOPO.Cp = Cp;
            end
            save(ftopo,'TOPO','-v7.3');
            %save NIRS
            if Z.FFX || nS==1
                NIRS.flags.session_groupOK = 1;
            else
                NIRS.flags.groupOK = 1;
            end
            save(newNIRSlocation,'NIRS');
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp('Could not do group analysis');
    end
end
out.NIRSmat = job.NIRSmat;