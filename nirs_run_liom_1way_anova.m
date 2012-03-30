function out = nirs_run_liom_1way_anova(job)
%try to get factorial design specification
try
    %this is not coded up at all
    out = nirs_spm_run_factorial_design(job);
end
Z = get_contrast_group_common_options(job);
%Get anova info
Z.anova_level = job.anova_level;
for l1=1:Z.anova_level
    level_name{l1} = job.level(l1).level_name;
    level_subj{l1} = job.level(l1).level_subj;
end
%Run simple group level analysis as a one sample t-test
Z.group_session_to_average = job.group_session_to_average;
%Generate contrasts for inverted responses
Z.GInv = 1;
Z.GFIS = 1;
Z.anova_dir_name = job.anova_dir_name;
number_dir_to_remove = job.number_dir_to_remove;
%Structure for passing more generic data
min_s = 2;
Z.StatStr = 'EC';
Z.LKC = job.StatMethod;
Z.CorrectionMethod = job.CorrectionMethod;
nS = size(job.NIRSmat,1);
%RFX - loop over subjects done later
nl = 1;
big_TOPO{nS} =[];
for Idx=1:nS
    %Load NIRS.mat information
    NIRS = [];
    load(job.NIRSmat{Idx,1});
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
    end
    TOPO = [];
    load(ftopo);
    %large structure
    big_TOPO{Idx} = TOPO;
    big_TOPO{Idx}.rendered_MNI = rendered_MNI;
end
%create a new TOPO at the group level
TOPO = [];

%Load NIRS.mat information
try
    [dir0,dummy,dummy2] = fileparts(job.NIRSmat{1});
    %extract previous directory
    tmp = strfind(dir0,filesep);
    dir_root = dir0(1:tmp(end-number_dir_to_remove));
    dir_group = fullfile(dir_root,Z.anova_dir_name);
    if ~exist(dir_group,'dir'), mkdir(dir_group); end
    %store in same directory as first subject
    ftopo = fullfile(dir_group,'TOPO.mat');
    %save a NIRS structure for the group
    newNIRSlocation = fullfile(dir_group,'NIRS.mat');
    NIRS.TOPO = ftopo;
    save(newNIRSlocation,'NIRS');
    job.NIRSmat{nl,1} = newNIRSlocation;
    Z.dir1 = dir_group;    
    %Big loop over views
    for v1=1:6
        view_estimated = 0;
        try
            if isfield(big_TOPO{1}.v{v1},'s1')
                s1 = big_TOPO{1}.v{v1}.s1;
                s2 = big_TOPO{1}.v{v1}.s2;
                ns = nS; %number of subjects
                view_estimated = 1;
            end
        catch
            view_estimated = 0;
        end
        
        if view_estimated
            [side_hemi spec_hemi] = nirs_get_brain_view(v1);
            %View dependent info for figures
            %brain = rend{v1}.ren;
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
            
            %Structure for passing GLM and interpolation data
            W = [];
            W.brain = brain;
            W.s1 = size(brain, 1);
            W.s2 = size(brain, 2);
            W.spec_hemi = spec_hemi;
            W.side_hemi = side_hemi;
            %Contrasts
            TOPO.v{v1}.group.ns = ns;
            TOPO.v{v1}.group.min_s = min_s;
            TOPO.v{v1}.group.s1 = s1;
            TOPO.v{v1}.group.s2 = s2;
            %Contrasts -- assume same contrasts for all subjects
            try
                xCon = big_TOPO{1}.SSxCon; % big_TOPO{1}.xCon}; %????
            catch
                xCon = big_TOPO{1}.xCon;
            end
            %Contrasts - assume same contrasts for all subjects
            %xCon = big_TOPO{1}.xCon;
            cbeta = zeros(ns,s1*s2);
            ccov_beta = zeros(ns,s1*s2);
            tmp = zeros(s1,s2);
            nC = length(xCon);
            H = initialize_assembled_figure_handles;
            H = initialize_assembled_figures(Z,H,0,'Group');
            load Split
            F.split = split;
            F.pathn = Z.dir1;
            %CF: copy figure structure
            CF.GInv = Z.GInv;
            CF.split = split;
            CF.nC = nC;
            
            %Loop over chromophores
            for h1=1:3 %including HbT
                hb = get_chromophore(h1);
                for c1=1:nC
                    try
                        %Skip F contrasts for now
                        if xCon{Z.group_session_to_average}(c1).STAT == 'T'
                            %fill in cbeta and ccov_beta
                            for f1=1:ns
                                try
                                    if ~isfield(big_TOPO{f1}.v{v1},'s')
                                        %group analysis of a group of
                                        %sessions analysis
                                        if isfield(big_TOPO{f1}.v{v1}.g.hb{h1},'beta_map')
                                            tmp = squeeze(big_TOPO{f1}.v{v1}.g.hb{h1}.beta_map(c1,:,:));
                                            cbeta(f1,:) = tmp(:);
                                            tmp = (tmp./squeeze(big_TOPO{f1}.v{v1}.g.hb{h1}.stat_map(c1,:,:))).^2;
                                            ccov_beta(f1,:) = tmp(:);
                                        else
                                            tmp = squeeze(big_TOPO{f1}.v{v1}.g.hb{h1}.c_interp_beta(c1,:,:));
                                            cbeta(f1,:) = tmp(:);
                                            tmp = squeeze(big_TOPO{f1}.v{v1}.g.hb{h1}.c_cov_interp_beta(c1,:,:));
                                            ccov_beta(f1,:) = tmp(:);
                                        end
                                    else
                                        %for is1=1:length(big_TOPO{f1}.v{v1}.s)
                                        is1 = Z.group_session_to_average;
                                        if isfield(big_TOPO{f1}.v{v1}.s{is1}.hb{h1},'beta_map')
                                            %do each session separately
                                            tmp = squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.beta_map(c1,:,:));
                                            cbeta(f1,:) = tmp(:);
                                            tmp = (tmp./squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.stat_map(c1,:,:))).^2;
                                            ccov_beta(f1,:) = tmp(:);
                                        else
                                            %do each session separately
                                            tmp = squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_interp_beta(c1,:,:));
                                            cbeta(f1,:) = tmp(:);
                                            tmp = squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_cov_interp_beta(c1,:,:));
                                            ccov_beta(f1,:) = tmp(:);
                                        end
                                    end
                                catch
                                    disp(['No data for subject ' int2str(f1) ' contrast ' int2str(c1)  ' chromophore ' hb ' and view ' int2str(v1)]);
                                end
                            end
                        else %quick fix for F stats
                            for f1=1:ns
                                try
                                    if ~isfield(big_TOPO{f1}.v{v1},'s')
                                        %group analysis of a group of
                                        %sessions analysis
                                        tmp = squeeze(big_TOPO{f1}.v{v1}.g.hb{h1}.c_interp_F(c1,:,:));
                                        cbeta(f1,:) = tmp(:);
                                        sz = size(squeeze(big_TOPO{f1}.v{v1}.g.hb{h1}.c_interp_F(c1,:,:)));
                                        tmp = ones(sz(1), sz(2));
                                        ccov_beta(f1,:) = tmp(:);
                                    else
                                        %do each session separately
                                        is1 = Z.group_session_to_average;
                                        %for is1=1:length(big_TOPO{f1}.v{v1}.s)
                                        tmp = squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_interp_F(c1,:,:));
                                        cbeta(f1,:) = tmp(:);
                                        sz = size(squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_interp_F(c1,:,:)));
                                        tmp = ones(sz(1), sz(2));
                                        ccov_beta(f1,:) = tmp(:);
                                        %end
                                    end
                                catch
                                    disp(['No data for subject ' int2str(f1) ' contrast ' int2str(c1)  ' chromophore ' hb ' and view ' int2str(v1)]);
                                end
                            end
                        end
                        if isfield(big_TOPO{1}.rendered_MNI{v1},'view_mask_2d')
                            for i0=1:length(big_TOPO)
                                cbeta(i0,:) = cbeta(i0,:).*big_TOPO{i0}.rendered_MNI{v1}.view_mask_2d(:)';
                            end
                        end
                        try
                            A = liom_anova(cbeta,ccov_beta,s1,s2,ns,min_s,level_subj);
                            %A = calc_hfgg(cbeta,A,ns,level_subj,1);
                        catch exception2
                            disp(exception2.identifier);
                            disp(exception2.stack(1));
                        end
                        
                        TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.A = A;
                        TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.c = xCon(c1);
                        
                        filestr = [num2str(Z.p_value) '_' spec_hemi '_' hb];
                        filestr_fig = [num2str(Z.p_value) ' ' spec_hemi ' ' hb];
                        info1 = [filestr xCon{Z.group_session_to_average}(c1).name];
                        info_for_fig1 = [filestr_fig xCon{Z.group_session_to_average}(c1).name];
                        F.contrast_info = info1;
                        F.contrast_info_for_fig = info_for_fig1;
                        F.contrast_info_both = [filestr xCon{Z.group_session_to_average}(c1).name]; %same for Pos and Neg, used for combined figures
                        F.contrast_info_both_for_fig = [filestr_fig xCon{Z.group_session_to_average}(c1).name]; %same for Pos and Neg, used for combined figures
                        
                        F.s_map = A.F;
                        
                        F.erdf = A.df;
                        F.eidf = A.dfbetween;
                        F.tstr = 'F'; %tstr;
                        F.hb = hb;
                        TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.F.erdf = F.erdf;
                        TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.F.eidf = F.eidf;
                        try
                            if Z.output_unc
                                DF = nirs_draw_figure(9,F,W,Z,[]);
                                H = nirs_copy_figure(H,DF,CF,c1,hb,1,F.tstr,0,0);
                            end
                            
                            if Z.LKC
                                DF = nirs_draw_figure(8,F,W,Z,A.LKC);
                                H = nirs_copy_figure(H,DF,CF,c1,hb,1,F.tstr,1,0);
                            end
                            
                            %A = calc_hfgg(cbeta,A,ns,B,nfac);
                            %For post-hoc contrasts only -- not done yet
                            switch Z.CorrectionMethod
                                case 1 %Huynh-Feldt
                                    F.erdf = F.erdf;
                                    F.eidf = F.eidf;
                                case 2 %Bonferroni
                                    p_value = Z.p_value/length(level_subj-1);
                                case 3 %Greenhouse-Gasser
                                    F.erdf = F.erdf;
                                    F.eidf = F.eidf;
                            end
                            %Calculate Huynh-Feldt and Greenhouse-Gasser corrections
                            %Where: at the site of minimal or maximal
                            %activation? or even at the highest of the two in absolute value
                            
                            %[EpsHF EpsList EpsGG]=GenCalcHFEps(Y,BTFacs,WInFacs,S)
                            
                        catch exception2
                            disp(exception2.identifier);
                            disp(exception2.stack(1));
                        end
                    catch exception
                        disp(exception.identifier);
                        disp(exception.stack(1));
                        disp(['Problem with a specific contrast ' int2str(c1) ' and chromophore ' hb ' for view ' int2str(v1)]);
                    end
                end
            end
            %save assembled figures
            call_save_assembled_figures(Z,W,H,0);
        end %if view_estimated
    end %end for v1
    save(ftopo,'TOPO');
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Could not do 1-way anova analysis');
end
out.NIRSmat = job.NIRSmat;
