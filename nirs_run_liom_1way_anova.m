function out = nirs_run_liom_1way_anova(job)
%try to get factorial design specification
try
    %this is not coded up at all
    out = nirs_spm_run_factorial_design(job);
end
%Get anova info
anova_level = job.anova_level;
for l1=1:anova_level
    level_name{l1} = job.level(l1).level_name;
    level_subj{l1} = job.level(l1).level_subj;
end
Pu = [];
Nu = [];
Cu = [];
%Run simple group level analysis as a one sample t-test
p_value = job.contrast_p_value;
try
    GroupColorbars = job.GroupColorbars;
catch
    GroupColorbars = 0;
end
try
    SmallFigures = job.SmallFigures;
catch
    SmallFigures = 0;
end
write_neg_pos = 0;
try
    group_session_to_average = job.group_session_to_average;
catch
    group_session_to_average = 1;
end
%Booleans to choose which figures to write to disk, if any
switch job.contrast_figures
    case 0
        gen_fig = 0;
        gen_tiff = 0;
    case 1
        gen_fig = 1;
        gen_tiff = 1;
    case 2
        gen_fig = 1;
        gen_tiff = 0;
    case 3
        gen_fig = 0;
        gen_tiff = 1;
end
try
    cbar.c_min = job.override_colorbar.colorbar_override.colorbar_min;
    cbar.c_max = job.override_colorbar.colorbar_override.colorbar_max;
    cbar.c_min2 = job.override_colorbar.colorbar_override.colorbar_min2;
    cbar.c_max2 = job.override_colorbar.colorbar_override.colorbar_max2;
    cbar.colorbar_override = 1;
catch
    cbar.colorbar_override = 0;
end

output_unc = 1;
try
    switch job.figures_visible
        case 1
            cbar.visible = 'on';
        case 0
            cbar.visible = 'off';
    end
catch
    cbar.visible = 'off';
end
%Generate contrasts for inverted responses
GInv = 1;
GFIS = 1;
try
    anova_dir_name = job.anova_dir_name;
catch
    anova_dir_name = 'Anova';
end
%Structure for passing more generic data
Z = [];
Z.gen_fig = gen_fig;
Z.gen_tiff = gen_tiff;
Z.p_value = p_value;
Z.GroupColorbars = GroupColorbars;
Z.cbar = cbar;
Z.GInv = GInv;
Z.GFIS = GFIS;
Z.output_unc = output_unc;
Z.SmallFigures = SmallFigures;
Z.write_neg_pos = write_neg_pos;
min_s = 2;

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
    if Idx == 1 %assume same configuration for each subject - could be generalized
        fname_ch = NIRS.Dt.ana.rend;
        load(fname_ch);
    end
    ftopo = fullfile(dir1,'TOPO.mat');
    TOPO = [];
    load(ftopo);
    %large structure
    big_TOPO{Idx} = TOPO;
end
%create a new TOPO at the group level
TOPO = [];

%Load NIRS.mat information
try
    
    [dir0,dummy,dummy2] = fileparts(job.NIRSmat{1});
    %extract previous directory
    tmp = strfind(dir0,filesep);
    dir_root = dir0(1:tmp(end-2));
    dir_group = fullfile(dir_root,anova_dir_name);
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
            if GFIS
                Pu = figure('Visible',cbar.visible,'Name',['Group' '_' num2str(p_value) '_Pos'],'NumberTitle','off');
                Nu = figure('Visible',cbar.visible,'Name',['Group' '_' num2str(p_value) '_Neg'],'NumberTitle','off');
                Cu = figure('Visible',cbar.visible,'Name',['Group' '_' num2str(p_value)],'NumberTitle','off');
            end
            load Split
            F.split = split;
            F.pathn = Z.dir1;
            %CF: copy figure structure
            CF.GInv = GInv;
            CF.split = split;
            CF.nC = nC;
            
            %Loop over chromophores
            for h1=1:3 %including HbT
                hb = get_chromophore(h1);
                for c1=1:nC
                    try
                        %Skip F contrasts for now
                        if xCon(c1).STAT == 'T'
                            %fill in cbeta and ccov_beta
                            for f1=1:ns
                                try
                                    if ~isfield(big_TOPO{f1}.v{v1},'s')
                                        %group analysis of a group of
                                        %sessions analysis
                                        
                                        tmp = squeeze(big_TOPO{f1}.v{v1}.g.hb{h1}.c_interp_beta(c1,:,:));
                                        cbeta(f1,:) = tmp(:);
                                        tmp = squeeze(big_TOPO{f1}.v{v1}.g.hb{h1}.c_cov_interp_beta(c1,:,:));
                                        ccov_beta(f1,:) = tmp(:);
                                        
                                    else
                                        %for is1=1:length(big_TOPO{f1}.v{v1}.s)
                                        is1 = group_session_to_average;
                                        %do each session separately
                                        tmp = squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_interp_beta(c1,:,:));
                                        cbeta(f1,:) = tmp(:);
                                        tmp = squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_cov_interp_beta(c1,:,:));
                                        ccov_beta(f1,:) = tmp(:);
                                        %end
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
                                        is1 = group_session_to_average;
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
                        try
                            A = liom_anova(cbeta,ccov_beta,s1,s2,ns,min_s,level_subj);
                        catch exception2
                            disp(exception2.identifier);
                            disp(exception2.stack(1));
                        end
                        
                        TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.A = A;
                        TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.c = xCon(c1);
                        
                        filestr = [num2str(p_value) '_' spec_hemi '_' hb];
                        filestr_fig = [num2str(p_value) ' ' spec_hemi ' ' hb];
                        info1 = [filestr xCon(c1).name];
                        info_for_fig1 = [filestr_fig xCon(c1).name];
                        F.contrast_info = info1;
                        F.contrast_info_for_fig = info_for_fig1;
                        F.contrast_info_both = [filestr xCon(c1).name]; %same for Pos and Neg, used for combined figures
                        F.contrast_info_both_for_fig = [filestr_fig xCon(c1).name]; %same for Pos and Neg, used for combined figures
                        
                        F.T_map = A.F;
                        F.erdf = A.df;
                        F.eidf = A.dfbetween;
                        F.tstr = 'F'; %tstr;
                        F.hb = hb;
                        try
                            DF = nirs_draw_figure(4,F,W,Z);
                        catch exception2
                            disp(exception2.identifier);
                            disp(exception2.stack(1));
                        end
                        try
                            if GFIS, [Pu,Nu,Cu] = nirs_copy_figure(Pu,Nu,Cu,DF,CF,c1,hb,1,F.tstr); end;
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
            save_assembled_figures(Z,W,Cu,'','unc',0);
            try close(Pu); end
            try close(Nu); end
            try close(Cu); end
        end %if view_estimated
    end %end for v1
    save(ftopo,'TOPO');
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp(['Could not do anova analysis']);
end
out.NIRSmat = job.NIRSmat;
end

function hb = get_chromophore(h1)
switch h1
    case 1
        hb = 'HbO';
    case 2
        hb = 'HbR';
    case 3
        hb = 'HbT';
end
end

