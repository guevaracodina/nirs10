function out = nirs_run_liom_group(job)
%try to get factorial design specification
try
    %this is not coded up at all
    out = nirs_spm_run_factorial_design(job);
end
Pu = [];
Nu = [];
Cu = [];
%Run simple group level analysis as a one sample t-test
FFX = job.FFX_or_RFX;
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
try
    write_neg_pos = job.write_neg_pos;
catch
    write_neg_pos = 0;
end
try
    group_session_to_average = job.group_session_to_average;
catch
    group_session_to_average = 1;
end
try
    group_dir_name = job.group_dir_name;
catch
    group_dir_name = 'Group';
end
try
    save_nifti_contrasts = job.save_nifti_contrasts;
catch
    save_nifti_contrasts = 0;
end
try
    simple_sum = job.simple_sum;
catch
    simple_sum = 1;
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
try %not used currently as false positives correction not implemented yet
    output_unc = 1; %job.output_unc;
catch
    output_unc = 1;
end
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
try
    GInv = job.GenerateInverted;
catch
    GInv = 1;
end
try
    GFIS = job.GroupFiguresIntoSubplots;
catch
    GFIS = 1;
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
Z.save_nifti_contrasts = save_nifti_contrasts;
Z.min_s = 2;
Z.FFX = FFX;
nS = size(job.NIRSmat,1);
Z.nS = nS;
Z.simple_sum = simple_sum;
Z.group_session_to_average = group_session_to_average;
if FFX || nS==1
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
        end
        TOPO = [];
        load(ftopo);
        %large structure
        big_TOPO{Idx} = TOPO;
    end
    %create a new TOPO at the group level
    TOPO = [];
end


%Loop over all subjects for FFX but go through only once for RFX
for Idx=1:nl
    %Load NIRS.mat information
    try
        if FFX || nS==1
            NIRS = [];
            load(job.NIRSmat{Idx,1});
            %load topographic information (formerly known as preproc_info)
            fname_ch = NIRS.Dt.ana.rend;
            load(fname_ch);
            %load(fullfile(dir1,'SPM.mat'));
            ftopo = NIRS.TOPO; %fullfile(dir1,'TOPO.mat');
            [dir1 fil1 ext1] = fileparts(ftopo);
            TOPO = [];
            load(ftopo);
            Z.dir1 = dir1;
        else
            [dir0,dummy,dummy2] = fileparts(job.NIRSmat{1});
            %extract previous directory
            tmp = strfind(dir0,filesep);
            dir_root = dir0(1:tmp(end-3));
            dir_group = fullfile(dir_root, group_dir_name);
            if ~exist(dir_group,'dir'), mkdir(dir_group); end
            %store in same directory as first subject
            ftopo = fullfile(dir_group,'TOPO.mat');
            %save a NIRS structure for the group
            newNIRSlocation = fullfile(dir_group,'NIRS.mat');
            NIRS.TOPO = ftopo;
            save(newNIRSlocation,'NIRS');
            job.NIRSmat{nl,1} = newNIRSlocation;
            Z.dir1 = dir_group;
         end
        
        %Big loop over views
        for v1=1:6
            view_estimated = 0;
            try
                if FFX || nS==1
                    if isfield(TOPO.v{v1},'s1')
                        s1 = TOPO.v{v1}.s1;
                        s2 = TOPO.v{v1}.s2;
                        ns = length(TOPO.v{v1}.s); %number of sessions
                        view_estimated = 1;
                    end
                else
                    if isfield(big_TOPO{1}.v{v1},'s1')
                        s1 = big_TOPO{1}.v{v1}.s1;
                        s2 = big_TOPO{1}.v{v1}.s2;
                        ns = nS; %number of subjects
                        view_estimated = 1;
                    end
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
                clear xCon
                if FFX || nS==1
                    %Contrasts
                    %xCon = TOPO.SSxCon; %need to generalize
                    xConS = TOPO.SSxConS;
                    TOPO.v{v1}.g.ns = ns;
                    TOPO.v{v1}.g.min_s = Z.min_s;
                    TOPO.v{v1}.g.s1 = s1;
                    TOPO.v{v1}.g.s2 = s2;
                else
                    TOPO.v{v1}.group.ns = ns;
                    TOPO.v{v1}.group.min_s = Z.min_s;
                    TOPO.v{v1}.group.s1 = s1;
                    TOPO.v{v1}.group.s2 = s2;
                    %Contrasts -- assume same contrasts for all subjects
                    try
                        xCon = big_TOPO{1}.SSxCon; % big_TOPO{1}.xCon}; %????
                        if isempty(xCon) %for contrasts defined over several sessions
                            xCon = big_TOPO{1}.xCon;
                        end
                    catch
                        xCon = big_TOPO{1}.xCon;
                    end
                    %Contrasts - assume same contrasts for all subjects
                    %xCon = big_TOPO{1}.xCon;
                end
                
                tmp = zeros(s1,s2);
                if exist('xCon','var')
                    nC = length(xCon);
                    cbeta = zeros(ns,s1*s2);
                    ccov_beta = zeros(ns,s1*s2);
                    Sess = [];
                    Cp = [];
                else
                    %build up list of contrasts and of their names
                    Clist = xConS{1};
                    Nlist = {}; for t2=1:length(xConS{1}), Nlist = [Nlist; xConS{1}(t2).name]; end
                    for t1=2:ns
                        for t2=1:length(xConS{t1})
                            if ~any(strcmp(xConS{t1}(t2).name,Nlist))
                                %add to the list
                                Clist(end+1) = xConS{t1}(t2);
                                Nlist{end+1} = xConS{t1}(t2).name;
                            end
                        end
                    end
                    nC = length(Clist);
                    xCon = Clist;
                    %for each contrast, build list of available sessions,
                    %and position of these contrasts in SSxCon
                    for c1=1:nC
                        Sess{c1} = [];
                        %Cp{c1} = []; %contrast position in the SSxCon list of that session
                        for t1=1:ns
                            Nlist = {}; for t2=1:length(xConS{t1}), Nlist = [Nlist; xConS{t1}(t2).name]; end
                            if any(strcmp(xCon(c1).name, Nlist)),
                                Sess{c1} = [Sess{c1} t1];
                                Cp{c1,t1} = find(strcmp(xCon(c1).name, Nlist)==1);
                            end
                        end
                    end
                end
                if GFIS
                    Pu = figure('Visible',cbar.visible,'Name',['Group' '_' num2str(p_value) '_Pos'],'NumberTitle','off');
                    %subplot(fh0P,nC,3,1);
                    if GInv
                        Nu = figure('Visible',cbar.visible,'Name',['Group' '_' num2str(p_value) '_Neg'],'NumberTitle','off');
                        Cu = figure('Visible',cbar.visible,'Name',['Group' '_' num2str(p_value)],'NumberTitle','off');
                    end
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
                    for c1=1:nC %Loop over contrasts
                        %Positive stats
                        [Pu,Nu,Cu,TOPO,big_TOPO] = fill_group(Pu,Nu,Cu,TOPO,big_TOPO,v1,c1,h1,Z,W,F,CF,xCon,ns,Sess,Cp,1);
                        %Negative stats
                        [Pu,Nu,Cu,TOPO,big_TOPO] = fill_group(Pu,Nu,Cu,TOPO,big_TOPO,v1,c1,h1,Z,W,F,CF,xCon,ns,Sess,Cp,0);
                    end
                end
                %save assembled figures
                if GFIS
                    if Z.write_neg_pos || ~GInv
                        save_assembled_figures(Z,W,Pu,'Pos','unc',0);
                    else
                        try close(Pu); end
                    end
                    if GInv
                        if Z.write_neg_pos
                            save_assembled_figures(Z,W,Nu,'Neg','unc',0);
                        else
                            try close(Nu); end
                        end
                        save_assembled_figures(Z,W,Cu,'','unc',0);
                    else
                        try close(Nu); end
                        try close(Cu); end
                    end
                end
            end %if view_estimated
        end %end for v1
        TOPO.xCon = xCon; %would not work if new contrasts are later added
        if exist('Sess','var')
            TOPO.Sess = Sess;
            TOPO.Cp = Cp;
        end
        save(ftopo,'TOPO');
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not do group analysis']);
    end
    %NIRS.TOPO = ftopo;
    %save(job.NIRSmat{Idx,1});
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
                                              
function [Pu,Nu,Cu,TOPO,big_TOPO] = fill_group(Pu,Nu,Cu,TOPO,big_TOPO,v1,c1,h1,Z,W,F,CF,xCon,ns,Sess,Cp,shb)
try
    hb = get_chromophore(h1);
    if shb
        strA = 'Pos';
        strB = 'Positive';
        sign_hb = 1;
    else
        strA = 'Neg';
        strB = 'Negative';
        sign_hb = -1;
    end
    if xCon(c1).STAT == 'T'
        fc = 0; %used only for FFX || nS==1
        %fill in cbeta and ccov_beta
        for f1=1:ns
            if Z.FFX || Z.nS==1
                %select sessions which had the contrast
                if any(f1==Sess{c1})
                    fc = fc+1;
                    %now use Cp{c1,f1} to access the required c_interp_beta instead of c1
                    tmp = sign_hb*squeeze(TOPO.v{v1}.s{f1}.hb{h1}.c_interp_beta(Cp{c1,f1},:,:));
                    cbeta(fc,:) = tmp(:);
                    tmp = squeeze(TOPO.v{v1}.s{f1}.hb{h1}.c_cov_interp_beta(Cp{c1,f1},:,:));
                    ccov_beta(fc,:) = tmp(:);
                end
            else
                if ~isfield(big_TOPO{f1}.v{v1},'s')
                    %group analysis of a group of sessions analysis
                    tmp = sign_hb*squeeze(big_TOPO{f1}.v{v1}.g.hb{h1}.c_interp_beta(c1,:,:));
                    cbeta(f1,:) = tmp(:);
                    tmp = squeeze(big_TOPO{f1}.v{v1}.g.hb{h1}.c_cov_interp_beta(c1,:,:));
                    ccov_beta(f1,:) = tmp(:);
                else
                    %for is1=1:length(big_TOPO{f1}.v{v1}.s)
                    is1 = Z.group_session_to_average;
                    %do each session separately
                    tmp = sign_hb*squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_interp_beta(c1,:,:));
                    cbeta(f1,:) = tmp(:);
                    tmp = squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_cov_interp_beta(c1,:,:));
                    ccov_beta(f1,:) = tmp(:);
                    %end
                end
            end
        end
    else %quick fix for F stats
        fc = 0; %used only for FFX || nS==1
        for f1=1:ns
            if Z.FFX || Z.nS==1
                %select sessions which had the contrast
                if any(f1==Sess{c1})
                    fc = fc+1;
                    tmp = squeeze(TOPO.v{v1}.s{f1}.hb{h1}.c_interp_F(Cp{c1,f1},:,:));
                    cbeta(fc,:) = tmp(:);
                end
            else
                if ~isfield(big_TOPO{f1}.v{v1},'s')
                    %group analysis of a group of sessions analysis
                    tmp = squeeze(big_TOPO{f1}.v{v1}.g.hb{h1}.c_interp_F(c1,:,:));
                    cbeta(f1,:) = tmp(:);
                else
                    %do each session separately
                    is1 = Z.group_session_to_average;
                    %for is1=1:length(big_TOPO{f1}.v{v1}.s)
                    tmp = squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_interp_F(c1,:,:));
                    cbeta(f1,:) = tmp(:);
                    %end
                end
            end
        end
    end
    %Generate group result as t-stat
    try
        if strcmp(xCon(c1).STAT,'T')
            G = liom_group(cbeta,ccov_beta,W.s1,W.s2,Z.min_s,Z.FFX,Z.simple_sum);
        else
            G = liom_group_F(cbeta,W.s1,W.s2);
        end
    catch exception2
        disp(exception2.identifier);
        disp(exception2.stack(1));
    end
    if Z.FFX || Z.nS==1
        TOPO.v{v1}.g.hb{h1}.c{2*c1-shb}.Tmap = G.tmap_group;
        TOPO.v{v1}.g.hb{h1}.c{2*c1-shb}.erdf = G.erdf_group;
        TOPO.v{v1}.g.hb{h1}.c{2*c1-shb}.beta_group = G.beta_group;
        TOPO.v{v1}.g.hb{h1}.c{2*c1-shb}.std_group = G.std_group;
        TOPO.v{v1}.g.hb{h1}.c{2*c1-shb}.type = strB;
        TOPO.v{v1}.g.hb{h1}.c{2*c1-shb}.var_bs = G.var_bs;
        TOPO.v{v1}.g.hb{h1}.c{2*c1-shb}.c = xCon(c1);
    else
        TOPO.v{v1}.group.hb{h1}.c{2*c1-shb}.Tmap = G.tmap_group;
        TOPO.v{v1}.group.hb{h1}.c{2*c1-shb}.erdf = G.erdf_group;
        TOPO.v{v1}.group.hb{h1}.c{2*c1-shb}.beta_group = G.beta_group;
        TOPO.v{v1}.group.hb{h1}.c{2*c1-shb}.std_group = G.std_group;
        TOPO.v{v1}.group.hb{h1}.c{2*c1-shb}.type = strB;
        TOPO.v{v1}.group.hb{h1}.c{2*c1-shb}.var_bs = G.var_bs;
        TOPO.v{v1}.group.hb{h1}.c{2*c1-shb}.c = xCon(c1);
    end
    erdf_group = max(G.erdf_group(:)); %quick fix...
    filestr = [num2str(Z.p_value) '_' W.spec_hemi '_' hb];
    filestr_fig = [num2str(Z.p_value) ' ' W.spec_hemi ' ' hb];
    info1 = [filestr '_' strA xCon(c1).name];
    info_for_fig1 = [filestr_fig ' ' strA xCon(c1).name];
    F.contrast_info = info1;
    F.contrast_info_for_fig = info_for_fig1;
    F.contrast_info_both = [filestr xCon(c1).name]; %same for Pos and Neg, used for combined figures
    F.contrast_info_both_for_fig = [filestr_fig xCon(c1).name]; %same for Pos and Neg, used for combined figures
    
    F.T_map = G.tmap_group;
    F.erdf = erdf_group;
    F.eidf = xCon(c1).eidf;
    F.tstr = 'T'; %xCon(c1).STAT; %tstr;
    F.hb = hb;
    if strcmp(F.tstr,'T')
        F.con = G.beta_group;
        F.ess = [];
    else
        F.con = [];
        F.ess = G.beta_group;
    end
    DF = [];
    try
        DF = nirs_draw_figure(4,F,W,Z);
    catch exception2
        disp(exception2.identifier);
        disp(exception2.stack(1));
    end
    try
        if Z.GFIS, [Pu,Nu,Cu] = nirs_copy_figure(Pu,Nu,Cu,DF,CF,c1,hb,shb,F.tstr); end;
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