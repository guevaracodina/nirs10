function out = nirs_run_liom_2way_anova(job)
%try to get factorial design specification
try
    %this is not coded up at all
    out = nirs_spm_run_factorial_design(job);
end
%Get anova info
% for l1=1:2
%     nlevel(l1) = job.fact(l1).levels;
%     names_level{l1} = job.fact(l1).name;
% end
anova2_sessions = job.anova2_sessions;
anova2_contrasts = job.anova2_contrasts;
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

% try
%     group_session_to_average = job.group_session_to_average;
% catch
%     group_session_to_average = 1;
% end
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
    anova_dir_name = 'AnovaTwoWay';
end
try 
    includeSubjectEffects = job.includeSubjectEffects;
catch
    includeSubjectEffects = 1;
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
Z.save_nifti_contrasts = 0;
Z.includeSubjectEffects = includeSubjectEffects;
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

%Load NIRS.mat information
try
    
    [dir0,dummy,dummy2] = fileparts(job.NIRSmat{1});
    %extract previous directory
    tmp = strfind(dir0,filesep);
    dir_root = dir0(1:tmp(end-3));
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
%             cbeta = zeros(ns,s1*s2);
%             ccov_beta = zeros(ns,s1*s2);
%             tmp = zeros(s1,s2);
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
            nC0 = length(anova2_contrasts);
            nS0 = length(anova2_sessions);
            
            %Loop over chromophores
            for h1=1:3 %including HbT
                hb = get_chromophore(h1);
                %Fill cbeta with session by contrast information
                sC = 0; %session counter
                Ns = length(big_TOPO{1}.v{v1}.s); %number of sessions
                for s1=1:Ns
                    %only selected sessions
                    if any(s1==anova2_sessions)
                        sC = sC + 1;
                        cC = 0; %contrast counter
                        for c1=1:nC
                            %only selected contrasts
                            if any(c1==anova2_contrasts)
                                cC = cC + 1;
                                %add to design matrix
                                %Skip F contrasts for now
                                if xCon(c1).STAT == 'T'
                                    fC = 0;
                                    %fill in cbeta
                                    for f1=1:ns
                                        try
                                            fC = fC+1;
                                            tmp = squeeze(big_TOPO{f1}.v{v1}.s{s1}.hb{h1}.c_interp_beta(c1,:,:));
                                            %                                             %assign space - all cells
                                            %                                             if cC == 1 && f1 == 1 && sC == 1
                                            %                                                 nC0 = length(anova2_contrasts);
                                            %                                                 nS0 = length(anova2_sessions);
                                            %                                                 cbeta{nS0,nC0} = zeros(ns,length(tmp(:)));
                                            %                                             end
                                            %                                             %further assign space - each cell
                                            %                                             if f1 == 0
                                            %                                                 cbeta{sC,cC} = zeros(ns,length(tmp(:)));
                                            %                                             end
                                            %now fill cbeta
                                            cbeta(fC,sC,cC,:) = tmp(:);
                                        catch exception %exception for difficulty filling cbeta
                                            disp(exception.identifier);
                                            disp(exception.stack(1));
                                            fC = fC-1; %remove
                                            %cbeta{sC,cC}(f1,:) = NaN(length(tmp(:)),1);
                                            disp(['No data for subject ' int2str(f1) ', contrast ' int2str(c1) ', session ' int2str(s1) ', chromophore ' hb ', and view ' int2str(v1)]);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end %end for s1
                
                X0 = []; %Design matrix of reduced model                
                ns0 = size(cbeta,1);
                X = zeros(ns0*nS0*nC0,nS0*nC0); %Design matrix of 2-way anova (full model)
                sC = 0; %session counter
                for s1=1:Ns
                    %only selected sessions
                    if any(s1==anova2_sessions)
                        sC = sC + 1;
                        cC = 0; %contrast counter
                        for c1=1:nC
                            %only selected contrasts
                            if any(c1==anova2_contrasts)
                                cC = cC + 1;
                                X0= [X0; eye(ns0)]; %intra-subject effects
                                rs = (1+ (cC-1)*ns0+nC0*(sC-1)*ns0);
                                re = rs+ns0-1;
                                X(rs:re,cC+nC0*(sC-1)) = ones(ns0,1); %treatment effects                                
                            end
                        end
                    end
                end
                if includeSubjectEffects
                    X = [X X0];
                else
                    X0 = [];
                end
           
                try
                    A = liom_group_2A(cbeta,X,X0,W.s1,W.s2,Z); %careful, s1 (session counter) not same as W.s1 (size of image)!
                catch exception2
                    disp(exception2.identifier);
                    disp(exception2.stack(1));
                end
                
                TOPO.v{v1}.group.hb{h1}.A = A;
                %TOPO.v{v1}.group.hb{h1}.A.c = xCon(c1);
                
                filestr = [num2str(p_value) '_' spec_hemi '_' hb];
                filestr_fig = [num2str(p_value) ' ' spec_hemi ' ' hb];
                info1 = [filestr ];
                info_for_fig1 = [filestr_fig ];
                F.contrast_info = info1;
                F.contrast_info_for_fig = info_for_fig1;
                F.contrast_info_both = [filestr ]; %same for Pos and Neg, used for combined figures
                F.contrast_info_both_for_fig = [filestr_fig ]; %same for Pos and Neg, used for combined figures
                
                F.T_map = A.tmap_group;
                F.erdf = A.erdf;
                F.eidf = A.eidf;
                F.tstr = 'F'; %tstr;
                F.hb = hb;
                try
                    DF = nirs_draw_figure(5,F,W,Z);
                catch exception2
                    disp(exception2.identifier);
                    disp(exception2.stack(1));
                end
                try
                    if GFIS, [Pu,Nu,Cu] = nirs_copy_figure(Pu,Nu,Cu,DF,CF,1,hb,1,F.tstr); end;
                catch exception2
                    disp(exception2.identifier);
                    disp(exception2.stack(1));
                end
                
            end %hb
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

