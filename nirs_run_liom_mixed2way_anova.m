function out = nirs_run_liom_mixed2way_anova(job)
%try to get factorial design specification
try
    %this is not coded up at all
    out = nirs_spm_run_factorial_design(job);
end
Z = get_contrast_group_common_options(job);
%Get between factor
Z.anova_level = job.anova_level;
for l1=1:Z.anova_level
    level_name{l1} = job.level(l1).level_name;
    level_subj{l1} = job.level(l1).level_subj;
end
Z.anova2_sessions = unique(job.anova2_sessions); %might be needed for Sarah
Z.level_subj = level_subj;
Z.anova2_contrasts = unique(job.anova2_contrasts);
Z.LKC = job.StatMethod;
Z.StatStr = 'EC';
Z.CorrectionMethod = job.CorrectionMethod;
Z.GInv = 1; %enforce
Z.GFIS = 1; %enforce
Z.anova_dir_name = job.anova_dir_name;
Z.includeSubjectEffects = job.includeSubjectEffects;
Z.AvgInterpBetaMode = 0; %A boolean to invoke the mode where whole sessions are averaged and compared
number_dir_to_remove = job.number_dir_to_remove;
min_s = 2; %not used
p_value0 = Z.p_value;
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
    [NIRS newNIRSlocation]= nirs_load(newNIRSlocation,job.NIRSmatCopyChoice,job.force_redo);
    NIRS.TOPO = ftopo;
    NIRSgroup = NIRS;
    save(newNIRSlocation,'NIRS');
    job.NIRSmat{nl,1} = newNIRSlocation;
    Z.dir1 = dir_group;
    if ~isfield(NIRS.flags,'anovamixed2_OK') || job.force_redo
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
                if isfield(rendered_MNI0{v1},'view_mask_2d')
                    W.brain_view_mask_2d = rendered_MNI0{v1}.view_mask_2d;
                end
                W = nirs_get_common_brain_mask(W,big_TOPO,v1);
                W = nirs_get_boundary(W,job);
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
                try
                    nC = length(xCon{1});
                catch
                     nC = length(xCon);
                     xCon0 = xCon; 
                     clear xCon;
                     
                     xCon{1} = xCon0;
                end
                %Add loop over effects to look at:
                %1: interaction A*B,
                %2: main A,
                %3: main B
                %4: effect of A at each level of B
                %5: effect of B at each level of A
                clear cbeta Eps
                for z1=1:5
                    %Handles for assembled figures
                    clear H
                    nC0 = length(Z.anova2_contrasts);
                    nS0 = length(Z.anova2_sessions);
                    if nC0 > 1
                        if nS0 >1
                            disp(['Problem, there are too many sessions or ' ...
                                'too many conditions: either there should be only ' ...
                                'one session, or only one condition, with '...
                                'the other being the 2nd anova factor']);
                        else
                            Beffect = 'Conditions';
                            nB = nC0;
                        end
                    else
                        if nS0 > 1
                            Beffect = 'Sessions';
                            nB = nS0;
                        else
                            disp(['Only one session and one condition: ' ...
                                'cannot do a mixed-2-anova, missing the second factor!']);
                        end
                    end
                    
                    switch z1
                        case {1,2,3}
                            H = initialize_assembled_figure_handles;
                            H = initialize_assembled_figures(Z,H,0,'Group');
                        case 4
                            for y1 = 1:Z.anova_level
                                H{y1} = initialize_assembled_figure_handles;
                                H{y1} = initialize_assembled_figures(Z,H{y1},0,'Group');
                            end
                        case 5
                            for y1 = 1:nB
                                H{y1} = initialize_assembled_figure_handles;
                                H{y1} = initialize_assembled_figures(Z,H{y1},0,'Group');
                            end
                    end
                    
                    load Split
                    F.split = split;
                    F.pathn = Z.dir1;
                    %CF: copy figure structure
                    CF.GInv = Z.GInv;
                    CF.split = split;
                    CF.nC = nC;
                    
                    Z.p_value = job.contrast_p_value; %reset because of redefinition of p_value later for post-hoc contrasts
                    %Loop over chromophores
                    
                    for h1=1:3 %including HbT
                        hb = get_chromophore(h1);
                        %Fill cbeta with session by contrast information
                        sC = 0; %session counter
                        if isfield(big_TOPO{1}.v{v1},'s') 
                            Ns = length(big_TOPO{1}.v{v1}.s); %number of sessions
                            group_sessions = 0;
                        else
                            if isfield(big_TOPO{1}.v{v1},'g') 
                                Ns = 1;
                                group_sessions = 1;
                            else
                                disp('Corrupted data, 2-anova will break');
                            end
                        end
                        for s1=1:Ns
                            %only selected sessions
                            if any(s1==Z.anova2_sessions)
                                sC = sC + 1;
                                cC = 0; %contrast counter
                                for c1=1:nC
                                    %only selected contrasts
                                    if any(c1==Z.anova2_contrasts)
                                        cC = cC + 1;
                                        %add to design matrix
                                        %Skip F contrasts for now
                                        if xCon{s1}(c1).STAT == 'T'
                                            fC = 0;
                                            %fill in cbeta
                                            if z1==1
                                                for f1=1:ns
                                                    try
                                                        fC = fC+1;
                                                        if group_sessions
                                                            %if isfield(big_TOPO{f1}.v{v1}.group.hb{h1},'beta_map')
                                                                %tmp = squeeze(big_TOPO{f1}.v{v1}.group.hb{h1}.c{2*c1-1}.beta_group);
                                                                tmp = squeeze(big_TOPO{f1}.v{v1}.g{1}.hb{h1}.beta_map(c1,:,:));
                                                            %else
                                                            %    tmp = squeeze(big_TOPO{f1}.v{v1}.group.hb{h1}.c_interp_beta(c1,:,:));
                                                            %end
                                                        else
                                                            if isfield(big_TOPO{f1}.v{v1}.s{s1}.hb{h1},'beta_map')
                                                                tmp = squeeze(big_TOPO{f1}.v{v1}.s{s1}.hb{h1}.beta_map(c1,:,:));
                                                            else
                                                                tmp = squeeze(big_TOPO{f1}.v{v1}.s{s1}.hb{h1}.c_interp_beta(c1,:,:));
                                                            end
                                                        end
                                                        %now fill cbeta
                                                        cbeta{h1}(fC,sC,cC,:) = tmp(:);
                                                        if isfield(big_TOPO{f1}.rendered_MNI{v1},'view_mask_2d')
                                                            cbeta{h1}(fC,sC,cC,:) = squeeze(cbeta{h1}(fC,sC,cC,:))'.*big_TOPO{f1}.rendered_MNI{v1}.view_mask_2d(:)';
                                                        end
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
                            end
                        end %end for s1
                        
                        if z1==1 && h1 == 1%no need to repeat calculation
                            cbeta{h1} = squeeze(cbeta{h1}); %squeeze not needed
                            %Construct each part of design matrix separately
                            Xs = []; %subject effects
                            Xa = []; %main effect of A (Subjects)
                            Xb = []; %main effect of B (could be condition or session, depending on the protocol)
                            Xab = []; %main effect of A*B interaction
                            ns0 = size(cbeta{h1},1);
                            
                            sX = ns0*nS0*nC0;                           
                            %subject effects
                            if Z.includeSubjectEffects
                                for s1=1:nS0
                                    for c1=1:nC0
                                        Xs= [Xs; eye(ns0)]; %intra-subject effects
                                    end
                                end
                            end
                            if nC0 > 1
                                nC0X = nC0;
                            else
                                %sessions rather than contrasts
                                nC0X = nS0;
                            end
                            %main effect of subject
                            for s1=2:length(level_subj)
                                x = zeros(sX,1);
                                for c1=1:nC0X
                                    x(level_subj{s1}+((c1-1)*ns0)) = 1;
                                    x(level_subj{s1-1}+((c1-1)*ns0)) = -1;
                                end
                                Xa = [Xa x];
                            end
                            
                            %main effect of condition 
                            for c1=1:nC0
                                if c1 > 1
                                    x2 = x; %previous x
                                end
                                x = zeros(sX,1);
                                for s1=1:nS0
                                    rs = (1+ (c1-1)*ns0+nC0*(s1-1)*ns0);
                                    re = rs+ns0-1;
                                    x(rs:re) = ones(ns0,1);
                                end
                                if c1 > 1
                                    Xb = [Xb -x+x2];
                                end
                            end
                            
                            %main effect of session 
                            for S1=1:nS0
                                if S1 > 1
                                    x2 = x; %previous x
                                end
                                x = zeros(sX,1);
                                rs = (1+ (S1-1)*ns0);
                                re = rs+ns0-1;
                                x(rs:re) = ones(ns0,1);
                                if S1 > 1
                                    Xb = [Xb -x+x2];
                                end
                            end
                            
                            %interaction of factor 1 and intensity (contrast)
                            for s1=2:length(level_subj)
                                for c1=2:nC0
                                    x = zeros(sX,1);
                                    x(level_subj{s1}+(c1-2)*ns0) = 1;
                                    x(level_subj{s1-1}+(c1-1)*ns0) = 1;
                                    x(level_subj{s1}+(c1-1)*ns0) = -1;
                                    x(level_subj{s1-1}+(c1-2)*ns0) = -1;
                                    Xab = [Xab x];
                                end
                            end
                            
                            %interaction of factor 1 and intensity (session)
                            for s1=2:length(level_subj)
                                for S1=2:nS0
                                    x = zeros(sX,1);
                                    x(level_subj{s1}+(S1-2)*ns0) = 1;
                                    x(level_subj{s1-1}+(S1-1)*ns0) = 1;
                                    x(level_subj{s1}+(S1-1)*ns0) = -1;
                                    x(level_subj{s1-1}+(S1-2)*ns0) = -1;
                                    Xab = [Xab x];
                                end
                            end
                            
                            %effect of A at each level of B
                            XeA = [];
                            if nS0 == 1
                                for c1=1:nC0 %B
                                    for s1=2:length(level_subj) %B
                                        x = zeros(sX,1);
                                        x(level_subj{s1}+(c1-1)*ns0) = 1;
                                        x(level_subj{s1-1}+(c1-1)*ns0) = -1;
                                        XeA = [XeA x];
                                    end
                                end
                            else
                                for S1=1:nS0 %B
                                    for s1=2:length(level_subj) %B
                                        x = zeros(sX,1);
                                        x(level_subj{s1}+(S1-1)*ns0) = 1;
                                        x(level_subj{s1-1}+(S1-1)*ns0) = -1;
                                        XeA = [XeA x];
                                    end
                                end
                            end
                            %effect of B at each level of A
                            XeB = [];
                            for s1=1:length(level_subj) %A
                                for c1=2:nC0 %B
                                    x = zeros(sX,1);
                                    x(level_subj{s1}+(c1-1)*ns0) = 1;
                                    x(level_subj{s1}+(c1-2)*ns0) = -1;
                                    XeB = [XeB x];
                                end
                            end
                            
                            for s1=1:length(level_subj) %A
                                for S1=2:nS0 %B
                                    x = zeros(sX,1);
                                    x(level_subj{s1}+(S1-1)*ns0) = 1;
                                    x(level_subj{s1}+(S1-2)*ns0) = -1;
                                    XeB = [XeB x];
                                end
                            end
                            %overall mean
                            M = ones(sX,1);
                        end
                        %Various tests:
                        switch z1
                            case 1
                                %Interaction A*B
                                strA = 'intAB';
                                X0 = [Xa Xb M Xs]; %should we remove Xs???
                                X = [Xab X0];
                                A = liom_group_2A(cbeta{h1},X,X0,W.s1,W.s2,Z); %careful, s1 (session counter) not same as W.s1 (size of image)!
                                %Mauchly test
                                %A.Mauchly = nirs_Mauchly(cbeta{h1},Z.p_value);
                                %fill B:
                                clear B
                                B.WInFacs = zeros(sX,1);
                                %                                 %1st factor: Subjects
                                %                                 for i0=1:nS0
                                %                                     B.WInFacs((i0-1)*ns0*nC0+(1:ns0*nC0),1) = i0;
                                %                                 end
                                %2nd factor: Contrasts or sessions
                                tmp0 = zeros(nC0*ns0,nS0);
                                for i0=1:nC0
                                    for S1=1:nS0
                                        tmp0((i0-1)*ns0+(1:ns0),S1) = i0*S1;
                                    end
                                end
                                B.WInFacs(:,1) = tmp0(:);
                                B.S = zeros(ns0*nB,1);
                                for is0 = 1:ns0
                                    for nB0 = 1:nB
                                        B.S(is0+ns0*(nB0-1)) = is0;
                                    end
                                end
                                %                                 for is0 = 1:ns0
                                %                                     for l1=1:Z.anova_level
                                %                                         if any(is0 == level_subj{l1})
                                %                                             for nB0 = 1:nB
                                %                                                 B.S(is0+ns0*(nB0-1)) = l1;
                                %                                             end
                                %                                         end
                                %                                     end
                                %                                 end
                                %
                                if nC0 > 1
                                    x = zeros(sX,1);
                                    for s1=1:length(level_subj)
                                        for c1=1:nC0
                                            x(level_subj{s1}+(c1-1)*ns0) = s1;
                                        end
                                    end
                                else
                                    x = zeros(sX,1);
                                    for s1=1:length(level_subj)
                                        for S1=1:nS0
                                            x(level_subj{s1}+(S1-1)*ns0) = s1;
                                        end
                                    end
                                end
                                B.BTFacs = x;
                                %The number of subjects at each level of
                                %the between-factor must be equal!
                                A = calc_hfgg(cbeta{h1},A,ns0,B,2);
                                
                                %Store for cases 2 to 5
                                Eps = A.Eps;
                                %A is output only to get strA -- not clean
                                [TOPO H A] = call_figure_2anova(TOPO,H,Z,W,F,A,CF,v1,h1,hb,strA,z1,A.LKC);
                            case 2
                                %Main effect of A (Subjects)
                                strA = 'mainA';
                                X0 = [Xb M]; %remove subject effects
                                X = [Xa X0];
                                A = liom_group_2A(cbeta{h1},X,X0,W.s1,W.s2,Z); %careful, s1 (session counter) not same as W.s1 (size of image)!
                                A.Eps = Eps;
                                [TOPO H A] = call_figure_2anova(TOPO,H,Z,W,F,A,CF,v1,h1,hb,strA,z1,A.LKC);
                            case 3
                                %Main effect of B (Conditions or perhaps Sessions)
                                strA = 'mainB';
                                X0 = [Xa M Xs];
                                X = [Xb X0];
                                A = liom_group_2A(cbeta{h1},X,X0,W.s1,W.s2,Z); %careful, s1 (session counter) not same as W.s1 (size of image)!
                                A.Eps = Eps;
                                [TOPO H A] = call_figure_2anova(TOPO,H,Z,W,F,A,CF,v1,h1,hb,strA,z1,A.LKC);
                            case 4
                                %Effect of A within levels of B
                                %loop over levels of B
                                for y1=1:Z.anova_level
                                    strA = ['effAonB' int2str(y1)];
                                    ind = 1:size(XeA,2);
                                    ind(y1) = [];
                                    X0 = [XeA(:,ind) Xb M]; %remove subject effects
                                    X  = [XeA Xb M];
                                    %need to adjust p value:
                                    Z.p_value = p_value0/Z.anova_level;
                                    A = liom_group_2A(cbeta{h1},X,X0,W.s1,W.s2,Z); %careful, s1 (session counter) not same as W.s1 (size of image)!
                                    A.Eps = Eps;
                                    A.y1 = y1;
                                    [TOPO H{y1} A] = call_figure_2anova(TOPO,H{y1},Z,W,F,A,CF,v1,h1,hb,strA,z1,A.LKC);
                                end
                            case 5
                                %Effect of B within levels of A
                                for y1=1:nB
                                    strA = ['effBonA' int2str(y1)];
                                    ind = 1:size(XeB,2);
                                    ind(y1) = [];
                                    X0 = [XeB(:,ind) Xa M Xs];
                                    X  = [XeB Xa M Xs];
                                    %need to adjust p value:
                                    Z.p_value = p_value0/nB;
                                    A = liom_group_2A(cbeta{h1},X,X0,W.s1,W.s2,Z); %careful, s1 (session counter) not same as W.s1 (size of image)!
                                    A.Eps = Eps;
                                    A.y1 = y1;
                                    [TOPO H{y1} A] = call_figure_2anova(TOPO,H{y1},Z,W,F,A,CF,v1,h1,hb,strA,z1,A.LKC);
                                end
                        end
                    end %hb
                    %save assembled figures
                    switch z1
                        case {1,2,3}
                            Z.strA = A.strA;
                            call_save_assembled_figures(Z,W,H,0);
                        case 4
                            for y1 = 1:nC0
                                Z.strA = ['effAonB' int2str(y1)];
                                call_save_assembled_figures(Z,W,H{y1},0);
                            end
                        case 5
                            for y1 = 1:nB
                                Z.strA = ['effBonA' int2str(y1)];
                                call_save_assembled_figures(Z,W,H{y1},0);
                            end
                    end
                end %for z1
            end %if view_estimated
        end %end for v1
        save(ftopo,'TOPO');
        save(fullfile(dir_group,'big_TOPO.mat'),'big_TOPO','-v7.3');
        NIRSgroup.flags.anovamixed2_OK = 1;
        NIRS = NIRSgroup;
        save(newNIRSlocation,'NIRS');
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Could not do mixed 2-anova analysis');
end
out.NIRSmat = job.NIRSmat;