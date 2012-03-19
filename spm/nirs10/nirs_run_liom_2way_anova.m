function out = nirs_run_liom_2way_anova(job)
%try to get factorial design specification
try
    %this is not coded up at all
    out = nirs_spm_run_factorial_design(job);
end
Z = get_contrast_group_common_options(job);
Z.anova2_sessions = job.anova2_sessions;
Z.anova2_contrasts = job.anova2_contrasts;
Z.LKC = job.StatMethod;
Z.StatStr = 'EC';
Z.CorrectionMethod = job.CorrectionMethod;
Z.GInv = 1; %enforce
Z.GFIS = 1; %enforce
Z.anova_dir_name = job.anova_dir_name;
Z.includeSubjectEffects = job.includeSubjectEffects;
number_dir_to_remove = job.number_dir_to_remove;
min_s = 2;
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
            nC = length(xCon);
            
            %Add loop over effects to look at:
            %1: interaction A*B,
            %2: main A,
            %3: main B
            %4: effect of A at each level of B
            %5: effect of B at each level of A
            for z1=1:5
                %Handles for assembled figures
                H = initialize_assembled_figure_handles;
                H = initialize_assembled_figures(Z,H,0,'Group');
                
                load Split
                F.split = split;
                F.pathn = Z.dir1;
                %CF: copy figure structure
                CF.GInv = Z.GInv;
                CF.split = split;
                CF.nC = nC;
                nC0 = length(Z.anova2_contrasts);
                nS0 = length(Z.anova2_sessions);
                Z.p_value = job.contrast_p_value; %reset because of redefinition of p_value later for post-hoc contrasts
                %Loop over chromophores
                for h1=1:3 %including HbT
                    hb = get_chromophore(h1);
                    %Fill cbeta with session by contrast information
                    sC = 0; %session counter
                    Ns = length(big_TOPO{1}.v{v1}.s); %number of sessions
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
                                                    if isfield(big_TOPO{f1}.v{v1}.s{s1}.hb{h1},'beta_map')
                                                        tmp = squeeze(big_TOPO{f1}.v{v1}.s{s1}.hb{h1}.beta_map(c1,:,:));
                                                    else
                                                        tmp = squeeze(big_TOPO{f1}.v{v1}.s{s1}.hb{h1}.c_interp_beta(c1,:,:));
                                                    end
                                                    %now fill cbeta
                                                    cbeta{h1}(fC,sC,cC,:) = tmp(:);
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
                        
                        %Construct each part of design matrix separately
                        Xs = []; %subject effects
                        Xa = []; %main effect of A (Sessions)
                        Xb = []; %main effect of B (Intensity -- contrast)
                        Xab = []; %main effect of A*B interaction
                        ns0 = size(cbeta{h1},1);
                        
                        sX = ns0*nS0*nC0;
                        %subject effects
                        for s1=1:nS0
                            for c1=1:nC0
                                Xs= [Xs; eye(ns0)]; %intra-subject effects
                            end
                        end
                        %main effect of session
                        for s1=1:nS0
                            if s1 > 1
                                x2 = x; %previous x
                            end
                            x = zeros(sX,1);
                            for c1=1:nC0
                                rs = (1+ (c1-1)*ns0+nC0*(s1-1)*ns0);
                                re = rs+ns0-1;
                                x(rs:re) = ones(ns0,1);
                            end
                            if s1 > 1
                                Xa = [Xa -x+x2];
                            end
                        end
                        
                        %main effect of intensity (contrast)
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
                        
                        %interaction of session and intensity (contrast)
                        for s1=2:nS0
                            for c1=2:nC0
                                x = zeros(sX,1);
                                %1
                                rs = (1+ (c1-1)*ns0+nC0*(s1-1)*ns0);
                                re = rs+ns0-1;
                                x(rs:re) = ones(ns0,1);
                                %2 (c1-1)
                                rs = (1+ (c1-1-1)*ns0+nC0*(s1-1)*ns0);
                                re = rs+ns0-1;
                                x(rs:re) = -ones(ns0,1);
                                %3 (s1-1)
                                rs = (1+ (c1-1)*ns0+nC0*(s1-1-1)*ns0);
                                re = rs+ns0-1;
                                x(rs:re) = -ones(ns0,1);
                                %4 c1-1 and s1-1
                                rs = (1+ (c1-1-1)*ns0+nC0*(s1-1-1)*ns0);
                                re = rs+ns0-1;
                                x(rs:re) = ones(ns0,1);
                                Xab = [Xab x];
                            end
                        end
                        %effect of A at each level of B
                        XeA = [];
                        for c1=1:nC0 %B
                            for s1=2:nS0 %B
                                x = zeros(sX,1);
                                %1
                                rs = (1+ (c1-1)*ns0+nC0*(s1-1)*ns0);
                                re = rs+ns0-1;
                                x(rs:re) = ones(ns0,1);
                                %3 (s1-1)
                                rs = (1+ (c1-1)*ns0+nC0*(s1-1-1)*ns0);
                                re = rs+ns0-1;
                                x(rs:re) = -ones(ns0,1);
                                XeA = [XeA x];
                            end
                        end
                        %effect of B at each level of A
                        XeB = [];
                        for s1=1:nS0 %A
                            for c1=2:nC0 %B
                                x = zeros(sX,1);
                                %1
                                rs = (1+ (c1-1)*ns0+nC0*(s1-1)*ns0);
                                re = rs+ns0-1;
                                x(rs:re) = ones(ns0,1);
                                %2 (c1-1)
                                rs = (1+ (c1-1-1)*ns0+nC0*(s1-1)*ns0);
                                re = rs+ns0-1;
                                x(rs:re) = -ones(ns0,1);
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
                            X0 = [Xa Xb M Xs];
                            X = [Xab X0];
                            A = liom_group_2A(cbeta{h1},X,X0,W.s1,W.s2,Z); %careful, s1 (session counter) not same as W.s1 (size of image)!
                            A = calc_hfgg(cbeta{h1},A,ns0,W.s1,W.s2,B,2);
                            [TOPO H] = call_figure_2anova(TOPO,H,Z,W,F,A,CF,v1,h1,hb,strA,z1,A.LKC);
                        case 2
                            %Main effect of A (Sessions)
                            strA = 'mainA';
                            X0 = [Xab Xb M Xs];
                            X = [Xa X0];
                            A = liom_group_2A(cbeta{h1},X,X0,W.s1,W.s2,Z); %careful, s1 (session counter) not same as W.s1 (size of image)!
                            [TOPO H] = call_figure_2anova(TOPO,H,Z,W,F,A,CF,v1,h1,hb,strA,z1,A.LKC);
                            
                        case 3
                            %Main effect of B (Intensity -- task)
                            strA = 'mainB';
                            X0 = [Xab Xa M Xs];
                            X = [Xb X0];
                            A = liom_group_2A(cbeta{h1},X,X0,W.s1,W.s2,Z); %careful, s1 (session counter) not same as W.s1 (size of image)!
                            [TOPO H] = call_figure_2anova(TOPO,H,Z,W,F,A,CF,v1,h1,hb,strA,z1,A.LKC);
                        case 4
                            %Effect of A within levels of B
                            %loop over levels of B
                            for y1=1:nC0
                                strA = ['effAonB' int2str(y1)];
                                rs = 1+(y1-1)*(nS0-1);
                                re = rs + nS0-2;
                                ind = 1:size(XeA,2);
                                ind(rs:re) = [];
                                X0 = [XeA(:,ind) Xb M Xs];
                                X  = [XeA Xb M Xs];
                                %need to adjust p value:
                                Z.p_value = p_value0/nC0;
                                A = liom_group_2A(cbeta{h1},X,X0,W.s1,W.s2,Z); %careful, s1 (session counter) not same as W.s1 (size of image)!
                                A.includeSubjectEffects = Z.includeSubjectEffects;
                                B = [];
                                B.Xs = Xs;
                                B.Xbw = [Xa Xb]; %Xb; %?
                                A = calc_hfgg(cbeta{h1},A,ns0,B,2);
                                [TOPO H] = call_figure_2anova(TOPO,H,Z,W,F,A,CF,v1,h1,hb,strA,z1,A.LKC);
                            end
                        case 5
                            %Effect of B within levels of A
                            for y1=1:nS0
                                strA = ['effBonA' int2str(y1)];
                                rs = 1+(y1-1)*(nC0-1);
                                re = rs + nC0-2;
                                ind = 1:size(XeB,2);
                                ind(rs:re) = [];
                                X0 = [XeB(:,ind) Xa M Xs];
                                X  = [XeB Xa M Xs];
                                %need to adjust p value:
                                Z.p_value = p_value0/nS0;
                                A = liom_group_2A(cbeta{h1},X,X0,W.s1,W.s2,Z); %careful, s1 (session counter) not same as W.s1 (size of image)!
                                A.includeSubjectEffects = Z.includeSubjectEffects;
                                B = [];
                                B.Xs = Xs;
                                B.Xbw = [Xa Xb]; %Xa; %?
                                A = calc_hfgg(cbeta{h1},A,ns0,B,2);
                                [TOPO H] = call_figure_2anova(TOPO,H,Z,W,F,A,CF,v1,h1,hb,strA,z1,A.LKC);
                            end
                    end
                end %hb
                %save assembled figures
                if ~(z1 == 4) && ~(z1 == 5) %does not work for these cases
                    Z.strA = strA;
                    call_save_assembled_figures(Z,W,H,0);
                end
            end %for z1
        end %if view_estimated
    end %end for v1
    save(ftopo,'TOPO');
    save(fullfile(dir_group,'big_TOPO.mat'),'big_TOPO','-v7.3');
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Could not do 2-anova analysis');
end
out.NIRSmat = job.NIRSmat;