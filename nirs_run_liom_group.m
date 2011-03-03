function out = nirs_run_liom_group(job)
%Run simple group level analysis as a one sample t-test
FFX = job.FFX_or_RFX;
p_value = job.contrast_p_value;
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
    cbar.colorbar_override = 1;
catch
    cbar.colorbar_override = 0;
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
if FFX || size(job.NIRSmat,1)==1 
    %fixed effects: loop over subjects first, as they are treated
    %separately
    %Loop over all subjects
    for Idx=1:size(job.NIRSmat,1)
        %Load NIRS.mat information
        try
            NIRS = [];
            load(job.NIRSmat{Idx,1});
            %load SPM - first GLM - might want to generalize 
            dir1 = NIRS.SPM{1};
            %load topographic information (formerly known as preproc_info)
            fname_ch = NIRS.Dt.ana.rend;
            load(fname_ch);
            %load(fullfile(dir1,'SPM.mat'));
            ftopo = fullfile(dir1,'TOPO.mat');
            TOPO = [];
            load(ftopo); 

            %Big loop over views 
            for v1=1:6
                view_estimated = 0;
                try
                    s1 = TOPO.v{v1}.s1;
                    s2 = TOPO.v{v1}.s2;
                    view_estimated = 1;
                catch
                end
                if view_estimated
                    switch v1
                        case 1 % 'ventral'
                            spec_hemi = 'ventral';
                            side_hemi = 1;
                        case 2 % 'dorsal'
                            spec_hemi = 'dorsal';
                            side_hemi = 2;
                        case 3 %'right_lateral'
                            spec_hemi = 'right';
                            side_hemi = 3;
                        case 4 %'left_lateral'
                            spec_hemi = 'left';
                            side_hemi = 4;
                        case 5 %'frontal'
                            spec_hemi = 'frontal';
                            side_hemi = 5;
                        case 6 %'occipital'
                            spec_hemi = 'occipital';
                            side_hemi = 6;
                    end

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
                    %Contrasts
                    xCon = TOPO.xCon;              

                    ns = length(TOPO.v{v1}.s);
                    min_s = 2;
                    TOPO.v{v1}.group.ns = ns;
                    TOPO.v{v1}.group.min_s = min_s;
                    TOPO.v{v1}.group.s1 = s1;
                    TOPO.v{v1}.group.s2 = s2;
                    cbeta = zeros(ns,s1*s2);
                    ccov_beta = zeros(ns,s1*s2);
                    tmp = zeros(s1,s2);
                    nC = length(xCon);
                    if GFIS
                        fh0P = figure('Visible',cbar.visible,'Name',['Group' '_' num2str(p_value) '_Pos'],'NumberTitle','off');
                        %subplot(fh0P,nC,3,1);
                        if GInv
                            fh0N = figure('Visible',cbar.visible,'Name',['Group' '_' num2str(p_value) '_Neg'],'NumberTitle','off');
                            %subplot(fh0N,nC,3,1);
                        else
                            fh0N = [];
                        end
                    end
                    load Split
                    if ns > 1
                        %Loop over chromophores
                        for h1=1:3 %including HbT
                            hb = get_chromophore(h1);
                            for c1=1:nC   
                                %fill in cbeta and ccov_beta                            
                                for f1=1:ns
                                    tmp = squeeze(TOPO.v{v1}.s{f1}.hb{h1}.c_interp_beta(c1,:,:));
                                    cbeta(f1,:) = tmp(:);
                                    tmp = squeeze(TOPO.v{v1}.s{f1}.hb{h1}.c_cov_interp_beta(c1,:,:));
                                    ccov_beta(f1,:) = tmp(:);
                                end
                                %Positive contrasts
                                %Generate group result as t-stat 
                                if GInv || strcmp(hb,'HbO') || strcmp(hb,'HbT')
                                    [tmap_group, erdf_group, var_bs, beta_group] = liom_group(...
                                                cbeta,ccov_beta,s1,s2,ns,min_s,FFX);

                                    TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.Tmap = tmap_group;
                                    TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.erdf = erdf_group;
                                    TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.beta_group = beta_group;
                                    TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.type = 'Positive';
                                    TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.var_bs = var_bs;
                                    TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.c = xCon(c1);

                                    info1 = [num2str(p_value) '_' spec_hemi '_' hb '_Pos' xCon(c1).n];
                                    info_for_fig1 = [num2str(p_value) ' ' spec_hemi ' ' hb ' Pos' xCon(c1).n];
                                    erdf_group = max(erdf_group(:)); %quick fix...
                                    [fh1 ax1 hc1] = nirs_draw_figure(4,brain,tmap_group,info1,...
                                        info_for_fig1,split,dir1,erdf_group,[],p_value,gen_fig,gen_tiff,cbar);
                                    if GFIS, [fh0P, fh0N] = nirs_copy_figure(fh0P,fh0N,c1,nC,hb,GInv,fh1,ax1,hc1,1,split); end                                   
                                end
                                %Negative contrasts
                                for f1=1:ns
                                    tmp = -squeeze(TOPO.v{v1}.s{f1}.hb{h1}.c_interp_beta(c1,:,:));
                                    cbeta(f1,:,:) = tmp(:);
                                end
                                if GInv || strcmp(hb,'HbR') 
                                    %Generate group result as t-stat 
                                    [tmap_group, erdf_group, var_bs, beta_group] = liom_group(...
                                                cbeta,ccov_beta,s1,s2,ns,min_s,FFX);

                                    TOPO.v{v1}.group.hb{h1}.c{2*c1}.Tmap = tmap_group;
                                    TOPO.v{v1}.group.hb{h1}.c{2*c1}.erdf = erdf_group;
                                    TOPO.v{v1}.group.hb{h1}.c{2*c1}.beta_group = beta_group;
                                    TOPO.v{v1}.group.hb{h1}.c{2*c1}.type = 'Negative';
                                    TOPO.v{v1}.group.hb{h1}.c{2*c1}.var_bs = var_bs;
                                    TOPO.v{v1}.group.hb{h1}.c{2*c1}.c = xCon(c1);

                                    info1 = [num2str(p_value) '_' spec_hemi '_' hb '_Neg' xCon(c1).n];
                                    info_for_fig1 = [num2str(p_value) ' ' spec_hemi ' ' hb ' Neg' xCon(c1).n];
                                    erdf_group = max(erdf_group(:)); %quick fix...
                                    [fh1 ax1 hc1] = nirs_draw_figure(4,brain,tmap_group,info1,...
                                        info_for_fig1,split,dir1,erdf_group,[],p_value,gen_fig,gen_tiff,cbar);
                                    if GFIS, [fh0P, fh0N] = nirs_copy_figure(fh0P,fh0N,c1,nC,hb,GInv,fh1,ax1,hc1,0,split); end
                                end
                            end
                        end
                        %save assembled figures
                        if GFIS
                            save_group_figures(fh0P, gen_fig,gen_tiff,p_value,'Pos',dir1,spec_hemi,cbar);
                            if GInv
                                save_group_figures(fh0N, gen_fig,gen_tiff,p_value,'Neg',dir1,spec_hemi,cbar);
                            end
                        end                         
                    end
                end %if view_estimated
            end %end for v1
            %TOPO.xCon = xCon; %would not work if new contrasts are later added        
            save(ftopo,'TOPO');
        catch
            disp(['Could not do FFX group analysis for subject' int2str(Idx)]);
        end
        %NIRS.TOPO = ftopo;
        %save(job.NIRSmat{Idx,1});
    end
else
    %RFX    
    %Loop over all subjects - load a large amount of data - might be too
    %much for many subjects 
    %number of subjects
    ns = size(job.NIRSmat,1);
    %minimum number of subjects for thresholding tmaps
    min_s = 2;
    big_TOPO{ns} =[];
    for Idx=1:ns
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
    %Contrasts - assume same contrasts for all subjects
    xCon = big_TOPO{1}.xCon;              

    %create a new TOPO at the group level
    TOPO = [];
    try
        %Big loop over views 
        for v1=1:6
            view_estimated = 0;
            try
                s1 = big_TOPO{1}.v{v1}.s1;
                s2 = big_TOPO{1}.v{v1}.s2;
                view_estimated = 1;
            catch
            end
            if view_estimated
                switch v1
                    case 1 % 'ventral'
                        spec_hemi = 'ventral';
                        side_hemi = 1;
                    case 2 % 'dorsal'
                        spec_hemi = 'dorsal';
                        side_hemi = 2;
                    case 3 %'right_lateral'
                        spec_hemi = 'right';
                        side_hemi = 3;
                    case 4 %'left_lateral'
                        spec_hemi = 'left';
                        side_hemi = 4;
                    case 5 %'frontal'
                        spec_hemi = 'frontal';
                        side_hemi = 5;
                    case 6 %'occipital'
                        spec_hemi = 'occipital';
                        side_hemi = 6;
                end

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

                TOPO.v{v1}.group.ns = ns;
                TOPO.v{v1}.group.min_s = min_s;
                TOPO.v{v1}.group.s1 = s1;
                TOPO.v{v1}.group.s2 = s2;
                cbeta = zeros(ns,s1*s2);
                ccov_beta = zeros(ns,s1*s2);
                tmp = zeros(s1,s2);
                nC = length(xCon);
                if GFIS
                    fh0P = figure('Visible',cbar.visible,'Name',['A_Group' '_' num2str(p_value) '_Pos'],'NumberTitle','off');
                    %subplot(fh0P,nC,3,1);
                    if GInv
                        fh0N = figure('Visible',cbar.visible,'Name',['A_Group' '_' num2str(p_value) '_Neg'],'NumberTitle','off');
                        %subplot(fh0N,nC,3,1);
                    else
                        fh0N = [];
                    end
                end
                load Split
                %Loop over chromophores
                for h1=1:3 %exclude HbT for now
                    hb = get_chromophore(h1);
                    %Loop over contrasts
                    for c1=1:nC   
                        %fill in cbeta and ccov_beta   
                        %Loop over subjects
                        for f1=1:ns 
                            %assume only one session
                            tmp = squeeze(big_TOPO{f1}.v{v1}.s{1}.hb{h1}.c_interp_beta(c1,:,:));
                            cbeta(f1,:) = tmp(:);
                            tmp = squeeze(big_TOPO{f1}.v{v1}.s{1}.hb{h1}.c_cov_interp_beta(c1,:,:));
                            ccov_beta(f1,:) = tmp(:);
                        end
                        %Positive contrasts
                        %Generate group result as t-stat 
                        if GInv || strcmp(hb,'HbO') || strcmp(hb,'HbT')
                            [tmap_group, erdf_group, var_bs, beta_group] = liom_group(...
                                        cbeta,ccov_beta,s1,s2,ns,min_s,FFX);

                            TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.Tmap = tmap_group;
                            TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.erdf = erdf_group;
                            TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.beta_group = beta_group;
                            TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.type = 'Positive';
                            TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.var_bs = var_bs;
                            TOPO.v{v1}.group.hb{h1}.c{2*c1-1}.c = xCon(c1);

                            info1 = [num2str(p_value) '_' spec_hemi '_' hb '_Pos' xCon(c1).n];
                            info_for_fig1 = [num2str(p_value) ' ' spec_hemi ' ' hb ' Pos' xCon(c1).n];
                            erdf_group = max(erdf_group(:)); %quick fix...
                            [fh1 ax1 hc1] = nirs_draw_figure(4,brain,tmap_group,info1,...
                                info_for_fig1,split,dir1,erdf_group,[],p_value,gen_fig,gen_tiff,cbar);
                            if GFIS, [fh0P, fh0N] = nirs_copy_figure(fh0P,fh0N,c1,nC,hb,GInv,fh1,ax1,hc1,1,split); end
                        end
                        %Negative contrasts
                        for f1=1:ns
                            tmp = -squeeze(big_TOPO{f1}.v{v1}.s{1}.hb{h1}.c_interp_beta(c1,:,:));
                            cbeta(f1,:,:) = tmp(:);
                        end
                        if GInv || strcmp(hb,'HbR') 
                            %Generate group result as t-stat 
                            [tmap_group, erdf_group, var_bs, beta_group] = liom_group(...
                                        cbeta,ccov_beta,s1,s2,ns,min_s,FFX);

                            TOPO.v{v1}.group.hb{h1}.c{2*c1}.Tmap = tmap_group;
                            TOPO.v{v1}.group.hb{h1}.c{2*c1}.erdf = erdf_group;
                            TOPO.v{v1}.group.hb{h1}.c{2*c1}.beta_group = beta_group;
                            TOPO.v{v1}.group.hb{h1}.c{2*c1}.type = 'Negative';
                            TOPO.v{v1}.group.hb{h1}.c{2*c1}.var_bs = var_bs;
                            TOPO.v{v1}.group.hb{h1}.c{2*c1}.c = xCon(c1);

                            info1 = [num2str(p_value) '_' spec_hemi '_' hb '_Neg' xCon(c1).n];
                            info_for_fig1 = [num2str(p_value) ' ' spec_hemi ' ' hb ' Neg' xCon(c1).n];
                            erdf_group = max(erdf_group(:)); %quick fix...
                            [fh1 ax1 hc1] = nirs_draw_figure(4,brain,tmap_group,info1,...
                                info_for_fig1,split,dir1,erdf_group,[],p_value,gen_fig,gen_tiff,cbar);
                            if GFIS, [fh0P, fh0N] = nirs_copy_figure(fh0P,fh0N,c1,nC,hb,GInv,fh1,ax1,hc1,0,split); end
                        end
                    end
                end
                if GFIS
                    save_group_figures(fh0P, gen_fig,gen_tiff,p_value,'Pos',dir1,spec_hemi,cbar);
                    if GInv
                        save_group_figures(fh0N, gen_fig,gen_tiff,p_value,'Neg',dir1,spec_hemi,cbar);
                    end
                end 
            end %if view_estimated
        end %end for v1
        TOPO.xCon = xCon; %would not work if new contrasts are later added 
        [dir0,~,~] = fileparts(job.NIRSmat{1});
        %store in same directory as first subject
        ftopo = fullfile(dir0,'TOPO.mat');
        save(ftopo,'TOPO');
    catch
        disp(['Could not do FFX group analysis for subject' int2str(Idx)]);
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

function save_group_figures(fh0,gen_fig,gen_tiff,p_value,Inv,dir1,spec_hemi,cbar) 
filen1 = fullfile(dir1,['Group_' num2str(p_value) '_' spec_hemi '_' Inv '.fig']);
if gen_fig
    saveas(fh0,filen1,'fig');
end
if gen_tiff
    filen2 = fullfile(dir1,['A_Group_' num2str(p_value) '_' spec_hemi '_' Inv '.tiff']);
    print(fh0, '-dtiffn', filen2);
end
if strcmp(cbar.visible, 'off')
    try close(fh0); end
end
end