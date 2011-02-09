function out = nirs_run_liom_group(job)
%Run simple group level analysis as a one sample t-test
FFX_or_RFX = job.FFX_or_RFX;
p_value = job.p_value;
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
if ~FFX_or_RFX
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
                brain = brain(end:-1:1,:);
                brain = brain * 0.5;
                %Loop over chromophores
                for h1=1:2 %exclude HbT for now
                    hb = get_chromophore(h1);
                    contrasts = TOPO.xCon;
                    for c1=1:length(contrasts)
                        %run group analysis
                        ns = length(TOPO.v{v1}.s);
                        if ns > 1
                            %Generate group result as t-stat 
                            min_s = 2;
                            [tmap_group, erdf_group] = nirs_spm_group(TOPO.v{v1}.s,ns,min_s);
                            %fill in SPM_nirs fields

                            TOPO.v{v1}.group.ns = ns;
                            TOPO.v{v1}.group.min_s = min_s;
                            switch h1
                                case 1
                                    TOPO.v{v1}.group.HbO.Tmap = tmap_group;
                                    TOPO.v{v1}.group.HbO.erdf = erdf_group;
                                case 2
                                    TOPO.v{v1}.group.HbR.Tmap = tmap_group;
                                    TOPO.v{v1}.group.HbR.erdf = erdf_group;
                            end
                                                                        
                            load split
                            nirs_draw_figure(2,brain,tmap_group,info1,...
                                info_for_fig1,split,pathn,erdf_group,[],p_value,gen_fig,gen_tiff)

                        end
                    end
                end



%Find keywords in the file names
%Chromophores: HbO, HbR, HbT
%View: ventral, dorsal, right, left, frontal, occipital
%T-stat side: positive or negative
%Subject or session number: not required (take all)
fn = job.Contrast_files;
%contrasts=[{'positive'};{'negative'}];
contrasts=[{'_Can_'};{'_negCanXcan'};{'_canXcan'};{'_negCan_'}]; 
[dir_spm fil1] = fileparts(fn{1});

%Create a tensorial array (chromophores by views by contrasts) to store
%file names found
fs{3,6,size(contrasts,1)}={};
%Loop over Chromophores
for h1=1:3
    hb = get_chromophore(h1);
    %Loop over views
    for v1=1:6
        vi = get_view(v1);
        %Loop over contrasts
        for c1=1:size(contrasts,1)
            %Loop over user-specified files
            for f1=1:size(fn,1)
                if findstr(lower(hb),lower(fn{f1}))
                    if findstr(lower(vi),lower(fn{f1}))
                        if findstr(lower(contrasts{c1}),lower(fn{f1}))
                        %if findstr(contrasts{c1},fn{f1})
                            %add the file name to the fs cell tensorial array
                            if ~iscell(fn{f1})
                                fs{h1,v1,c1}=[fs{h1,v1,c1};{fn{f1}}];
                            else
                                fs{h1,v1,c1}=[fs{h1,v1,c1};fn{f1}];
                            end
                        end
                    end
                end
            end %end for f1
        end %end for c1
    end %end for v1
end %end for h1

%To generate figures
load([spm('dir') filesep 'rend' filesep 'render_single_subj.mat']);
load Split
ResultsFile = fullfile(dir_spm,'ResultsGroupAll.ps');

for v1=1:6
    vi = get_view(v1);
    %View dependent info for figures
    p_value = 0.05;
    flag_figure = 1;        
    brain = rend{v1}.ren;
    if issparse(brain),
        d = size(brain);
        B1 = spm_dctmtx(d(1),d(1));
        B2 = spm_dctmtx(d(2),d(2));
        brain = B1*brain*B2';
    end;
    msk = brain>1;brain(msk)=1;
    msk = brain<0;brain(msk)=0;
    brain = brain(end:-1:1,:);
    brain = brain * 0.5;
    %Loop over chromophores
    for h1=1:3
        hb = get_chromophore(h1);
        for c1=1:size(contrasts,1)
            SPM_nirs = [];
            %run group analysis
            nsubj = size(fs{h1,v1,c1},1);
            if nsubj > 1
                %Generate group result as t-stat 
                min_subj = 2;
                [tmap_group, erdf_group] = nirs_spm_group(fs{h1,v1,c1},nsubj,min_subj);
                %fill in SPM_nirs fields
                SPM_nirs.nirs.level = 'group';
                SPM_nirs.nirs.tstat = tmap_group;
                SPM_nirs.nirs.erdf = erdf_group;
                SPM_nirs.nirs.nsubj = nsubj;
                SPM_nirs.nirs.min_subj = min_subj;
                SPM_nirs.nirs.Hb = hb;
                SPM_nirs.nirs.spec_hemi = vi;
                SPM_nirs.nirs.side_hemi = v1;
                info1 = ['T_SPM_group_' hb '_' vi '_' contrasts{c1}];
                info_for_fig1 = ['T SPM group ' hb ' ' vi ' ' contrasts{c1}];
                save(fullfile(dir_spm,[info1 '.mat']), 'SPM_nirs');
                %Generate figures                
                T_map = SPM_nirs.nirs.tstat;
                erdf = SPM_nirs.nirs.erdf;                
                %Uncorrected
                str_cor = 'unc';
                th_z = spm_invTcdf(1-p_value, erdf);
                fh1 = nirs_SPM_NIRS_draw_figure(th_z,brain,info1,...
                    info_for_fig1,T_map,flag_figure,str_cor,split);
                filen1 = fullfile(dir_spm,[info1 '_' str_cor '.fig']);
                print(fh1,'-dpsc2','-append', ResultsFile);
                saveas(fh1,filen1,'fig');
                try close(fh1); end
            end
        end
    end
end
        
out = [];
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

function vi = get_view(v1)
    switch v1
        case 1
            vi='ventral';
        case 2
            vi='dorsal';
        case 3
            vi='right';
        case 4
            vi='left';
        case 5
            vi='frontal';
        case 6
            vi='occipital';
    end
end

