function out = nirs_run_NIRS_SPM_group(job)
%Run simple group level analysis as a one sample t-test

%To generate F-test at group of sessions level for canXcan contrast
canXcanOn = 1; %1;

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
                %Generate group results as F-stat - use for canXcan ONLY!
                if canXcanOn
                    [Fmap_group erdf_group] = group_Fstat_canXcan(fs{h1,v1,c1},nsubj,min_subj);                  
                    %fill in SPM_nirs fields
                    SPM_nirs.nirs.level = 'group';
                    SPM_nirs.nirs.Fstat = Fmap_group;
                    SPM_nirs.nirs.erdf = erdf_group;
                    SPM_nirs.nirs.nsubj = nsubj;
                    SPM_nirs.nirs.min_subj = min_subj;
                    SPM_nirs.nirs.Hb = hb;
                    SPM_nirs.nirs.spec_hemi = vi;
                    SPM_nirs.nirs.side_hemi = v1;
                    info1 = ['F_SPM_group_' hb '_' vi '_' contrasts{c1}];
                    info_for_fig1 = ['F SPM group ' hb ' ' vi ' ' contrasts{c1}];
                    save(fullfile(dir_spm,[info1 '.mat']), 'SPM_nirs');
                    %Generate figures                
                    F_map = SPM_nirs.nirs.Fstat;
                    erdf = SPM_nirs.nirs.erdf;                
                    %Uncorrected
                    str_cor = 'unc';
                    th_z = spm_invFcdf(1-p_value, nsubj, erdf); %PP
                    fh1 = nirs_SPM_NIRS_draw_figure(th_z,brain,info1,...
                        info_for_fig1,F_map,flag_figure,str_cor,split);
                    filen1 = fullfile(dir_spm,[info1 '_' str_cor '.fig']);
                    print(fh1,'-dpsc2','-append', ResultsFile);
                    saveas(fh1,filen1,'fig');
                    try close(fh1); end
                end
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

function [Fmap_group erdf_group] = group_Fstat_canXcan(fname_cinterp_SPM_nirs, nsubj,min_subj)

for kk = 1:nsubj
    load(fname_cinterp_SPM_nirs{kk});
    if kk == 1
        msk = zeros(1, cinterp_SPM_nirs.s1 * cinterp_SPM_nirs.s2);
    end
    index = find(cinterp_SPM_nirs.cbeta ~= 0);
    msk(index) = msk(index) + 1;
    clear index
end

s1 = cinterp_SPM_nirs.s1;
s2 = cinterp_SPM_nirs.s2;

index_msk = find(msk > min_subj-1);
indiv_beta2 = zeros(nsubj, length(index_msk));
indiv_cov  = zeros(nsubj, length(index_msk));

for kk = 1:nsubj
    load(fname_cinterp_SPM_nirs{kk});
    indiv_beta2(kk,:) = cinterp_SPM_nirs.cbeta(index_msk).^2;
    indiv_cov(kk,:) = cinterp_SPM_nirs.ccov_beta(index_msk);
end

%load number of degrees of freedom
load(cinterp_SPM_nirs.fname_SPM);
erdf = SPM_nirs.xX.erdf;

                    
%avg_beta = sum(indiv_beta)./msk(index_msk);

% h_wait = waitbar(0, 'Please wait... ');
% for i = 1:length(index_msk)
%     waitbar(i/length(index_msk))
%     %X = ones(msk(index_msk(i)),1);
%     %index = find(indiv_beta(:, i) ~= 0);
%     %beta_tmp = indiv_beta(index, i);
%     %res = beta_tmp - X *avg_beta(i);
%     %if length(X) == 1
%         %var_bs = 0;
%     %    disp('err');
%     %    return
%     %else
%         %erdf(1, i) = msk(index_msk(i)) - rank(X);
%         %var_bs = sum(res.^2)./erdf(1,i); % variance between subjects
%     %end
%     %Extra-sum-of-squares of full model: sum of each session ESS
%     denum = sum(indiv_cov(:, i)); %indiv_cov(index, i) + var_bs; %%% covariance individual + inter- subject variance
%     %Degrees of freedom
%     %n = ;
%     %p = ;
%     %p0 = ;
%     
%     numer = sum(indiv_beta(:, i).^2);
%     %./tmp_denum;
%     %numer = sum(numer);
%     %denum = sqrt(sum(1./tmp_denum));
%     Fmap(1,i) = numer/(denum/erdf);
% end
Fmap = sum(indiv_beta2,1)./sum(indiv_cov,1);
Fmap = Fmap/nsubj;
%Fmap = Fmap*erdf;
%close(h_wait);

%clear indiv_beta2
%clear indiv_cov

Fmap_group = zeros(1, s1*s2);
Fmap_group(index_msk) = Fmap(:);
Fmap_group = reshape(Fmap_group, s1, s2);

erdf_group = erdf;

% erdf_group = zeros(1, s1*s2);
% erdf_group(index_msk) = erdf(:);
% erdf_group = reshape(erdf_group, s1, s2);



end
