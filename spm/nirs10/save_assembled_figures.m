function save_assembled_figures(Z,W,fh0,Inv,str_cor,f1)
%All this does is save a previously constructed arrangement of figures
if f1 > 0
    strf = ['Sess' int2str(f1) '_'];
    str_grp = '';
else
    strf = '';
    str_grp = 'Group_';
end
if ~isempty(Inv)
    Inv = ['_' Inv];
end
if isfield(Z,'scon')
    scon = ['_' Z.scon];
else
    scon = '';
end
%'Sess'
if isfield(Z,'Idx')
    subj_str = ['S' gen_num_str(Z.Idx,2) '_'];
    try
        if isfield(Z,'subj_id')
            %This is useful if the subject numbers are not preserved due to
            %some subjects being dropped for various reasons
            subj_str = ['ID' Z.subj_id subj_str];
        end
    end
else
    subj_str = ''; %group analysis
end
filestr = [subj_str str_grp strf str_cor '_' num2str(Z.p_value) '_' W.spec_hemi Inv scon];

filen1 = fullfile(Z.dir1,[filestr '.fig']);
if Z.gen_fig
    saveas(fh0,filen1,'fig'); %save as .fig
end
if Z.gen_tiff %label 'A' conventiently brings figures at beginning of folder list
    filen2 = fullfile(Z.dir1,['A_' filestr '.tiff']); %save as .tiff
    print(fh0, '-dtiffn', filen2);
end
if strcmp(Z.cbar.visible, 'off')
    nirs_close_figure(fh0);
end