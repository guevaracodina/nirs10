function save_assembled_figures(Z,W,fh0,Inv,str_cor,f1) 
%All this does is save a previously constructed arrangement of figures
if f1 > 0
    strf = ['_S' int2str(f1)];
    str_grp = '';
else
    strf = '';
    str_grp = 'Group_';
end
if ~isempty(Inv)
    Inv = ['_' Inv];
end
filestr = [str_grp str_cor '_' num2str(Z.p_value) '_' W.spec_hemi strf Inv];

filen1 = fullfile(Z.dir1,[filestr '.fig']);
if Z.gen_fig
    saveas(fh0,filen1,'fig'); %save as .fig
end
if Z.gen_tiff %label 'A' conventiently brings figures at beginning of folder list  
    filen2 = fullfile(Z.dir1,['A_' filestr '.tiff']); %save as .tiff
    print(fh0, '-dtiffn', filen2);
end
if strcmp(Z.cbar.visible, 'off')
    try close(fh0); end
end