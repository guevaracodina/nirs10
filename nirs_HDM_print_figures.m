function nirs_HDM_print_figures(HDM,h,fname)
fname = [fname 'Su' HDM.subj_id 'S' int2str(HDM.s1) 'Ch' gen_num_str(HDM.c1,3)];
if isfield(HDM,'print_it') && HDM.print_it
    fname = [fname 'I' gen_num_str(HDM.it.k,3)];
end
fig_dir = fullfile(HDM.dir1,'fig');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
filen1 = fullfile(fig_dir,[fname '.fig']); 
saveas(h,filen1,'fig'); %save as .fig
filen2 = fullfile(fig_dir,[fname '.png']);
print(h,'-dpng', filen2,'-r300');
if ~HDM.DO.generate_figures
    try
        close(h); 
    catch
        disp('Problem closing figure');
    end
end