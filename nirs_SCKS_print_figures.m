function nirs_SCKS_print_figures(SCKS,h,fname)
fname = [fname 'Su' SCKS.subj_id 'S' int2str(SCKS.s1) 'Ch' gen_num_str(SCKS.c1,3)];
fig_dir = fullfile(SCKS.dir1,'fig');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
filen1 = fullfile(fig_dir,[fname '.fig']); 
saveas(h,filen1,'fig'); %save as .fig
filen2 = fullfile(fig_dir,[fname '.png']);
print(h,'-dpng', filen2,'-r300');
if ~SCKS.DO.generate_figures
    try
        close(h); 
    catch
        disp('Problem closing figure');
    end
end