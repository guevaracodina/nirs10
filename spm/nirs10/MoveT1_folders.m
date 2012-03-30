%move t1 folders
for i0=2:89
    try
        oldpath = ['W:\Claudine\SPMDataNT\S' gen_num_str(i0,3) '\T1'];
        newpath = ['W:\Claudine\SPMDataL\S' gen_num_str(i0,3)];
        if ~exist(newpath,'dir'), mkdir(newpath); end
        copyfile(oldpath,fullfile(newpath,'\T1'))
    catch
        disp(['problem for subject ' int2str(i0)])
    end
end

