function nirs_write_log(flog,mlog)
    %write to screen and to log file
    disp(mlog);
    for i=1:size(mlog,1);
        fprintf(flog,'%s\n',mlog(i,:));
    end
end