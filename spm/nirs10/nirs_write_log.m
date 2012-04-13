function nirs_write_log(flog,mlog)
    %write to screen and to log file
    disp(mlog);
    for i=1:size(mlog,1);
        fprintf(flog,'%s\n',mlog(i,:));
    end
end
% % % 
% % % %Check if log file exists
% % % if ~exist(DirLog,'dir'), mkdir(DirLog); end
% % % temp_log = fullfile(DirAnalysis,DirLog,['epifMRI_log' date]);
% % % k = 0;
% % % if exist(temp_log,'file')
% % %     k = 2;
% % %     %test = 1;
% % %     while 1
% % %         temp_log2 = [temp_log '_' int2str(k)];
% % %         if exist(temp_log2,'file')
% % %             k = k+1;
% % %         else
% % %             break;
% % %         end
% % %     end
% % % end
% % % 
% % % if k
% % %     log_file = temp_log2;
% % % else
% % %     log_file = temp_log;
% % % end
% % % 
% % % %Open log file
% % % try
% % %     flog = fopen(log_file,'wt');