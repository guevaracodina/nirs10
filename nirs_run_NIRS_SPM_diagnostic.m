function out = nirs_run_NIRS_SPM_diagnostic(job)
% 'Contrast estimates and 90% C.I.'
%data_files = job.data_files;
ch_num = job.ch_num;

%SPM_nirs = [];
%nirs_data = [];
%loop over data files
for f1=1:size(job.data_files,1)
    data_file = job.data_files{f1,1};
    try
        load(data_file);
        SPM_nirs;
    catch
        disp('Could not load SPM_NIRS data file');
    end
    try
        load(SPM_nirs.nirs.fname);
        nirs_data;
    catch
        disp('Could not load HbO/HbR concentration file');
    end

    if findstr('HbO',data_file)
        hb = 1;
        hbstr = 'HbO';
    else
        if findstr('HbR',data_file)
            hb = 2;
            hbstr = 'HbR';
        else
            if findstr('HbT',data_file)
                hb = 3;
                hbstr = 'HbT';
            else
                disp('Unrecognized data type');
                out = [];
                return
            end
        end
    end
    switch hb
        case 1    
            d = nirs_data.oxyData(:,ch_num);        
        case 2
            d = nirs_data.dxyData(:,ch_num);
        case 3
            d = nirs_data.oxyData(:,ch_num) + nirs_data.dxyData(:,ch_num);
    end
    %Filtered data - not valid for HbT, for which oxy & dxy need to be filtered
    %separately
    for i2=1:size(d,2)
         ky(:,i2) = spm_filter_HPF_LPF_WMDL(SPM_nirs.xX.K.KL, d(:,i2));
    end
    n_U = size(SPM_nirs.Sess.U,2);
    for i1=1:n_U
        U(:,i1) = SPM_nirs.Sess.U(1,i1).u;
    end
    %seconds
    l1 = linspace(0,size(d,1)*SPM_nirs.xY.RT,size(d,1));
    l2 = linspace(0,size(U,1)*SPM_nirs.Sess.U(1,1).dt,size(U,1));
    %Plot results
    for i2=1:size(d,2)
        figure;
        title([hbstr]);
        %Black: raw data
        plot(l1,d(:,i2),'k'); hold on
        %Red: detrended data
        plot(l1,ky(:,i2),'r'); hold on
        %Blue: first onsets
        plot(l2,U(:,1),'b'); hold on
        plot(l2,-U(:,1),'b'); hold on
        %if more than one type of onsets, only plot the first three
        if size(U,2) > 1 %green: second onsets
            plot(l2,U(:,2),'g'); 
            plot(l2,-U(:,2),'g'); 
            if size(U,2) > 2 %yellow: third onsets
                plot(l2,U(:,3),'y'); 
                plot(l2,-U(:,3),'y'); 
            end
        end
        hold off
    end
end
%------------------------------------------------------------------
out = [];
end