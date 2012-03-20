function [th_z str_cor index_over index_over2] = nirs_get_threshold(fign,F,tstr,erdf,eidf,...
    s_map,GInv,p_value,StatStr,StatStr2,GroupColorbars,G,UseCorrelRes,sum_kappa,nchn,LKCflag)
    if ~isstruct(G)
        LKC = G;
    else
        LKC = G.LKC;
    end
switch fign
    %     case 1
    %         %no threshold
    %         str_cor = 'all';
    %         th_z = 0;
    %         index_over = find(T_map);
    case 2
        %uncorrected
        str_cor = 'unc';
        %threshold
        if tstr == 'T'
            th_z = spm_invTcdf(1-p_value, erdf);
        else
            th_z = spm_invFcdf(1-p_value, eidf,erdf);
        end
        index_over = find(s_map > th_z);
        %if GInv, index_over2 = find(-T_map > th_z); end
        if GInv,
            if GroupColorbars
                index_over2 = [find(s_map > th_z); find(-s_map > th_z)];
            else
                index_over2 = find(-s_map > th_z);
            end
        end
    case 3
        %%%%%%%%%%%%%%%%%
        %Get threshold
        if LKCflag || UseCorrelRes
            th_z = calc_EC(LKC,p_value,tstr,[eidf,erdf]);
            str_cor = StatStr;
        else
            if tstr == 'T'
                %Tube
                str_cor = StatStr;
                %find threshold
                z_value = 1:0.0001:7;
                p_value_tube = ((sum_kappa * gamma(3/2))/(2*(pi^(3/2))))*...
                    (1-gammainc((z_value(:).^2)/2, 3/2));
                index_z = [];
                ini_ran = 10^(-10);
                n = 0;
                while isempty(index_z) == 1
                    ran = ini_ran * (10^n);
                    n = n+1;
                    index_z = find(p_value_tube > p_value - ran & p_value_tube < p_value + ran);
                end
                index_z = index_z(end);
                th_z = z_value(index_z);
            else
                %For F-contrast, do a Bonferroni correction for now
                str_cor = StatStr2;
                p_value = p_value/nchn;
                th_z = spm_invFcdf(1-p_value, eidf,erdf);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%
        %Mask by threshold value
        if tstr == 'T'
            index_over = find(s_map > th_z);
            %if GInv, index_over2 = find(-T_map > th_z); end
            if GInv,
                if GroupColorbars
                    index_over2 = [find(s_map > th_z); find(-s_map > th_z)];
                else
                    index_over2 = find(-s_map > th_z);
                end
            end
        else
            index_over = find(s_map > th_z);
        end
        
    case 4
        %group
        
        if strcmp(tstr,'F')
            %if length(erdf) > 1 %?
            str_cor = [StatStr '_Anova'];
            xth_z = zeros(1,length(erdf));
            index_over = [];
            xerdf = erdf(:);
            udf = unique(xerdf);
            xth = zeros(1,length(udf));
            for i1=2:length(udf) %exclude 0
                if LKCflag || UseCorrelRes
                    xth(i1) = calc_EC(LKC,p_value,tstr,[eidf,udf(i1)]);
                else
                    xth(i1) = spm_invFcdf(1-p_value, eidf, udf(i1));
                end
            end
            s_map = s_map(:); %??
            for i1 = 1:length(xerdf)
                if xerdf(i1) > 0
                    xth_z(i1) = xth(xerdf(i1)==udf);
                    if s_map(i1) > xth_z(i1)
                        index_over = [index_over i1];
                    end
                end
            end
            th_z = min(xth_z(xth_z>0));
            %end
        else
            str_cor = [StatStr '_Group'];
            if erdf == 0
                disp([F.contrast_info ' : Problem with number of degrees of freedom']);
                index_over = [];
                if GInv
                    index_over2 = [];
                end
                th_z = 0;
            else
                if LKCflag
                    test_LKC_by_pixel = 0;
                    if test_LKC_by_pixel && isfield(G,'ns')
                        th_zi = zeros(1,G.ns-G.min_subj+1);
                        for i0=1:(G.ns-G.min_subj+1)
                            th_zi(i0) = calc_EC(G.LKCi{i0},p_value,tstr,[G.eidfi(i0),G.erdfi(i0)]);
                        end
                        th_z = max(th_zi(1:end/2)); %completely heuristic
                        %a=[]; for i0=1:(G.ns-G.min_subj+1), a = [a length(G.idxi{i0})]; end, a
                    else
                        th_z = calc_EC(LKC,p_value,tstr,[eidf,erdf]);
                    end
                    str_cor = StatStr;
                else
                    th_z = spm_invTcdf(1-p_value, erdf);
                end
                index_over = find(s_map > th_z);
                %if GInv, index_over2 = find(-T_map > th_z); end
                if GInv,
                    if GroupColorbars
                        index_over2 = [find(s_map > th_z); find(-s_map > th_z)];
                    else
                        index_over2 = find(-s_map > th_z);
                    end
                end
            end
        end
    case 5
        %for 2-anova -- LKC
        str_cor = [StatStr '_2A' int2str(eidf) '_' int2str(erdf-eidf)];
        if strcmp(F.tstr,'F') % should always be the case
            th_z = calc_EC(LKC,p_value,tstr,[eidf,erdf]);
            index_over = find(s_map > th_z);
            index_over2 = []; %not used
        end
    case 7
        %for 2-anova -- uncorrected
        str_cor = ['unc' '_2A' int2str(eidf) '_' int2str(erdf-eidf)];
        if strcmp(F.tstr,'F') % should always be the case
            th_z = spm_invFcdf(1-p_value, eidf, erdf);
            index_over = find(s_map > th_z);
            index_over2 = []; %not used
        end
    case 6 %group, uncorrected
        str_cor = 'unc_Group';
        %threshold
        if tstr == 'T'
            th_z = spm_invTcdf(1-p_value, erdf);
        else
            th_z = spm_invFcdf(1-p_value, eidf,erdf);
        end
        if erdf == 0
            disp([F.contrast_info ' : Problem with number of degrees of freedom']);
            index_over = [];
            if GInv
                index_over2 = [];
            end
            th_z = 0;
        else
            index_over = find(s_map > th_z);
            %if GInv, index_over2 = find(-T_map > th_z); end
            if GInv,
                if GroupColorbars
                    index_over2 = [find(s_map > th_z); find(-s_map > th_z)];
                else
                    index_over2 = find(-s_map > th_z);
                end
            end
        end
    case 8
        if strcmp(F.tstr,'F')
            %if length(erdf) > 1 %?
            str_cor = [StatStr '_Anova'];
            xth_z = zeros(1,length(erdf));
            index_over = [];
            xerdf = erdf(:);
            udf = unique(xerdf);
            xth = zeros(1,length(udf));
            for i1=2:length(udf) %exclude 0
                if LKCflag %|| Z.UseCorrelRes
                    xth(i1) = calc_EC(LKC,p_value,tstr,[eidf,udf(i1)]);
                else
                    xth(i1) = spm_invFcdf(1-p_value, eidf, udf(i1));
                end
            end
            s_map = s_map(:); %??
            for i1 = 1:length(xerdf)
                if xerdf(i1) > 0
                    xth_z(i1) = xth(xerdf(i1)==udf);
                    if s_map(i1) > xth_z(i1)
                        index_over = [index_over i1];
                    end
                end
            end
            th_z = min(xth_z(xth_z>0));
            index_over2 = [];
        end        
    case 9
        if strcmp(F.tstr,'F')
            %if length(erdf) > 1 %?
            str_cor = [StatStr '_Anova'];
            xth_z = zeros(1,length(erdf));
            index_over = [];
            xerdf = erdf(:);
            udf = unique(xerdf);
            xth = zeros(1,length(udf));
            for i1=2:length(udf) %exclude 0
                xth(i1) = spm_invFcdf(1-p_value, eidf, udf(i1));                
            end
            s_map = s_map(:); %??
            for i1 = 1:length(xerdf)
                if xerdf(i1) > 0
                    xth_z(i1) = xth(xerdf(i1)==udf);
                    if s_map(i1) > xth_z(i1)
                        index_over = [index_over i1];
                    end
                end
            end
            th_z = min(xth_z(xth_z>0));
            index_over2 = [];
        end
end