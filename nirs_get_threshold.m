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
        if LKCflag % || UseCorrelRes
            th_z = calc_EC(LKC,p_value,tstr,[eidf,erdf]);
%             ps = sort();
%             [P] = spm_P_FDR(Z,[eidf,erdf],tstr,p_value*,Ps)
            %Bonferroni
            th_z2 = spm_invTcdf(1-p_value/nchn, erdf);
            str_cor = StatStr;
            %Take lower threshold of EC or Bonferroni
            if th_z2 < th_z
                th_z = th_z2;
                str_cor = [StatStr 'Bonf'];                
            end           
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
                    switch StatStr
                        case 'EC'
                            test_LKC_by_pixel = 0;
                            if test_LKC_by_pixel && isfield(G,'ns')
                                th_zi = zeros(1,G.ns-G.min_subj+1);
                                for i0=1:(G.ns-G.min_subj+1)
                                    th_zi(i0) = calc_EC(G.LKCi{i0},p_value,tstr,[G.eidfi(i0),G.erdfi(i0)]);
                                end
                                th_z = max(th_zi(1:end/2)); %completely heuristic
                                %a=[]; for i0=1:(G.ns-G.min_subj+1), a = [a length(G.idxi{i0})]; end, a
                                str_cor = StatStr;
                            else
                                th_z = calc_EC(LKC,p_value,tstr,[eidf,erdf]);
                                %Bonferroni
                                th_z2 = spm_invTcdf(1-p_value/nchn, erdf);% KP ATTENTION! ERROR! nchn is always 0; warnings will be produced. Bonf is not working well 
                                str_cor = StatStr;
                                %Take lower threshold of EC or Bonferroni
                                if th_z2 < th_z
                                    th_z = th_z2;
                                    str_cor = [StatStr 'Bonf'];
                                end
                            end
                        case 'Bonf'
                        case '2DpFDR'
                            switch F.hb
                                case 'HbR'
                                    ch = 2;
                                case 'HbO'
                                    ch = 1;
                                case 'HbT'
                                    ch = 3;
                                otherwise
                                    ch = 0;
                            end
                            OP.fdrtype = 1; %FDR-BH
                            OP.u_thz = 2.5; %initial threshold
                            OP.ttype = 1; %one-tailed t-test
                            [p_thz t_thz fdrflag first_k stop_k] = nirs_2D_Bonferroni_FDR_threshold(3,s_map,ones(size(s_map)),ch,p_value,erdf,OP);
                            if fdrflag % Found a valid threhold
                                th_z = abs(t_thz);
                            else % All null hypothesis have to be accepted
                                disp('2D pFDR did not find a threshold. Infinite number will be imposed to produce a null image.');
                                th_z = Inf; % Give a high number so that a null image will be produced
                            end
                            str_cor = StatStr;
                        otherwise
                    end
                else
                    th_z = spm_invTcdf(1-p_value, erdf);
                end
                %index_over = find(s_map > th_z);
                index_over = find(s_map >= th_z); %KP changed
                %if GInv, index_over2 = find(-T_map > th_z); end
                if GInv,
                    if GroupColorbars
                        %index_over2 = [find(s_map > th_z); find(-s_map > th_z)];
                        index_over2 = [find(s_map >= th_z); find(-s_map >= th_z)]; %KP changed
                    else
                        %index_over2 = find(-s_map > th_z);
                        index_over2 = find(-s_map >= th_z); %KP changed
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
            if isnan(th_z)                
                th_z = spm_invTcdf(1-p_value, 10);
                disp('Problem with erdf: Imposing threshold')
            end
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
                    min_udf = 2;
                    if udf(i1) > min_udf
                        xth(i1) = calc_EC(LKC,p_value,tstr,[eidf,udf(i1)]);
                    else
                        xth(i1) =  spm_invFcdf(1-p_value, eidf, udf(i1));
                    end
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
            %take median or max of udf values??? why do we have different
            %udf values??? problem with mask???
            th_z = min(xth_z(xth_z> spm_invFcdf(1-p_value, eidf, max(udf)))); %exceed threshold of unc. F stat
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