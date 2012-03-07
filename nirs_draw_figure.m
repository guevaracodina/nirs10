function DF = nirs_draw_figure(fign,F,W,Z,LKC)
%SORRY!!! THIS IS VERY POORLY CODED UP NOW - too many options and cases
%code is now very obscure and hard to maintain
try
    DF = [];
    %NOTE: as the figures are generated in the Visible = Off mode, when opening
    %them in Matlab, one needs to set(gcf,'Visible','on') in order to see them
    %at all!
    GInv = Z.GInv;
    p_value = Z.p_value;
    erdf = F.erdf;
    eidf = F.eidf;
    tstr = F.tstr;
    if isfield(W,'brain_view_mask_2d') % for back-compatibility
        F.s_map=F.s_map.*W.brain_view_mask_2d;
    end
    s_map = F.s_map;
    try
        sum_kappa = F.sum_kappa;
    end
    if ~(fign == 4) && ~(fign == 5) && ~(fign == 6)
        switch F.hb
            case 'HbO'
                nchn = length(W.ch_HbO);
            case 'HbT'
                nchn = length(W.ch_HbT);
            case 'HbR'
                nchn = length(W.ch_HbR);
        end
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
                if Z.GroupColorbars
                    index_over2 = [find(s_map > th_z); find(-s_map > th_z)];
                else
                    index_over2 = find(-s_map > th_z);
                end
            end
        case 3
            %%%%%%%%%%%%%%%%%
            %Get threshold
            if Z.LKC || Z.UseCorrelRes
                th_z = calc_EC(LKC,p_value,tstr,[eidf,erdf]);
                str_cor = Z.StatStr;
            else
                if tstr == 'T'
                    %Tube
                    str_cor = Z.StatStr;
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
                    str_cor = Z.StatStr2;
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
                    if Z.GroupColorbars
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
            
            if strcmp(F.tstr,'F')
                %if length(erdf) > 1 %?
                str_cor = 'Anova';
                xth_z = zeros(1,length(erdf));
                index_over = [];
                xerdf = erdf(:);
                udf = unique(xerdf);
                xth = zeros(1,length(udf));
                for i1=2:length(udf) %exclude 0
                    if Z.LKC || Z.UseCorrelRes
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
                str_cor = 'Group';
                if erdf == 0
                    disp([F.contrast_info ' : Problem with number of degrees of freedom']);
                    index_over = [];
                    if GInv
                        index_over2 = [];
                    end
                    th_z = 0;
                else
                    if Z.LKC  
                        th_z = calc_EC(LKC,p_value,tstr,[eidf,erdf]);
                        str_cor = Z.StatStr;
                    else
                        th_z = spm_invTcdf(1-p_value, erdf);
                    end
                    index_over = find(s_map > th_z);
                    %if GInv, index_over2 = find(-T_map > th_z); end
                    if GInv,
                        if Z.GroupColorbars
                            index_over2 = [find(s_map > th_z); find(-s_map > th_z)];
                        else
                            index_over2 = find(-s_map > th_z);
                        end
                    end
                end
            end
        case 5
            %for 2-anova
            str_cor = ['2A' int2str(eidf) '_' int2str(erdf-eidf)];
            if strcmp(F.tstr,'F') % should always be the case
                if Z.LKC || Z.UseCorrelRes
                    th_z = calc_EC(LKC,p_value,tstr,[eidf,erdf]);
                else
                    th_z = spm_invFcdf(1-p_value, eidf, erdf);
                end
                index_over = find(s_map > th_z);
                index_over2 = []; %not used
            end
        case 6 %group, uncorrected
            str_cor = 'Group_unc';
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
                    if Z.GroupColorbars
                        index_over2 = [find(s_map > th_z); find(-s_map > th_z)];
                    else
                        index_over2 = find(-s_map > th_z);
                    end
                end
            end
    end
    I = [];
    I.index_over = index_over;
    Y1 = nirs_make_figure(I,F,W,Z,str_cor,th_z,0);
    if ~isempty(Y1)
        DF.fh1 = Y1.fh1;
        DF.ax1 = Y1.ax1;
        DF.hc1 = Y1.hc1;
        DF.split1 = Y1.split1;
        DF.hc1_min = Y1.hc1_min;
        DF.hc1_max = Y1.hc1_max;
        DF.tick_number = Y1.tick_number;
        DF.fontsize_choice = Y1.fontsize_choice;
        DF.GroupColorbars = Z.GroupColorbars;
    else
        DF = [];
    end
    
    if strfind(F.contrast_info,'Pos')
        run_other_chart = 1;
    else
        run_other_chart = 0;
    end
    if (GInv && Z.GroupColorbars)||(~Z.GroupColorbars && run_other_chart)
        if tstr == 'T'
            I.index_over2 = index_over2;
            Y2 = nirs_make_figure(I,F,W,Z,str_cor,th_z,1);
            if ~isempty(Y2)
                DF.fh2 = Y2.fh1;
                DF.ax2 = Y2.ax1;
                if Z.GroupColorbars
                    DF.hc2 = Y2.hc1;
                    DF.hc2_min = Y2.hc1_min;
                    DF.hc2_max = Y2.hc1_max;
                    DF.split2 = Y2.split1;
                else
                    if isfield(Y2,'hc1')
                        DF.hc1 = Y2.hc1;
                        DF.hc1_min = Y2.hc1_min;
                        DF.hc1_max = Y2.hc1_max;
                    end
                    if isfield(Y2,'hc2')
                        DF.hc2 = Y2.hc2;
                        DF.hc2_min = Y2.hc2_min;
                        DF.hc2_max = Y2.hc2_max;
                    end
                    DF.split1 = Y2.split1;
                end
                DF.tick_number = Y2.tick_number;
                DF.fontsize_choice = Y2.fontsize_choice;
                DF.GroupColorbars = Z.GroupColorbars;
            end
        else %'F' - careful - as this does not create a new figure, need to not close the figure
            %later on
            DF.fh2 = Y1.fh1;
            DF.ax2 = Y1.ax1;
            if isfield(Y1,'hc1')
                DF.hc2 = Y1.hc1;
                DF.hc2_min = Y1.hc1_min;
                DF.hc2_max = Y1.hc1_max;
            end
            DF.split2 = Y1.split1;
            DF.tick_number = Y1.tick_number;
            DF.fontsize_choice = Y1.fontsize_choice;
            DF.GroupColorbars = Z.GroupColorbars;
        end
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem drawing figure');
end