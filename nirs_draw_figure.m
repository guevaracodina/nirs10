function DF = nirs_draw_figure(fign,F,W,Z,G)
try
    if ~isstruct(G)
        LKC = G;
    else
        LKC = G.LKC;
    end
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
    if isfield(F,'sum_kappa')
        sum_kappa = F.sum_kappa;
    else
        sum_kappa = 0;
    end
    nchn = 0;
    if ~(fign == 4) && ~(fign == 5) && ~(fign == 6) && ~(fign == 7)
        switch F.hb
            case 'HbO'
                nchn = length(W.ch_HbO);
            case 'HbT'
                nchn = length(W.ch_HbT);
            case 'HbR'
                nchn = length(W.ch_HbR);
        end
    end
    if ~isfield(Z,'StatStr2')
        Z.StatStr2 = '';
        Z.UseCorrelRes = 1;
    end
    %Get threshold for maps
    [th_z str_cor index_over index_over2] = nirs_get_threshold(fign,F,tstr,erdf,...
        eidf,s_map,GInv,p_value,Z.StatStr,Z.StatStr2,Z.GroupColorbars,G,Z.UseCorrelRes,sum_kappa,nchn);
    
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