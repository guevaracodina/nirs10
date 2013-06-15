function nirs_draw_figure_cine(F,W,Z)
try
    %NOTE: as the figures are generated in the Visible = Off mode, when opening
    %them in Matlab, one needs to set(gcf,'Visible','on') in order to see them at all!
    if isfield(W,'brain_view_mask_2d') % for back-compatibility
        F.s_map=F.s_map.*W.brain_view_mask_2d;
    end
    %Get threshold for maps
    [th_z str_cor index_over index_over2] = nirs_get_index_over_cine(F.erdf,F.s_map,Z.p_value);
    I = [];
    I.index_over = index_over;
%    Y1 = nirs_make_figure(I,F,W,Z,str_cor,th_z,0);
%     if ~isempty(Y1)
%         DF.fh1 = Y1.fh1;
%         DF.ax1 = Y1.ax1;
%         DF.hc1 = Y1.hc1;
%         DF.split1 = Y1.split1;
%         DF.hc1_min = Y1.hc1_min;
%         DF.hc1_max = Y1.hc1_max;
%         DF.tick_number = Y1.tick_number;
%         DF.fontsize_choice = Y1.fontsize_choice;
%         DF.GroupColorbars = Z.GroupColorbars;
%     else
%         DF = [];
%     end
        
    I.index_over2 = index_over2;
    Y2 = nirs_make_figure(I,F,W,Z,str_cor,th_z,1);
%     if ~isempty(Y2)
%         DF.fh2 = Y2.fh1;
%         DF.ax2 = Y2.ax1;
%         if Z.GroupColorbars
%             DF.hc2 = Y2.hc1;
%             DF.hc2_min = Y2.hc1_min;
%             DF.hc2_max = Y2.hc1_max;
%             DF.split2 = Y2.split1;
%         else
%             if isfield(Y2,'hc1')
%                 DF.hc1 = Y2.hc1;
%                 DF.hc1_min = Y2.hc1_min;
%                 DF.hc1_max = Y2.hc1_max;
%             end
%             if isfield(Y2,'hc2')
%                 DF.hc2 = Y2.hc2;
%                 DF.hc2_min = Y2.hc2_min;
%                 DF.hc2_max = Y2.hc2_max;
%             end
%             DF.split1 = Y2.split1;
%         end
%         DF.tick_number = Y2.tick_number;
%         DF.fontsize_choice = Y2.fontsize_choice;
%         DF.GroupColorbars = Z.GroupColorbars;
%     end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem drawing figure');
end