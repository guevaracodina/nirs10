function h = nirs_draw_figure_cine(F,W,Z)
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
    I.index_over2 = index_over2;
    h = nirs_make_figure_cine(I,F,W,Z,str_cor,th_z);
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem drawing figure');
end