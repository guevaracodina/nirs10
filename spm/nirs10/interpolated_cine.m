function interpolated_cine(Z,W,C,Q,f1,erdf,hb,u1,k0)
try
    stat_map = C.stat_map;
    index_mask = Q.index_mask;
    s1 = W.s1;
    s2 = W.s2;
    spec_hemi = W.spec_hemi;
    p_value = Z.p_value;
    nD = size(stat_map,1);
    
    %Structure F
    load Split
    F.split = split;
    F.pathn = Z.dir1;
    F.erdf = erdf;
    F.hb = hb;
    strf = ['_S' int2str(f1)];
    strf_fig = [' S' int2str(f1)];
    filestr = [num2str(p_value) '_' spec_hemi '_' hb strf];
    filestr_fig = [num2str(p_value) ' ' spec_hemi ' ' hb strf_fig];
    %loop over onset delays   
    for c1=1:nD
        s_map = zeros(s1, s2);
        tstr = 'T';
        tname = ['t' int2str(u1) 'o' gen_num_str(k0,3) 't' gen_num_str(c1,2)]; %onset delay
        F.tstr = tstr;
        s_map(index_mask) = squeeze(stat_map(c1,index_mask));
        F.s_map = s_map;
%        F.eidf = 1;
%        F.contrast_info = [filestr tname];
%        F.contrast_info_for_fig = [filestr_fig tname];
        F.contrast_info_both = [filestr tname]; %same for Pos and Neg, used for combined figures
        F.contrast_info_both_for_fig = [filestr_fig tname]; %same for Pos and Neg, used for combined figures
        nirs_draw_figure_cine(F,W,Z); %DF: draw figure structure: handles to figure, axes and colorbar
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem in interpolated_maps')
end