function [TOPO H] = call_figure_2anova(TOPO,H,Z,W,F,A,CF,v1,h1,hb,Astr,z1,LKC)
try
    if z1 <= 3
        TOPO.v{v1}.group.hb{h1}.c{2*z1-1} = A;
    else
        %not stored for now
    end
    filestr = [num2str(Z.p_value) '_' W.spec_hemi '_' hb '_' Astr];
    filestr_fig = [num2str(Z.p_value) ' ' W.spec_hemi ' ' hb ' ' Astr];
    info1 = filestr;
    info_for_fig1 = filestr_fig;
    F.contrast_info = info1;
    F.contrast_info_for_fig = info_for_fig1;
    F.contrast_info_both = filestr; %same for Pos and Neg, used for combined figures
    F.contrast_info_both_for_fig = filestr_fig; %same for Pos and Neg, used for combined figures
    F.s_map = A.Tmap;
    F.erdf = A.erdf;
    F.eidf = A.eidf;
    F.tstr = 'F'; %tstr;
    F.hb = hb;
    if Z.LKC
        DF = nirs_draw_figure(5,F,W,Z,LKC);
        H = nirs_copy_figure(H,DF,CF,1,hb,1,F.tstr,1,0);
    end
    if Z.output_unc
        DF = nirs_draw_figure(7,F,W,Z,[]);
        H = nirs_copy_figure(H,DF,CF,1,hb,1,F.tstr,0,0);
    end
catch exception2
    disp(exception2.identifier);
    disp(exception2.stack(1));
end