function [H,TOPO,big_TOPO] = fill_group(H,TOPO,big_TOPO,v1,c1,h1,Z,W,F,CF,xCon,ns,shb)
try
    hb = get_chromophore(h1);
    if shb
        strA = 'Pos';
        strB = 'Positive';
    else
        strA = 'Neg';
        strB = 'Negative';
    end
    %Get mask
    BV = [];
    if Z.FFX || Z.nS==1
        if isfield(TOPO.rendered_MNI{v1},'view_mask_2d')
            BV.brain_view_mask_2d{1} = TOPO.rendered_MNI{v1}.view_mask_2d;
        end
    else
        if isfield(big_TOPO{1}.rendered_MNI{v1},'view_mask_2d')
            for i0=1:length(big_TOPO)
                BV.brain_view_mask_2d{i0} = big_TOPO{i0}.rendered_MNI{v1}.view_mask_2d;
            end
        end
    end
    [cbeta ccov_beta new_version] = fill_group_arrays(TOPO,big_TOPO,v1,c1,h1,Z,xCon,ns,shb);
    if ~isempty(BV)
        if length(BV.brain_view_mask_2d) == 1
            cbeta = cbeta.*repmat(BV.brain_view_mask_2d{1}(:)',[size(cbeta,1) 1]);
        else
            for i0=1:size(cbeta,1)
                cbeta(i0,:) = cbeta(i0,:).*BV.brain_view_mask_2d{i0}(:)';
            end
        end
    end

    %Generate group result as t-stat
    try
        %This is a mess, with confusing storage of LKC0 and interplay with
        %Z.simple_sum
        if Z.two_samples
            G = liom_2samples_t_test(Z,cbeta,W.s1,W.s2);
        else
            if (Z.LKC || new_version)
                G = liom_group_new(cbeta,W.s1,W.s2);
                LKC0 = G.LKC; %store
            else
                if strcmp(xCon(c1).STAT,'T')
                    G = liom_group(cbeta,ccov_beta,W.s1,W.s2,Z.min_s,Z.FFX,Z.simple_sum);
                else
                    G = liom_group_F(cbeta,W.s1,W.s2);
                end
                G.LKC = [];
            end
            if (Z.LKC || new_version) && ~Z.simple_sum
                if strcmp(xCon(c1).STAT,'T')
                    G = liom_group(cbeta,ccov_beta,W.s1,W.s2,Z.min_s,Z.FFX,Z.simple_sum);
                else
                    G = liom_group_F(cbeta,W.s1,W.s2);
                end
                G.LKC = LKC0;
            end
        end
    catch exception2
        disp(exception2.identifier);
        disp(exception2.stack(1));
    end
    if Z.FFX || Z.nS==1
        fg = 'g';
    else
        fg = 'group';
    end
    if Z.two_samples
        str2t = '_2stt';
        str2t2 = ' 2stt';
    else
        str2t = '';
        str2t2 = '';
    end
    %Fill TOPO -- needs to be cleaned up
    %factor of 2 in c{} was for positive vs negative contrasts -- may want
    %to get rid of that -- but the extract_map module should be updated too
    TOPO.v{v1}.(fg).hb{h1}.c{2*c1-shb}.Tmap = G.tmap_group;
    TOPO.v{v1}.(fg).hb{h1}.c{2*c1-shb}.erdf = G.erdf_group;
    TOPO.v{v1}.(fg).hb{h1}.c{2*c1-shb}.beta_group = G.beta_group;
    TOPO.v{v1}.(fg).hb{h1}.c{2*c1-shb}.std_group = G.std_group;
    TOPO.v{v1}.(fg).hb{h1}.c{2*c1-shb}.type = strB;
    TOPO.v{v1}.(fg).hb{h1}.c{2*c1-shb}.LKC = G.LKC;
    TOPO.v{v1}.(fg).hb{h1}.c{2*c1-shb}.var_bs = G.var_bs;
    TOPO.v{v1}.(fg).hb{h1}.c{2*c1-shb}.c = xCon(c1);
    erdf_group = max(G.erdf_group(:)); %quick fix...
    F.s_map = G.tmap_group;
%     %******************Ke Peng, for one-tailed t-test, for display only
%     switch hb
%         case {'HbO','HbT'}
%             F.s_map(find(F.s_map < 0)) = 0;
%         case 'HbR'
%             F.s_map(find(F.s_map > 0)) = 0;
%         otherwise
%     end
%   *****
    TOPO.v{v1}.(fg).hb{h1}.c{2*c1-shb}.Tmap = F.s_map;
    filestr = [num2str(Z.p_value) '_' W.spec_hemi '_' hb str2t];
    filestr_fig = [num2str(Z.p_value) ' ' W.spec_hemi ' ' hb str2t2];
    info1 = [filestr '_' strA xCon(c1).name];
    info_for_fig1 = [filestr_fig ' ' strA xCon(c1).name];
    F.contrast_info = info1;
    F.contrast_info_for_fig = info_for_fig1;
    F.contrast_info_both = [filestr xCon(c1).name]; %same for Pos and Neg, used for combined figures
    F.contrast_info_both_for_fig = [filestr_fig xCon(c1).name]; %same for Pos and Neg, used for combined figures
    
    F.erdf = erdf_group;
    F.eidf = xCon(c1).eidf;
    F.tstr = 'T'; %xCon(c1).STAT; %tstr;
    F.hb = hb;
    if Z.output_unc
        try
            DF = nirs_draw_figure(6,F,W,Z,G);
            if Z.GFIS, H = nirs_copy_figure(H,DF,CF,c1,hb,shb,F.tstr,1-Z.output_unc,Z.write_neg_pos); end
        catch exception2
            disp(exception2.identifier);
            disp(exception2.stack(1));
        end
    end
    try
        DF = nirs_draw_figure(4,F,W,Z,G);
        if ~isfield(TOPO.v{v1}.(fg).hb{h1}.c{2*c1-shb},'th_z') || (isfield(TOPO.v{v1}.(fg).hb{h1}.c{2*c1-shb},'th_z') && isinf(TOPO.v{v1}.(fg).hb{h1}.c{2*c1-shb}.th_z))
            try
                TOPO.v{v1}.(fg).hb{h1}.c{2*c1-shb}.th_z = DF.th_z; % KP Record threshold value
            catch
                TOPO.v{v1}.(fg).hb{h1}.c{2*c1-shb}.th_z = Inf;
            end
        end

        if Z.GFIS, H = nirs_copy_figure(H,DF,CF,c1,hb,shb,F.tstr,Z.LKC,Z.write_neg_pos); end
    catch exception2
        disp(exception2.identifier);
        disp(exception2.stack(1));
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp(['Group: Problem with a specific contrast ' int2str(c1) ' and chromophore ' hb ' for view ' int2str(v1)]);
end