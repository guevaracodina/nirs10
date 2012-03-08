function [H,TOPO,big_TOPO] = fill_group(H,TOPO,big_TOPO,v1,c1,h1,Z,W,F,CF,xCon,ns,shb)
try
    Sess = TOPO.Sess;
    Cp = TOPO.Cp;
    hb = get_chromophore(h1);
    new_version = 0;
    if shb
        strA = 'Pos';
        strB = 'Positive';
        sign_hb = 1;
    else
        strA = 'Neg';
        strB = 'Negative';
        sign_hb = -1;
    end
    if xCon(c1).STAT == 'T'
        fc = 0; %used only for FFX || nS==1
        %fill in cbeta and ccov_beta
        for f1=1:ns
            if Z.FFX || Z.nS==1
                %select sessions which had the contrast
                if any(f1==Sess{c1})
                    fc = fc+1;
                    if isfield(TOPO.v{v1}.s{f1}.hb{h1},'stat_map')
                        tmp = sign_hb*squeeze(TOPO.v{v1}.s{f1}.hb{h1}.beta_map(Cp{c1,f1},:,:));
                        cbeta(fc,:) = tmp(:);
                        if ~Z.simple_sum
                            tmp = squeeze(TOPO.v{v1}.s{f1}.hb{h1}.beta_map(Cp{c1,f1},:,:))./squeeze(TOPO.v{v1}.s{f1}.hb{h1}.stat_map(Cp{c1,f1},:,:));
                        end
                        ccov_beta(fc,:) = tmp(:).^2;
                        new_version = 1;
                    else
                        %now use Cp{c1,f1} to access the required c_interp_beta instead of c1
                        tmp = sign_hb*squeeze(TOPO.v{v1}.s{f1}.hb{h1}.c_interp_beta(Cp{c1,f1},:,:));
                        cbeta(fc,:) = tmp(:);
                        tmp = squeeze(TOPO.v{v1}.s{f1}.hb{h1}.c_cov_interp_beta(Cp{c1,f1},:,:));
                        ccov_beta(fc,:) = tmp(:);
                    end
                end
            else
                if ~isfield(big_TOPO{f1}.v{v1},'s')
                    if isfield(big_TOPO{f1}.v{v1}.g{1}.hb{h1},'stat_map')
                        tmp = sign_hb*squeeze(big_TOPO{f1}.v{v1}.g{1}.hb{h1}.beta_map(c1,:,:));
                        cbeta(f1,:) = tmp(:);
                        if ~Z.simple_sum
                            tmp = squeeze(big_TOPO{f1}.v{v1}.g{1}.hb{h1}.beta_map(c1,:,:))./squeeze(big_TOPO{f1}.v{v1}.g{1}.hb{h1}.stat_map(c1,:,:));
                        end
                        ccov_beta(f1,:) = tmp(:).^2;
                        new_version = 1;
                    else
                        %group analysis of a group of sessions analysis
                        tmp = sign_hb*squeeze(big_TOPO{f1}.v{v1}.g{1}.hb{h1}.c_interp_beta(c1,:,:));
                        cbeta(f1,:) = tmp(:);
                        tmp = squeeze(big_TOPO{f1}.v{v1}.g{1}.hb{h1}.c_cov_interp_beta(c1,:,:));
                        ccov_beta(f1,:) = tmp(:);
                    end
                else
                    %for is1=1:length(big_TOPO{f1}.v{v1}.s)
                    is1 = Z.group_session_to_average;
                    %do each session separately
                    if isfield(big_TOPO{f1}.v{v1}.s{is1}.hb{h1},'stat_map')
                        tmp = sign_hb*squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.beta_map(c1,:,:));
                        cbeta(f1,:) = tmp(:);
                        if ~Z.simple_sum
                            tmp = squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.beta_map(c1,:,:))./squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.stat_map(c1,:,:));
                        end
                        ccov_beta(f1,:) = tmp(:).^2;
                        new_version = 1;
                    else
                        tmp = sign_hb*squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_interp_beta(c1,:,:));
                        cbeta(f1,:) = tmp(:);
                        tmp = squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_cov_interp_beta(c1,:,:));
                        ccov_beta(f1,:) = tmp(:);
                    end
                end
            end
        end
    else %quick fix for F stats
        fc = 0; %used only for FFX || nS==1
        for f1=1:ns
            if Z.FFX || Z.nS==1
                %select sessions which had the contrast
                if any(f1==Sess{c1})
                    fc = fc+1;
                    tmp = squeeze(TOPO.v{v1}.s{f1}.hb{h1}.c_interp_F(Cp{c1,f1},:,:));
                    cbeta(fc,:) = tmp(:);
                end
            else
                if ~isfield(big_TOPO{f1}.v{v1},'s')
                    %group analysis of a group of sessions analysis
                    tmp = squeeze(big_TOPO{f1}.v{v1}.g{1}.hb{h1}.c_interp_F(c1,:,:));
                    cbeta(f1,:) = tmp(:);
                else
                    %do each session separately
                    is1 = Z.group_session_to_average;
                    %for is1=1:length(big_TOPO{f1}.v{v1}.s)
                    tmp = squeeze(big_TOPO{f1}.v{v1}.s{is1}.hb{h1}.c_interp_F(c1,:,:));
                    cbeta(f1,:) = tmp(:);
                    %end
                end
            end
        end
    end
    %Generate group result as t-stat
    try
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
    catch exception2
        disp(exception2.identifier);
        disp(exception2.stack(1));
    end
    if Z.FFX || Z.nS==1
        fg = 'g';
    else
        fg = 'group';
    end
    
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
    filestr = [num2str(Z.p_value) '_' W.spec_hemi '_' hb];
    filestr_fig = [num2str(Z.p_value) ' ' W.spec_hemi ' ' hb];
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
    %     if strcmp(F.tstr,'T')
    %         F.con = G.beta_group;
    %         F.ess = [];
    %     else
    %         F.con = [];
    %         F.ess = G.beta_group;
    %     end
    DF = [];
    if Z.output_unc
        try
            DF = nirs_draw_figure(6,F,W,Z,G.LKC);
        catch exception2
            disp(exception2.identifier);
            disp(exception2.stack(1));
        end
        try
            if Z.GFIS, H = nirs_copy_figure(H,DF,CF,c1,hb,shb,F.tstr,1-Z.output_unc,Z.write_neg_pos); end
        catch exception2
            disp(exception2.identifier);
            disp(exception2.stack(1));
        end
    end
    DF = [];
    try
        DF = nirs_draw_figure(4,F,W,Z,G.LKC);
    catch exception2
        disp(exception2.identifier);
        disp(exception2.stack(1));
    end
    try
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