function H = initialize_assembled_figures(Z,H,f1,label)
%All this does is save a previously constructed arrangement of figures
if f1 > 0
    strf = ['_S' int2str(f1)];
    %str_grp = '';
else
    strf = '';
    %str_grp = 'Group';
end

if isfield(Z,'scon')
    scon = [Z.scon '_'];
else
    scon = '';
end
if isfield(Z,'Idx')
    subj_str = ['_S' gen_num_str(Z.Idx,2)];
    try 
    if isfield(Z,'subj_id')
        %This is useful if the subject numbers are not preserved due to
        %some subjects being dropped for various reasons
        subj_str = [Z.subj_id subj_str]; 
    end
    end
else
    subj_str = ''; %group analysis
end
if isfield(Z,'spec_hemi')
filestr = [subj_str '_' num2str(Z.p_value) '_' Z.spec_hemi strf scon];
else
filestr = [subj_str '_' num2str(Z.p_value) strf scon];   
end

vis = Z.cbar.visible;

if Z.GFIS
    if Z.write_neg_pos
        H.Pt = figure('Visible',vis,'Name',[label '_' Z.StatStr '_' filestr '_Pos'],'NumberTitle','off');
        if Z.output_unc
            H.Pu = figure('Visible',vis,'Name',[label '_unc_' filestr '_Pos'],'NumberTitle','off');
        end
    end
    if Z.GInv
        if Z.write_neg_pos
            H.Nt = figure('Visible',vis,'Name',[label '_' Z.StatStr '_' filestr '_Neg'],'NumberTitle','off');
        end
        H.Ct = figure('Visible',vis,'Name',[label '_' Z.StatStr '_' filestr],'NumberTitle','off');
        if Z.output_unc
            if Z.write_neg_pos
                H.Nu = figure('Visible',vis,'Name',[label '_unc_' filestr '_Neg'],'NumberTitle','off');
            end
            H.Cu = figure('Visible',vis,'Name',[label '_unc_' filestr],'NumberTitle','off');
        end
    end
end