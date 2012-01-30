function H = initialize_assembled_figures(Z,H,f1,label)
vis = Z.cbar.visible;
p = num2str(Z.p_value);
if f1 == 0
    str = ''; %group of sessions
else
    str = ['_S' int2str(f1)]; %individual sessions
end
if Z.GFIS
    if Z.write_neg_pos
        H.Pt = figure('Visible',vis,'Name',[label '_' Z.StatStr '_' p str '_Pos'],'NumberTitle','off');
        if Z.output_unc
            H.Pu = figure('Visible',vis,'Name',[label '_unc_' p str '_Pos'],'NumberTitle','off');
        end
    end
    if Z.GInv
        if Z.write_neg_pos
            H.Nt = figure('Visible',vis,'Name',[label '_' Z.StatStr '_' p str '_Neg'],'NumberTitle','off');
        end
        H.Ct = figure('Visible',vis,'Name',[label '_' Z.StatStr '_' p str],'NumberTitle','off');
        if Z.output_unc
            if Z.write_neg_pos
                H.Nu = figure('Visible',vis,'Name',[label '_unc_' p str '_Neg'],'NumberTitle','off');
            end
            H.Cu = figure('Visible',vis,'Name',[label '_unc_' p str],'NumberTitle','off');
        end
    end
end