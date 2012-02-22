function call_save_assembled_figures(Z,W,H,f1)
pos_label = 'pos';
neg_label = 'neg';
unc_label = 'unc';
try %works for both group and single sessions
    if Z.GFIS
        if Z.write_neg_pos || ~Z.GInv
            %Just the positive images
            save_assembled_figures(Z,W,H.Pt,pos_label,Z.StatStr,f1);
            if Z.output_unc
                save_assembled_figures(Z,W,H.Pu,pos_label,unc_label,f1);
            end
        end
        if Z.GInv
            if Z.write_neg_pos
                %Just the negative images
                save_assembled_figures(Z,W,H.Nt,neg_label,Z.StatStr,f1);
                if Z.output_unc
                    save_assembled_figures(Z,W,H.Nu,neg_label,unc_label,f1);
                end
            end
            %These are usually the only figures produced (FWE corrected and
            %showing both increases and decreases)
            save_assembled_figures(Z,W,H.Ct,'',Z.StatStr,f1);
            if Z.output_unc %uncorrected images
                save_assembled_figures(Z,W,H.Cu,'',unc_label,f1);
            end
        end
    end
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem assembling figures');
end