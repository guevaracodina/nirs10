function [TOPO H] = constrasts_core(Z,W,TOPO,xX,f1)
%Note that even for contrasts pertaining to individual
%sessions, the number of degrees of freedom it taken
%to be that for all the sessions - this is what is
%done in SPM - this looks like an error, however, it
%should be small under situations of interest, as the
%number of dof in each session is fairly large, even
%more so for fNIRS than for fMRI
if f1 ==0
    xCon = TOPO.xCon;
else
    xCon = TOPO.SSxCon{f1};
end
%Initialize figures
%Handles for assembled figures
H = initialize_assembled_figure_handles;
H = initialize_assembled_figures(Z,H,f1,'Z');
for h1=1:3
    %HbO
    try
        hb = get_chromophore(h1);
        %Note that var(ch_HbO) depends on HbO vs HbR
        Q = interpolation_kernel(h1,W,Z.LKC,Z.UseCorrelRes);
        if Z.LKC || Z.UseCorrelRes
            C = loop_contrasts_new(h1,Q,W,xCon,Z,f1,TOPO.cSigma,xX);
        else
            C = loop_contrasts(h1,Q,W,xCon,Z);
        end
        if Z.gen_fig || Z.gen_tiff
            if Z.LKC || Z.UseCorrelRes
                H = interpolated_maps_new(Z,W,C,Q,xCon,f1,xX.erdf,hb,H);
            else
                H = interpolated_maps(Z,W,C,Q,xCon,f1,xX.erdf,hb,H);
            end
        end
        %fill TOPO
        TOPO = fill_TOPO(TOPO,C,W.side_hemi,f1,hb);
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Problem with ' hb ' contrast']);
    end
end
call_save_assembled_figures(Z,W,H,f1);
