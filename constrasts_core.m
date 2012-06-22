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

%**************************************************************************
%Ke Peng
%Do not understand why here it is an "isfield" function. As the Z.Avg is
%defined previously at Line 4 (prepare_constrast_core_call.m which is pre-executed than the current function),
%isfield(Z,'Avg') will always be 1 so that the Z.LKC will be always set to
%0. We are not able to perform an Euler Correction.
%Try to change it to "if Z.Avg == 1"
%18/06/2012
%**************************************************************************
if Z.Avg == 1
%if isfield(Z,'Avg')
%**************************************************************************
    W.Avg = Z.Avg;
if Z.Avg == 1
    Z.UseCorrelRes = 0;
    Z.LKC = 0;
    xX.erdf = 1;
    Z.output_unc = 1;
end
else
    W.Avg = 0;
end
%Initialize figures
%Handles for assembled figures
try
    Z.use_nCloop = 1;
    if Z.use_nCloop
        %maximal number of contrasts to group in assembled
        %figures
        nC = length(xCon);
        Z.nCloop = 4;
        nCl = 0;
        for c1=1:nC %Loop over contrasts
            Z.c1eff = mod(c1,Z.nCloop);
            if Z.c1eff == 0
                Z.c1eff = Z.nCloop;
            end
            if Z.c1eff == 1
                nCl = nCl + 1;
                Z.scon = '';
                for k0=c1:(c1+Z.nCloop-1)
                    if k0 <= nC
                        Z.scon = [Z.scon '_' xCon(k0).name];
                    end
                end
                %Handles for assembled figures
                H{nCl} = initialize_assembled_figure_handles;
                H{nCl} = initialize_assembled_figures(Z,H{nCl},f1,'Z');
            end
        end
    else
        H{1} = initialize_assembled_figure_handles;
        H{1} = initialize_assembled_figures(Z,H{1},f1,'Z');
    end
    for h1=1:3
        %HbO
        try
            hb = get_chromophore(h1);
            %Note that var(ch_HbO) depends on HbO vs HbR
            Q = interpolation_kernel(h1,W,Z.LKC,Z.UseCorrelRes);
            if Z.LKC || Z.UseCorrelRes || Z.Avg
                C = loop_contrasts_new(h1,Q,W,xCon,Z,f1,TOPO.cSigma,xX);
            else
                C = loop_contrasts(h1,Q,W,xCon,Z);
            end
            if Z.gen_fig || Z.gen_tiff
                if Z.LKC || Z.UseCorrelRes || Z.Avg
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
    if Z.use_nCloop
        nCl = 0;
        for c1=1:nC
            Z.c1eff = mod(c1,Z.nCloop);
            if Z.c1eff == 0
                Z.c1eff = Z.nCloop;
            end
            if Z.c1eff == Z.nCloop || c1 == nC
                nCl = nCl + 1;
                Z.scon = '';
                for k0=c1:(c1+Z.nCloop-1)
                    if k0 <= nC
                        Z.scon = [Z.scon xCon(k0).name];
                    end
                end 
                call_save_assembled_figures(Z,W,H{nCl},f1);
            end           
        end
    else
        call_save_assembled_figures(Z,W,H{1},f1);
    end
catch exception2
    disp(exception2.identifier);
    disp(exception2.stack(1));
    disp(['Problem with ' hb ' contrast -- setting up figures or saving them']);
end
