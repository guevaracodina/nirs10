function [coreg xSPM opt] = nirs_orthPrep(coreg, sel, render_template, coreg_layer, m2v)
%Ke Peng, 2013-09-12, called in nirs_run_liom_orth_coreg

ftemplate = fullfile(fileparts(which('spm')),'toolbox','nirs10','xSPMtemplate.mat');
load(ftemplate);
xSPM.Z = [];
xSPM.XYZ = [];
xSPM.XYZmm = [];
xSPM.label = {};
coreg_projected_channels_on_cortex_rather_than_midpoint = 1;

%Add fiducial points
xSPM.Z = [15 15 15]; %Fiducial value
xSPM.XYZmm = coreg.corr.F.r.m.mm.p;
xSPM.label = {'F-N','F-L','F-R'};

Ns = coreg.corr.S.N;
Nd = coreg.corr.D.N;
Nc = size(coreg.corr.C.id,2)/2;
Q = coreg.Q;

%If no channel coordinates, add them
if render_template
    try
        Ch_wc1 = coreg.corr.C.w.m.mm.c1.p(1:3,:);
        Ch_wS = coreg.corr.C.w.m.mm.fp(1:3,:);
        opt.flag = 0;
    catch
        Ch_wS = zeros(4, Nc);
        if coreg_projected_channels_on_cortex_rather_than_midpoint
            for i0 = 1:Nc
                Si0 = coreg.corr.C.id(2,i0);
                Di0 = coreg.corr.C.id(3,i0)+Ns;
                Ch_wS(:,i0) = (coreg.corr.P.w.m.mm.p(:,Si0) + coreg.corr.P.w.m.mm.p(:,Di0))/2;
            end
            Ch_wc1 = projection_CS(Ch_wS); %Korean's method
        else
            Ch_wc1 = zeros(4, Nc);            
            for i0 = 1 : Nc
                Si0 = coreg.corr.C.id(2,i0);
                Di0 = coreg.corr.C.id(3,i0)+Ns;
                Ch_ws(:,i0) = (coreg.corr.P.w.m.mm.p(:,Si0) + coreg.corr.P.w.m.mm.p(:,Di0))/2;
                Ch_wc1(:,i0) = (coreg.corr.P.w.m.mm.c1.p(:,Si0) + coreg.corr.P.w.m.mm.c1.p(:,Di0))/2;
            end
        end
        coreg.corr.C.w.m.mm.fp = Ch_wS;
        coreg.corr.C.w.m.mm.c1.p = Ch_wc1;
        coreg.corr.C.w.m.vx.fp = m2v.w\Ch_wS;
        coreg.corr.C.w.m.vx.c1.p = m2v.w\Ch_wc1;
        coreg.corr.C.r.m.mm.fp = Q\Ch_wS;
        coreg.corr.C.r.m.mm.c1.p = Q\Ch_wc1;
        coreg.corr.C.r.m.vx.fp = m2v.m\coreg.corr.C.r.m.mm.fp;
        coreg.corr.C.r.m.vx.c1.p = m2v.m\coreg.corr.C.r.m.mm.c1.p;
        opt.flag = 1;
    end
else
    try
        Ch_mc1 = coreg.corr.C.r.m.mm.c1.p(1:3,:);
        Ch_mS = coreg.corr.C.r.m.mm.fp(1:3,:);
        opt.flag = 0;
    catch
        Ch_mS = zeros(3, Nc);
        if coreg_projected_channels_on_cortex_rather_than_midpoint
            for i0 = 1:Nc
                Si0 = coreg.corr.C.id(2,i0);
                Di0 = coreg.corr.C.id(3,i0)+Ns;
                Ch_mS(:,i0) = (coreg.corr.P.r.m.mm.p(:,Si0) + coreg.corr.P.r.m.mm.p(:,Di0))/2;
            end
            Ch_wS = Q * [Ch_mS; ones(1,size(Ch_mS,2))];
            Ch_wc1 = projection_CS(Ch_wS); %Korean's method
            Ch_mc1 = Q\Ch_wc1;
            Ch_mc1 = Ch_mc1(1:3,:);
        else
            Ch_mc1 = zeros(3, Nc);
            for i0 = 1 : Nc
                Si0 = coreg.corr.C.id(2,i0);
                Di0 = coreg.corr.C.id(3,i0)+Ns;
                Ch_mS(:,i0) = (coreg.corr.P.r.m.mm.p(:,Si0) + coreg.corr.P.r.m.mm.p(:,Di0))/2;
                Ch_mc1(:,i0) = (coreg.corr.P.r.m.mm.c1.p(:,Si0) + coreg.corr.P.r.m.mm.c1.p(:,Di0))/2;
            end
        end
        coreg.corr.C.r.m.mm.fp = [Ch_mS; ones(1,size(Ch_mS,2))];
        coreg.corr.C.r.m.mm.c1.p = [Ch_mc1; ones(1,size(Ch_mc1,2))];
        coreg.corr.C.r.m.vx.fp = m2v.m\coreg.corr.C.r.m.mm.fp;
        coreg.corr.C.r.m.vx.c1.p = m2v.m\coreg.corr.C.r.m.mm.c1.p;
        opt.flag = 1;
    end
end

%Add channels or optodes
if render_template %rendered on a spm template
    fd_wmm = Q*[xSPM.XYZmm; [1 1 1]]; %Fiducial coordinates in mm
    xSPM.XYZmm = fd_wmm(1:3,:); % Fiducial coordinates in mm (MNI, normalised)
    switch sel
        case 0
            % Add channels (projected onto scalp)
            % Add channels (projected onto cortex)
            for i0 = 1: Nc
                xSPM.Z = [xSPM.Z -10 10];
                xSPM.XYZmm = [xSPM.XYZmm Ch_wS(:,i0) Ch_wc1(:,i0)];
                xSPM.label = [xSPM.label ['C' int2str(i0) '_Slp'] ['C' int2str(i0) '_Ctx']];
            end
            opt.type = 'Channel';
        case 1
            for i0 = 1 : Ns
                % Add sources (projected onto scalp)
                % Add sources (projected onto cortex)
                OP_wS = coreg.corr.P.w.m.mm.p(1:3,:);
                OP_wc1 = coreg.corr.P.w.m.mm.c1.p(1:3,:);
                if (OP_wS(1,i0) ~= 0) || (OP_wS(2,i0) ~= 0) || (OP_wS(3,i0) ~= 0)
                    xSPM.Z = [xSPM.Z -10 10];
                    xSPM.XYZmm = [xSPM.XYZmm OP_wS(:,i0) OP_wc1(:,i0)];
                    xSPM.label = [xSPM.label ['S' int2str(i0) '_Slp'] ['S' int2str(i0) '_Ctx']];
                end
            end
            opt.type = 'Source';
            clear i0
        case 2
            for i0 = 1 : Nd
                % Add detectors (projected onto scalp)
                % Add detectors (projected onto cortex)
                OP_wS = coreg.corr.P.w.m.mm.p(1:3,:);
                OP_wc1 = coreg.corr.P.w.m.mm.c1.p(1:3,:);
                if (OP_wS(1,i0+Ns) ~= 0) || (OP_wS(2,i0+Ns) ~= 0) || (OP_wS(3,i0+Ns) ~= 0)
                    xSPM.Z = [xSPM.Z -10 10];
                    xSPM.XYZmm = [xSPM.XYZmm OP_wS(:,i0+Ns) OP_wc1(:,i0+Ns)];
                    xSPM.label = [xSPM.label ['D' int2str(i0) '_Slp'] ['D' int2str(i0) '_Ctx']];
                end
            end
            opt.type = 'Detector';
            clear i0
        case 3
            OP_wS = coreg.corr.P.w.m.mm.p(1:3,:);
            OP_wc1 = coreg.corr.P.w.m.mm.c1.p(1:3,:);
            if coreg_layer
                OPtmp = OP_wc1;
                Chtmp = Ch_wc1;
            else
                OPtmp = OP_wS;
                Chtmp = Ch_wS;
            end
            for i0 = 1 : Ns
                % Add sources (projected onto scalp)
                % Add sources (projected onto cortex)
                if (OPtmp(1,i0) ~= 0) || (OPtmp(2,i0) ~= 0) || (OPtmp(3,i0) ~= 0)
                    xSPM.Z = [xSPM.Z 0];
                    xSPM.XYZmm = [xSPM.XYZmm OPtmp(:,i0)];
                    xSPM.label = [xSPM.label ['S' int2str(i0)]];
                end
            end
            for i0 = 1 : Nd
                % Add detectors (projected onto scalp)
                % Add detectors (projected onto cortex)
                if (OPtmp(1,i0+Ns) ~= 0) || (OPtmp(2,i0+Ns) ~= 0) || (OPtmp(3,i0+Ns) ~= 0)
                    xSPM.Z = [xSPM.Z 10];
                    xSPM.XYZmm = [xSPM.XYZmm OPtmp(:,i0+Ns)];
                    xSPM.label = [xSPM.label ['D' int2str(i0)]];
                end
            end
            for i0 = 1: Nc
                xSPM.Z = [xSPM.Z -15];
                xSPM.XYZmm = [xSPM.XYZmm Chtmp(:,i0)];
                xSPM.label = [xSPM.label ['C' int2str(i0)]];
            end            
            opt.type = 'All Channels/Sources/Detectors';
        otherwise
            disp('Invalide choice. Must choose one of the following Channels/Sources/Detectors/All of the above');
            disp('Aborting...');
            out.NIRSmat = NIRSmat;
            return
    end    
else
    switch sel
        case 0
            % Add channels (projected onto scalp)
            % Add channels (projected onto cortex)
            for i0 = 1: Nc
                xSPM.Z = [xSPM.Z -10 10];
                xSPM.XYZmm = [xSPM.XYZmm Ch_mS(:,i0) Ch_mc1(:,i0)];
                xSPM.label = [xSPM.label ['C' int2str(i0) '_Slp'] ['C' int2str(i0) '_Ctx']];
            end
            opt.type = 'Channel';
        case 1
            OP_mS = coreg.corr.P.r.m.mm.p(1:3,:);
            OP_mc1 = coreg.corr.P.r.m.mm.c1.p(1:3,:);
            for i0 = 1 : Ns
                % Add sources (projected onto scalp)
                % Add sources (projected onto cortex)
                if (OP_mS(1,i0) ~= 0) || (OP_mS(2,i0) ~= 0) || (OP_mS(3,i0) ~= 0)
                    xSPM.Z = [xSPM.Z -10 10];
                    xSPM.XYZmm = [xSPM.XYZmm OP_mS(:,i0) OP_mc1(:,i0)];
                    xSPM.label = [xSPM.label ['S' int2str(i0) '_Slp'] ['S' int2str(i0) '_Ctx']];
                end
            end
            opt.type = 'Source';
            clear i0
        case 2
            OP_mS = coreg.corr.P.r.m.mm.p(1:3,:);
            OP_mc1 = coreg.corr.P.r.m.mm.c1.p(1:3,:);
            for i0 = 1 : Nd
                % Add detectors (projected onto scalp)
                % Add detectors (projected onto cortex)
                if (OP_mS(1,i0+Ns) ~= 0) || (OP_mS(2,i0+Ns) ~= 0) || (OP_mS(3,i0+Ns) ~= 0)
                    xSPM.Z = [xSPM.Z -10 10];
                    xSPM.XYZmm = [xSPM.XYZmm OP_mS(:,i0+Ns) OP_mc1(:,i0+Ns)];
                    xSPM.label = [xSPM.label ['D' int2str(i0) '_Slp'] ['D' int2str(i0) '_Ctx']];
                end
            end
            opt.type = 'Detector';
            clear i0
        case 3
            OP_mS = coreg.corr.P.r.m.mm.p(1:3,:);
            OP_mc1 = coreg.corr.P.r.m.mm.c1.p(1:3,:);
            if coreg_layer
                OPtmp = OP_mc1;
                Chtmp = Ch_mc1;
            else
                OPtmp = OP_mS;
                Chtmp = Ch_mS;
            end
            for i0 = 1 : Ns
                % Add sources (projected onto scalp)
                % Add sources (projected onto cortex)
                if (OPtmp(1,i0) ~= 0) || (OPtmp(2,i0) ~= 0) || (OPtmp(3,i0) ~= 0)
                    xSPM.Z = [xSPM.Z 0];
                    xSPM.XYZmm = [xSPM.XYZmm OPtmp(:,i0)];
                    xSPM.label = [xSPM.label ['S' int2str(i0)]];
                end
            end
            for i0 = 1 : Nd
                % Add detectors (projected onto scalp)
                % Add detectors (projected onto cortex)
                if (OPtmp(1,i0+Ns) ~= 0) || (OPtmp(2,i0+Ns) ~= 0) || (OPtmp(3,i0+Ns) ~= 0)
                    xSPM.Z = [xSPM.Z 10];
                    xSPM.XYZmm = [xSPM.XYZmm OPtmp(:,i0+Ns)];
                    xSPM.label = [xSPM.label ['D' int2str(i0)]];
                end
            end
            for i0 = 1: Nc
                xSPM.Z = [xSPM.Z -15];
                xSPM.XYZmm = [xSPM.XYZmm Chtmp(:,i0)];
                xSPM.label = [xSPM.label ['C' int2str(i0)]];
            end            
            opt.type = 'All Channels/Sources/Detectors';
        otherwise
            disp('Invalide choice. Must choose one of the following Channels/Sources/Detectors/All of the above');
            disp('Aborting...');
            out.NIRSmat = NIRSmat;
            return
    end
end