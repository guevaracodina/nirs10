function out = nirs_run_liom_orth_coreg(job)
%This function intends to co-register the optode position as well as 
%the activations onto orthogonal T1 image
%Ke Peng,
%2013-07-25, version 0.1, Function created

%select xSPM template

global opt_info

ftemplate = fullfile(fileparts(which('spm')),'toolbox','nirs10','xSPMtemplate.mat');

for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try       
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'OrthCoreg_OK') || job.force_redo)
            load(ftemplate);
            xSPM.Z = [];
            xSPM.XYZ = [];
            xSPM.XYZmm = [];
            xSPM.label = {};
            file_render = job.render_image{1};
            if isempty(file_render)
                file_render = NIRS.Dt.ana.T1;
                disp('Render image not specified. T1 image is chosen by default');
            end
            
            if isfield(job.coreg_type, 'NIRS_channels_optodes')
                
                V_render = spm_vol(file_render);
                                
                %Add fiducial points
                xSPM.Z = [15 15 15]; %Fiducial value
                xSPM.XYZmm = NIRS.Cf.H.F.r.m.mm.p;
                xSPM.label = {'F-N','F-L','F-R'};
                
                Ns = NIRS.Cf.H.S.N;
                Nd = NIRS.Cf.H.D.N;
                %Add channels or optodes  
                
                switch job.coreg_type.NIRS_channels_optodes.channel_optode_select
                    case 0
                        % Add channels (projected onto scalp)
                        % Add channels (projected onto cortex)
                        Nc = size(NIRS.Cf.H.C.id,2)/2;
                        for i0 = 1 : Nc
                            Si0 = NIRS.Cf.H.C.id(2,i0);
                            Di0 = NIRS.Cf.H.C.id(3,i0)+Ns;
                            Ci0_s = (NIRS.Cf.H.P.r.m.mm.p(:,Si0) + NIRS.Cf.H.P.r.m.mm.p(:,Di0))/2;
                            Ci0_c = (NIRS.Cf.H.P.r.m.mm.c1.p(:,Si0) + NIRS.Cf.H.P.r.m.mm.c1.p(:,Di0))/2;
                            xSPM.Z = [xSPM.Z -10 10];
                            xSPM.XYZmm = [xSPM.XYZmm Ci0_s Ci0_c];
                            xSPM.label = [xSPM.label ['C' int2str(i0)] ['C' int2str(i0)]];
                        end
                        opt_info.type = 'Channel';
                        opt_info.label = xSPM.label;
                        clear i0 Si0 Di0 Ci0
                    case 1
                        for i0 = 1 : Ns
                            % Add sources (projected onto scalp)
                            % Add sources (projected onto cortex)
                            if (NIRS.Cf.H.P.r.m.mm.p(1,i0) ~= 0) || (NIRS.Cf.H.P.r.m.mm.p(2,i0) ~= 0) || (NIRS.Cf.H.P.r.m.mm.p(3,i0) ~= 0)
                                xSPM.Z = [xSPM.Z -10 10];
                                xSPM.XYZmm = [xSPM.XYZmm NIRS.Cf.H.P.r.m.mm.p(:,i0) NIRS.Cf.H.P.r.m.mm.c1.p(:,i0)];
                                xSPM.label = [xSPM.label ['S' int2str(i0)] ['S' int2str(i0)]];
                            end                         
                        end
                        opt_info.type = 'Source';
                        opt_info.label = xSPM.label;
                        clear i0
                    case 2
                        for i0 = 1 : Nd
                            % Add detectors (projected onto scalp)
                            % Add detectors (projected onto cortex)
                            if (NIRS.Cf.H.P.r.m.mm.p(1,i0+Ns) ~= 0) || (NIRS.Cf.H.P.r.m.mm.p(2,i0+Ns) ~= 0) || (NIRS.Cf.H.P.r.m.mm.p(3,i0+Ns) ~= 0)
                                xSPM.Z = [xSPM.Z -10 10];
                                xSPM.XYZmm = [xSPM.XYZmm NIRS.Cf.H.P.r.m.mm.p(:,i0+Ns) NIRS.Cf.H.P.r.m.mm.c1.p(:,i0+Ns)];
                                xSPM.label = [xSPM.label ['D' int2str(i0)] ['D' int2str(i0)]];
                            end
                        end
                        opt_info.type = 'Detector';
                        opt_info.label = xSPM.label;
                        clear i0
                    case 3
                        Nc = size(NIRS.Cf.H.C.id,2)/2;
                        for i0 = 1 : Ns
                            if (NIRS.Cf.H.P.r.m.mm.p(1,i0) ~= 0) || (NIRS.Cf.H.P.r.m.mm.p(2,i0) ~= 0) || (NIRS.Cf.H.P.r.m.mm.p(3,i0) ~= 0)
                                xSPM.Z = [xSPM.Z 0];
                                switch job.coreg_layer
                                    case 0
                                        xSPM.XYZmm = [xSPM.XYZmm NIRS.Cf.H.P.r.m.mm.p(:,i0)];
                                    case 1
                                        xSPM.XYZmm = [xSPM.XYZmm NIRS.Cf.H.P.r.m.mm.c1.p(:,i0)];
                                    otherwise
                                        xSPM.XYZmm = [xSPM.XYZmm NIRS.Cf.H.P.r.m.mm.p(:,i0)];
                                end
                                xSPM.label = [xSPM.label ['S' int2str(i0)]];
                            end
                        end
                        for i0 = 1 : Nd
                            if (NIRS.Cf.H.P.r.m.mm.p(1,i0+Ns) ~= 0) || (NIRS.Cf.H.P.r.m.mm.p(2,i0+Ns) ~= 0) || (NIRS.Cf.H.P.r.m.mm.p(3,i0+Ns) ~= 0)
                                xSPM.Z = [xSPM.Z 10];
                                switch job.coreg_layer
                                    case 0 
                                        xSPM.XYZmm = [xSPM.XYZmm NIRS.Cf.H.P.r.m.mm.p(:,i0+Ns)];
                                    case 1
                                        xSPM.XYZmm = [xSPM.XYZmm NIRS.Cf.H.P.r.m.mm.c1.p(:,i0+Ns)];
                                    otherwise
                                        xSPM.XYZmm = [xSPM.XYZmm NIRS.Cf.H.P.r.m.mm.p(:,i0+Ns)];
                                end
                                xSPM.label = [xSPM.label ['D' int2str(i0)]];
                            end
                        end
                        
                        L.Z = [];
                        L.XYZmm = [];
                        L.label = [];
                        
                        for i0 = 1 : Nc
                            Si0 = NIRS.Cf.H.C.id(2,i0);
                            Di0 = NIRS.Cf.H.C.id(3,i0)+Ns;
                            switch job.coreg_layer
                                case 0
                                    Ci0 = (NIRS.Cf.H.P.r.m.mm.p(:,Si0) + NIRS.Cf.H.P.r.m.mm.p(:,Di0))/2;
                                case 1
                                    Ci0 = (NIRS.Cf.H.P.r.m.mm.c1.p(:,Si0) + NIRS.Cf.H.P.r.m.mm.c1.p(:,Di0))/2;
                                otherwise
                                    Ci0 = (NIRS.Cf.H.P.r.m.mm.p(:,Si0) + NIRS.Cf.H.P.r.m.mm.p(:,Di0))/2;
                            end
                            xSPM.Z = [xSPM.Z -15];
                            %xSPM.Z(1,i0+3) = -15;
                            xSPM.XYZmm = [xSPM.XYZmm Ci0];
                            %xSPM.XYZmm(:,i0+3) = Ci0;
                            xSPM.label = [xSPM.label ['C' int2str(i0) ': S' int2str(Si0) '-D' int2str(Di0-Ns)]];
                                              
                        end
                        
                        opt_info.type = 'All Channels/Sources/Detectors';
                        opt_info.label = xSPM.label;
                       
                    otherwise
                        disp('Invalide choice. Must choose one of the following Channels/Sources/Detectors/All of the above');
                        disp('Aborting...');
                        out.NIRSmat = NIRSmat;
                        return
                end
                
                %Transfer mm -> voxel
                tSPM.tXYZmm = [xSPM.XYZmm;ones(1,size(xSPM.XYZmm,2))];
                tSPM.tXYZ = V_render.mat\tSPM.tXYZmm;
                opt_info.XYZmm = tSPM.tXYZmm;
                
                xSPM.XYZ = tSPM.tXYZ;
                tSPM.Z = xSPM.Z;
                tSPM.label = xSPM.label;
                %xSPM.XYZ = tXYZ(1:3,:);
                
                %Interpolation, as the points are too small
                
                for x0 = -3 : 3
                    for y0 = -3 : 3
                        for z0 = -3 : 3
                            if x0 == 0 && y0 == 0 && z0 == 0
                                continue;
                            else
                                new_points = [tSPM.tXYZ(1,:)+x0; tSPM.tXYZ(2,:)+y0; tSPM.tXYZ(3,:)+z0; tSPM.tXYZ(4,:)];
                                xSPM.XYZ = [xSPM.XYZ new_points];
                                xSPM.Z = [xSPM.Z tSPM.Z];
                                xSPM.label = [xSPM.label tSPM.label];
                            end
                        end
                    end
                end
                
                clear new_points x0 y0 z0
                
                %Transfer back to mm
                xSPM.XYZmm = [];
                xSPM.XYZmm = V_render.mat * xSPM.XYZ;
                xSPM.XYZ = xSPM.XYZ(1:3,:);
                xSPM.XYZmm = xSPM.XYZmm(1:3,:);

                xSPM.M = V_render.mat;
                xSPM.DIM = V_render.dim;
            else
                    
            end
            
            %Codes from spm_sections.m
            Fgraph = spm_figure('GetWin','Graphics');
            spm_results_ui('Clear',Fgraph);
            nirs_orthviews('Reset');
            
            switch job.coreg_type.NIRS_channels_optodes.channel_optode_select
                case 0
                    color_ity = 'Red: Channels onto Scalp; Yellow: Channels onto cortex';
                    hCol = axes('Parent',Fgraph,'Position',[0.02 0.80 0.96 0.02],'Visible','off');
                    text(0.5,0,color_ity,'Parent',hCol,...
                    'HorizontalAlignment','center','VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);                
                case 1
                    color_ity = 'Red: Sources onto Scalp; Yellow: Sources onto cortex';
                    hCol = axes('Parent',Fgraph,'Position',[0.02 0.80 0.96 0.02],'Visible','off');
                    text(0.5,0,color_ity,'Parent',hCol,...
                    'HorizontalAlignment','center','VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                case 2
                    color_ity = 'Red: Sources onto Scalp; Yellow: Sources onto cortex';
                    hCol = axes('Parent',Fgraph,'Position',[0.02 0.80 0.96 0.02],'Visible','off');
                    text(0.5,0,color_ity,'Parent',hCol,...
                    'HorizontalAlignment','center','VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                case 3
                    color_ity = 'Red-Channels; Orange-Sources; Yellow-Detectors';
                    hCol = axes('Parent',Fgraph,'Position',[0.02 0.80 0.96 0.02],'Visible','off');
                    text(0.5,0,color_ity,'Parent',hCol,...
                    'HorizontalAlignment','center','VerticalAlignment','baseline',...
                    'FontWeight','Bold','FontSize',14);
                otherwise
            end

%             global st prevsect
%             st.Space = spm_matrix([0 0 0  0 0 -pi/2]) * st.Space;
%             prevsect = file_render;

            h = nirs_orthviews('Image', file_render, [0.05 0.05 0.9 0.7]);
            nirs_orthviews('AddContext', h); 
            nirs_orthviews('MaxBB');
            %if ~isempty(hReg), spm_orthviews('Register', hReg); end
            nirs_orthviews('AddBlobs', h, xSPM.XYZ, xSPM.Z, xSPM.M);
            nirs_orthviews('Redraw');
            
            %Add button
            if job.coreg_type.NIRS_channels_optodes.channel_optode_select ~= 3
                uicontrol(Fgraph,'Style','PushButton','String','Correspondence','FontSize',12,...
                'ToolTipString',...
                'find corresponding place on cortex/scalp of an optode',...
                'Callback','nirs_orthviews(''correspondence'');',...
                'Interruptible','on','Enable','on',...
                'Position',[500 850 150 050]);
            end
        end
        NIRS.flags.OrthCoreg_OK = 1;
        save(newNIRSlocation,'NIRS');
        
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not Coregister for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
    end
    
    
    
end
out.NIRSmat = job.NIRSmat;
end
