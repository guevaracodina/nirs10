function out = nirs_run_liom_orth_coreg(job)
%This function intends to co-register the optode position as well as 
%the activations onto orthogonal T1 image
%Ke Peng,
%2013-07-25, version 0.1, Function created

%select xSPM template
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
                    case 1
                        for i0 = 1 : Ns
                            if (NIRS.Cf.H.P.r.m.mm.p(1,i0) ~= 0) || (NIRS.Cf.H.P.r.m.mm.p(2,i0) ~= 0) || (NIRS.Cf.H.P.r.m.mm.p(3,i0) ~= 0)
                                xSPM.Z = [xSPM.Z 5];
                                xSPM.XYZmm = [xSPM.XYZmm NIRS.Cf.H.P.r.m.mm.p(:,i0)];
                                xSPM.label = [xSPM.label ['S' int2str(i0)]];
                            end
                        end
                        clear i0
                    case 2
                        for i0 = 1 : Nd
                            if (NIRS.Cf.H.P.r.m.mm.p(1,i0+Ns) ~= 0) || (NIRS.Cf.H.P.r.m.mm.p(2,i0+Ns) ~= 0) || (NIRS.Cf.H.P.r.m.mm.p(3,i0+Ns) ~= 0)
                                xSPM.Z = [xSPM.Z 5];
                                xSPM.XYZmm = [xSPM.XYZmm NIRS.Cf.H.P.r.m.mm.p(:,i0+Ns)];
                                xSPM.label = [xSPM.label ['D' int2str(i0)]];
                            end
                        end
                        clear i0
                    otherwise
                end
                
                %Transfer mm -> voxel
                tSPM.tXYZmm = [xSPM.XYZmm;ones(1,size(xSPM.XYZmm,2))];
                tSPM.tXYZ = V_render.mat\tSPM.tXYZmm;
                
                xSPM.XYZ = tSPM.tXYZ;
                tSPM.Z = xSPM.Z;
                tSPM.label = xSPM.label;
                %xSPM.XYZ = tXYZ(1:3,:);
                
                %Interpolation, as the points are too small
                
                for x0 = -2 : 2
                    for y0 = -2 : 2
                        for z0 = -2 : 2
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
                clear new_points
                
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
            spm_orthviews('Reset');
            
%             global st prevsect
%             st.Space = spm_matrix([0 0 0  0 0 -pi/2]) * st.Space;
%             prevsect = file_render;

            h = spm_orthviews('Image', file_render, [0.05 0.05 0.9 0.7]);
            spm_orthviews('AddContext', h); 
            spm_orthviews('MaxBB');
            %if ~isempty(hReg), spm_orthviews('Register', hReg); end
            spm_orthviews('AddBlobs', h, xSPM.XYZ, xSPM.Z, xSPM.M);
            spm_orthviews('Redraw');
            
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