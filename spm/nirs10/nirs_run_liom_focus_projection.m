function out = nirs_run_liom_focus_projection(job)
%This function intends to project the contrasts onto 2-D planes, try to
%locate the clinic epileptic focus.
%Ke Peng,
%2012-08-08, version 0.1, Function created
%**************************************************************************
%Loop over all subjects
%select template
ftemplate = fullfile(fileparts(which('spm')),'toolbox','nirs10','xSPMtemplate.mat');

%-Get brightness & colours
if job.render_style
    brt = NaN;
else
    brt = 1;
end
col = [];
if isfinite(brt)
    switch job.render_blobs
        case 0
            brt = 1;
        case 1
            brt = 0.75;
        case 2
            brt = 0.5;
        case 3
            brt = 0.25;
    end
    col = eye(3);
        % ask for custom colours & get rgb values
        %------------------------------------------------------------------
%         if spm_input('Which colours?','!+1','b',{'RGB','Custom'},[0 1],1)
%             for k = 1:num
%                 col(k,:) = uisetcolor(col(k,:),sprintf('Colour of blob set %d',k));
%             end
%         end
end
 
for Idx=1:size(job.NIRSmat,1)

    %Load NIRS.mat information
    try       
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'Proj_OK') || job.force_redo)
            load(ftemplate);
            xSPM.Z = [];
            xSPM.XYZ = [];
            xSPM.XYZmm = [];
            try
                if length(job.focus_radius) == 1
                    rf = job.focus_radius;
                else
                    rf = job.focus_radius(Idx);
                end
                lf = job.focus_location(Idx,:);
            catch
                disp('Problem with focus radius or coordinates for this subject');
                disp('Check that you entered the data correctly');
            end
            for x0 = -rf:rf
                for y0 = -rf:rf
                    for z0 = - rf:rf
                        if x0*x0+y0*y0+z0*z0 <= rf*rf
                            xSPM.Z = [xSPM.Z 1];
                            xSPM.XYZmm = [xSPM.XYZmm [lf(1)+x0; lf(2)+y0; lf(3)+z0]];
                        end
                    end
                end
            end
            T1 = NIRS.Dt.ana.T1;
            V1 = spm_vol(T1);
            tXYZmm = [xSPM.XYZmm;ones(1,size(xSPM.XYZmm,2))];
            tXYZ = V1.mat\tXYZmm;
            %**************************************************************
            %Normalise mm coordinates and voxel coordinates
            %Ke Peng, 2012-09-26
            %**************************************************************
            tXYZwmm = NIRS.Cf.H.P.Q * tXYZmm;
            %tXYZwvx = NIRS.Cf.H.P.Af * tXYZ;
            fwT1 = NIRS.Dt.ana.fwT1;
            V_wT1 = spm_vol(fwT1);
            tXYZwvx = V_wT1.mat\tXYZwmm;
            xSPM.XYZ = tXYZwvx(1:3,:);
            xSPM.M = V_wT1.mat;
            xSPM.DIM = V_wT1.dim;
            %**************************************************************
            %xSPM.XYZ = tXYZ(1:3,:);
            %xSPM.M = V1.mat;
            %xSPM.DIM = V1.dim;
            if ~isfield(NIRS.Dt.ana,'c1')
                [dir0 fil0] = fileparts(NIRS.Dt.ana.T1);
                [files,dirs] = spm_select('FPList',dir0,'c1*');
                found = [];
                for f0=1:size(files,1)
                    t1 = files(f0,:);
                    [dir1 fil1 ext1] = fileparts(t1);
                    if strcmp(deblank(ext1),'.nii')
                        %should be a good c1 file
                        found = deblank(t1);
                        break;
                    end
                end
                if ~isempty(found)
                    NIRS.Dt.ana.c1 = t1;
                else
                    %use T1 to compute segmentation, then calculate render file
                    NIRS.Dt.ana.c1 = t1;
                end
            end
            %check if render file on gray matter is available
            if ~isfield(NIRS.Dt.ana,'render_c1')
                [dir0 fil0] = fileparts(NIRS.Dt.ana.T1);
                [files,dirs] = spm_select('FPList',dir0,'render_c1*');%First try to find render_c1
                
                if isempty(files)
                    [files,dirs] = spm_select('FPList',dir0,'render_wc1*');%then try to find render_wc1
                end
                
                found = [];
                for f0=1:size(files,1)
                    t1 = files(f0,:);
                    [dir1 fil1 ext1] = fileparts(t1);
                    if strcmp(deblank(ext1),'.mat')
                        %should be a good c1 file
                        found = t1;
                        break;
                    end
                end
                if ~isempty(found)
                    NIRS.Dt.ana.render_c1 = deblank(t1);
                else
                    %use c1 to compute render file
                    NIRS.Dt.ana.render_c1 = t1;
                end
            end
            
            if isfield(job.proj_contrasts,'contrasts_enabled')
                %Load TOPO file

                try
                    Topo_file = NIRS.TOPO;
                    load(Topo_file);
                catch
                    disp(['Could not load TOPO file for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
                end
                
                 %Load TopoData file
                 file_topodata = fullfile(newNIRSlocation(1:(strfind(newNIRSlocation,'dataSPM')+7)), 'TopoData.mat');
                 try
                    load(file_topodata);
                 catch
                     disp(['TopoData.mat not found at ' file_topodata]);
                 end

            end
            
           
            
            %do projection
            if isfield(job.proj_contrasts,'contrasts_enabled')
                nirs_focus_contrast_render(xSPM,brt,job.thres_sel,NIRS.Dt.ana.render_c1,TOPO,rendered_MNI,job.proj_contrasts.contrasts_enabled.contrasts_session,job.proj_contrasts.contrasts_enabled.contrasts_views,NIRS.Cf.H.P.Q,NIRS.Cf.H.P.Af);
            else
                nirs_focus_render(xSPM,brt,NIRS.Dt.ana.render_c1);
            end
        end
        NIRS.flags.Proj_OK = 1;
        save(newNIRSlocation,'NIRS');
        
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not project on focus for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
    end
    
    
    
end
out.NIRSmat = job.NIRSmat;