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
            if rf > 15
                for x0 = -rf:2:rf
                    for y0 = -rf:2:rf
                        for z0 = - rf:2:rf
                            if x0*x0+y0*y0+z0*z0 <= rf*rf
                                xSPM.Z = [xSPM.Z 1];
                                xSPM.XYZmm = [xSPM.XYZmm [lf(1)+x0; lf(2)+y0; lf(3)+z0]];
                            end
                        end
                    end
                end
            else
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
            end
            
            %From the TopoData file, we can infer whether coregistration
            %was done on template or on subject
            file_topodata =  NIRS.Dt.ana.rend;
            try
                load(file_topodata);
            catch
                disp(['TopoData.mat not found at ' file_topodata]);
            end
            
            try
               [dirTD,filTD,extTD] = fileparts(file_topodata);
               rend_extracted = load(fullfile(dirTD, ['render_c1Extracted' '.mat']));
               for kk = 1 : length(rend_extracted.rend)
                    rend_extracted.rend{kk}.view_mask_2d = flipud(rendered_MNI{kk}.view_mask_2d);
                    rend_extracted.rend{kk}.render_template = rendered_MNI{kk}.render_template;
               end
               rendered_MNI = rend_extracted.rend;
               clear rend_extracted
            catch
                disp(['render_c1Extracted not found at ' dir_TD]);
            end

            render_template = rendered_MNI{1}.render_template;
            
            tXYZmm = [xSPM.XYZmm;ones(1,size(xSPM.XYZmm,2))];
            anatT1 = NIRS.Dt.ana.T1;
            [dirT1, fil, ext] = fileparts(anatT1);
            
            if render_template
                %Instead we want projection on the template -- need to redo
                %the spm_surf on the template.
                
                %**************************************************************
                %Normalise mm coordinates and voxel coordinates
                %Ke Peng, 2012-09-26
                %**************************************************************
                wT1 = NIRS.Dt.ana.wT1;
                Q = (wT1.VG.mat/wT1.Affine)/wT1.VF.mat;
                tXYZwmm = Q * tXYZmm;
                fwT1 = fullfile(dirT1,['w' fil ext(1:4)]);

                
                T1template = fullfile(fileparts(which('spm')),'templates','T1.nii');
                V_tp = spm_vol(T1template);
                tXYZwvx = V_tp.mat\tXYZwmm;
                xSPM.XYZ = tXYZwvx(1:3,:);
                xSPM.M = V_tp.mat;
                xSPM.DIM = V_tp.dim;
%                 V_wT1 = spm_vol(fwT1);
%                 tXYZwvx = V_wT1.mat\tXYZwmm;
%                 xSPM.XYZ = tXYZwvx(1:3,:);
%                 xSPM.M = V_wT1.mat;
%                 xSPM.DIM = V_wT1.dim;
%                 fwc1 =  fullfile(dirT1,['wc1' fil ext(1:4)]);               
                %**************************************************************
%                 %check if render file on gray matter is available
%                 if ~isfield(NIRS.Dt.ana,'render_wc1')
%                     %generate rendered image
%                     spm_surf(fwc1,1);
%                     NIRS.Dt.ana.render_wc1 = fullfile(dirT1,['render_wc1' fil '.mat']);                    
%                 end
%                 rend_file = NIRS.Dt.ana.render_wc1;
            else
                %                 fc1 =  fullfile(dirT1,['c1' fil ext(1:4)]);
%                 %check if segmentation is required
%                 if ~spm_existfile(fc1)
%                     disp('Need to run segmentation');
%                 end
%                 
%                 %check if render file on gray matter is available
%                 if ~isfield(NIRS.Dt.ana,'render_c1')
%                     %generate rendered image
%                     spm_surf(fc1,1);
%                     NIRS.Dt.ana.render_c1 = fullfile(dirT1,['render_c1' fil '.mat']);                  
%                 end
%                 rend_file = NIRS.Dt.ana.render_c1;
                file_c1 = fullfile(dirT1,['c1' fil ext(1:4)]);
                V_c1 = spm_vol(file_c1);
                
                tXYZmm = [xSPM.XYZmm;ones(1,size(xSPM.XYZmm,2))];
                tXYZ = V_c1.mat\tXYZmm;
                
                xSPM.XYZ = tXYZ(1:3,:);
                xSPM.M = V_c1.mat;
                xSPM.DIM = V_c1.dim;
            end
            
            
            
            if isfield(job.proj_contrasts,'contrasts_enabled')
                %Load TOPO file

                try
                    Topo_file = NIRS.TOPO;
                    load(Topo_file);
                catch
                    disp(['Could not load TOPO file for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
                end
            end
            
            
            
           
            
            %do projection
            if isfield(job.proj_contrasts,'contrasts_enabled')
                nirs_focus_contrast_render(xSPM,brt,job.thres_sel,TOPO,rendered_MNI,job.proj_contrasts.contrasts_enabled.contrasts_session,job.proj_contrasts.contrasts_enabled.contrasts_views);
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