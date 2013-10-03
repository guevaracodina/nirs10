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
                lf = job.focus_location(Idx,:)';
            catch
                disp('Problem with focus radius or coordinates for this subject');
                disp('Check that you entered the data correctly');
            end
            if rf > 30
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
%             bout = rf+5;bout2 = bout*bout;
%             bin = rf-5;bin2 = bin * bin;
%             for x0 = -bout : bout
%                 for y0 = -bout : bout
%                     for z0 = -bout : bout
%                         xyz2 = x0*x0+y0*y0+z0*z0;
%                         if  xyz2 <= bin2
%                             continue;
%                         else
%                             if xyz2 <= bout2
%                                 xSPM.Z = [xSPM.Z 1];
%                                 xSPM.XYZmm = [xSPM.XYZmm [lf(1)+x0; lf(2)+y0; lf(3)+z0]];
%                             end
%                         end
%                     end
%                 end
%             end
%             
%             x0 = -bout : bout;y0 = -bout:bout;z0 = -bout:bout;
%             lt = 2*bout+1;
%             CXYZmm = zeros(3,lt^3);
%             for t0 = 1 : lt
%                 CXYZmm(1,((t0-1)*lt^2+1):(t0*lt^2)) = x0(t0);
%                 for t00 = 1 : lt
%                     CXYZmm(2,((t00-1)*lt+1 + (t0-1)*lt^2):(t00*lt +(t0-1)*lt^2)) = y0(t00);
%                     CXYZmm(3,((t00-1)*lt+1 + (t0-1)*lt^2):(t00*lt + (t0-1)*lt^2)) = z0;
%                 end
%             end
%             disCXYZ = sqrt(CXYZmm(1,:).^2 + CXYZmm(2,:).^2 + CXYZmm(3,:).^2);
%             disidx = find((disCXYZ < bout)&(disCXYZ > bin));
%             tXYZmm = CXYZmm(:,disidx);
%             xSPM.XYZmm(1,:) = lf(1) + tXYZmm(1,:);
%             xSPM.XYZmm(2,:) = lf(2) + tXYZmm(2,:);
%             xSPM.XYZmm(3,:) = lf(3) + tXYZmm(3,:);
%             xSPM.Z = ones(1,length(disidx));
% %             for tidx = 1:length(disidx)
% %                 pidx = disidx(tidx);
% %                 xSPM.Z = [xSPM.Z 1];
% %                 xSPM.XYZmm = [xSPM.XYZmm [lf(1)+CXYZmm(1,pidx); lf(2)+CXYZmm(2,pidx); lf(3)+CXYZmm(3,pidx)]];
% %             end
            
            
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
                disp(['render_c1Extracted not found at ' dirTD]);
                disp('Render on template instead. Warning! The focus position may not be correct!');
            end

            render_template = rendered_MNI{1}.render_template;
            

            tXYZmm = [xSPM.XYZmm;ones(1,size(xSPM.XYZmm,2))];
            anatT1 = NIRS.Dt.ana.T1;
            [dirT1, fil, ext] = fileparts(anatT1);
            
            if render_template
                %Instead we want projection on the template -- need to redo
                %the spm_surf on the template.
                wT1 = NIRS.Dt.ana.wT1;
                Q = (wT1.VG.mat/wT1.Affine)/wT1.VF.mat;
                tXYZwmm = Q * tXYZmm;

                T1template = fullfile(fileparts(which('spm')),'templates','T1.nii');
                V_tp = spm_vol(T1template);
                tXYZwvx = V_tp.mat\tXYZwmm;
                xSPM.XYZ = tXYZwvx(1:3,:);
                xSPM.M = V_tp.mat;
                xSPM.DIM = V_tp.dim;
            else
                
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
                Topo_file = NIRS.TOPO;                
                try
                    load(Topo_file);
                catch
                    disp(['Could not load TOPO group file for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
                end
%                 if length(sess) == 1
%                     if sess == 0 %Group view
%                         Topo_file = NIRS.TOPO;
%                         [dir_TP,fil_TP,ext_TP] = fileparts(Topo_file);
%                         Topo_group_file = fullfile([dir_TP '\Group'],[fil_TP ext_TP]);
%                         try
%                             load(Topo_group_file);
%                         catch
%                             disp(['Could not load TOPO group file for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
%                         end
%                     else
%                         try
%                             Topo_file = NIRS.TOPO;
%                             load(Topo_file);
%                         catch
%                             disp(['Could not load TOPO file for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
%                         end
%                     end
%                 else
%                     try
%                         Topo_file = NIRS.TOPO;
%                         load(Topo_file);
%                     catch
%                         disp(['Could not load TOPO file for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
%                     end
%                 end
            end

            %do projection
            if isfield(job.proj_contrasts,'contrasts_enabled')
                disp_option.thres_sel = job.thres_sel;
                disp_option.sessions = job.proj_contrasts.contrasts_enabled.contrasts_session;
                disp_option.views = job.proj_contrasts.contrasts_enabled.contrasts_views;
                disp_option.chromophore = job.proj_contrasts.contrasts_enabled.chromophore_select;
                disp_option.activation = (job.proj_contrasts.contrasts_enabled.activation_select - 1) * 2 + 1;
                [NewNIRSDir NewNIRSFile NewNIRSExt] = fileparts(newNIRSlocation);
                disp_option.save_dir = NewNIRSDir;
                nirs_focus_contrast_render(NIRS,xSPM,brt,TOPO,rendered_MNI,disp_option);
            else
                nirs_focus_render(xSPM,brt,rendered_MNI);
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