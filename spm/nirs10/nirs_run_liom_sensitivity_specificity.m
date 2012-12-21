function out = nirs_run_liom_sensitivity_specificity(job)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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

total_NIRS = size(job.NIRSmat,1);
subject_sen_spe = cell(total_NIRS,1);

for Idx=1:size(job.NIRSmat,1)
    
    %Load NIRS.mat information
    try
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'Proj_OK') || job.force_redo)
            
            %**************************************************************
            [dir0 subject_name ext0] = fileparts(NIRS.Dt.s.p);
            subject_sen_spe{Idx,1}.subj_name = subject_name;
            
            if job.sessions_sel == 0
                sessions = 1:length(NIRS.Dt.fir.Sess);
            end
            views = [2 3 4 5];
            %**************************************************************
            
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
                if size(job.focus_location,2) == 1
                    lf = job.focus_location(:,Idx)';
                else
                    lf = job.focus_location(Idx,:);
                end
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
%             for kk=1:6
%             rendered_MNI{kk}.render_template = 1;
%             end
            render_template = rendered_MNI{1}.render_template;
                 
            T1 = NIRS.Dt.ana.T1;
            V1 = spm_vol(T1);
            tXYZmm = [xSPM.XYZmm;ones(1,size(xSPM.XYZmm,2))];
            tXYZ = V1.mat\tXYZmm;
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
                V_wT1 = spm_vol(fwT1);
                tXYZwvx = V_wT1.mat\tXYZwmm;
                xSPM.XYZ = tXYZwvx(1:3,:);
                xSPM.M = V_wT1.mat;
                xSPM.DIM = V_wT1.dim;
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
                xSPM.XYZ = tXYZ(1:3,:);
                xSPM.M = V1.mat;
                xSPM.DIM = V1.dim;
            end           
            %Load TOPO file
            try
                Topo_file = NIRS.TOPO;
                load(Topo_file);
            catch
                disp(['Could not load TOPO file for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
            end
            %sensitivity and specificity check for subject
            [sen_spe render_proj] = nirs_sensitivity_specificity_calcu(...
                job.thres_sel,xSPM,brt,TOPO,rendered_MNI,sessions,views);
            NIRS.sen_spe.sen_spe = sen_spe;
            NIRS.sen_spe.render_img = render_proj;
            
            %For HbR only
            subject_sen_spe{Idx,1}.sensitivity = sen_spe.sen_spe_over_all_sessions{1}.sensitivity;
            subject_sen_spe{Idx,1}.specificity = sen_spe.sen_spe_over_all_sessions{1}.specificity;
        end
        NIRS.flags.Proj_OK = 1;
        save(newNIRSlocation,'NIRS');
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not do sensitivity and/or specificity analysis for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
    end
end

sen_spe_NaN = 0;
sen_1 = 0;
spe_1 = 0;

disp('====================================================================')
disp('==================Sensitivity & Specificity Report==================');
disp('--------------------------------------------------------------------');
disp('Subject name                           Sensitivity       Specificity');
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');
for Idx = 1:total_NIRS
    if isnan(subject_sen_spe{Idx,1}.sensitivity)
        sen_spe_NaN = sen_spe_NaN + 1;
        disp([subject_sen_spe{Idx,1}.subj_name '                           Unknown          Unknown']);
    else
        disp([subject_sen_spe{Idx,1}.subj_name '     ' int2str(subject_sen_spe{Idx,1}.sensitivity) '                  ' int2str(subject_sen_spe{Idx,1}.specificity)]);
        if subject_sen_spe{Idx,1}.sensitivity == 1
            sen_1 = sen_1 + 1;
        end
        if subject_sen_spe{Idx,1}.specificity == 1
            spe_1 = spe_1 + 1;
        end
    end
    
end
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');
disp(['Subtotal:                                  ' int2str(sen_1) '                  ' int2str(spe_1)]);
disp('--------------------------------------------------------------------');
disp(['Coverage includes focus:                   ' int2str(total_NIRS-sen_spe_NaN) '                  ' int2str(total_NIRS-sen_spe_NaN)]);
disp(['Total subjects:                            ' int2str(total_NIRS) '                  ' int2str(total_NIRS)]);
disp(['Valid subjects ratio:                     ' num2str((1-sen_spe_NaN/total_NIRS)*100,'%4.2f\n') '%' '             ' num2str((1-sen_spe_NaN/total_NIRS)*100,'%4.2f\n') '%']);
disp(['Good ratio (over valid subjects):        ' num2str((sen_1/(total_NIRS-sen_spe_NaN))*100,'%4.2f\n') '%' '             ' num2str((spe_1/(total_NIRS-sen_spe_NaN))*100,'%4.2f\n') '%']);
disp('--------------------------------------------------------------------');
disp('=========================End of Report==============================');
disp('====================================================================');

out.NIRSmat = job.NIRSmat;


