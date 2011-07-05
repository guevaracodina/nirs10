function out = nirs_run_MCsegment(job)
% Complete image segmentation after New Segment
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% 2010-06

%outVsegmented = {};
outNIRSmat = {};
try
    outNIRSmat = job.NIRSmat_optional;%%%%%%%job.NIRSmat_optional;
    nsubj = length(outNIRSmat);
    NIRSok = 1;
catch
    NIRSok = 0;
    nsubj = length(job.image_in);
end

for Idx=1:nsubj
    if NIRSok
        try
            load(outNIRSmat{Idx});
        catch
            disp(['Could not load NIRS.mat for ' int2str(Idx) 'th subject.']);
        end
    end
    if ~isempty(job.image_in{1,1})
        try
            V.fname = job.image_in{Idx,:};
        catch
            V.fname = NIRS.Dt.ana.T1;
        end
    else
        try
            V.fname = NIRS.Dt.ana.T1;
        catch
            disp(['Could not find anatomical image for ' int2str(Idx) 'th subject.']);
        end
    end
    
    
    % USER OPTIONS %
    thresh_as = job.thresh_as;
    rebel_surrounding = job.rebel_surrounding;
    rebel_thresh_hs = job.rebel_thresh_hs;

    c1method = job.wtm.sorting_method;
    c2method = job.grm.sorting_method;
    c3method = job.csf.sorting_method;
    c4method = job.skl.sorting_method;
    c5method = job.skn.sorting_method;

    output_autonaming = job.output_autonaming;
    if output_autonaming==0
        output_prefix = strcat(int2str(c1method),...
            int2str(c2method),...
            int2str(c3method),...
            int2str(c4method),...
            int2str(c5method));
    else
        output_prefix = job.output_prefix;
    end

    se_size_pi     = job.process_image.se_size_pi;
    gaussfilt_size = job.process_image.gaussfilt_size;
    gaussfilt_sdev = job.process_image.gaussfilt_sdev;
    
    % Check if already run
    [dir1, file1, dummy] = fileparts(V.fname);
    %tmpf = spm_select('List',dir1,['_segmented_' file1]);
    %if ~isempty(tmpf)
    if spm_existfile(fullfile(dir1,[output_prefix,'_segmented_',file1,'.nii']))
        %         NIRS.Cs.mcs.seg = tmpf(1,:);
        %MC Segmentation already done, skipping     
        if NIRSok
            NIRS.Dt.ana.T1seg = fullfile(dir1,[output_prefix,'_segmented_',file1,'.nii']);
            disp(['MC Segment already run for file ' NIRS.Dt.ana.T1seg ' - skipping.'])
            save(outNIRSmat{Idx},'NIRS');
        end
        
    else
        % Obtain tissue segmentation with SPM's New Segment, if not already done
        nirs_spm_NewSegment(V.fname);

        % Then process image and save segmented image
        try
            [dir,name] = fileparts(V.fname);
            V = spm_vol(V.fname);

            spm_progress_bar('Init',50,'Completing MC Segment','Progress');

            %%% %%% %%% %%% %%% verifier que les fichiers existent %%% %%% %%% %%% %%%
            nirs_MCsegment_process_image(fullfile(dir,['c1',name,'.nii']),...
                c1method, se_size_pi, gaussfilt_size, gaussfilt_sdev);
            nirs_MCsegment_process_image(fullfile(dir,['c2',name,'.nii']),...
                c2method, se_size_pi, gaussfilt_size, gaussfilt_sdev);
            nirs_MCsegment_process_image(fullfile(dir,['c3',name,'.nii']),...
                c3method, se_size_pi, gaussfilt_size, gaussfilt_sdev);

            spm_progress_bar('Set',3);

            nirs_MCsegment_process_image(fullfile(dir,['c4',name,'.nii']),...
                c4method, se_size_pi, gaussfilt_size, gaussfilt_sdev);
            nirs_MCsegment_process_image(fullfile(dir,['c5',name,'.nii']),...
                c5method, se_size_pi, gaussfilt_size, gaussfilt_sdev);

            spm_progress_bar('Set',5);

            % gathering of layers
            % NewSegment is considered to provide a perfect segmentation
            % Control (sum of processed_ci) : according to the number of layers
            % to which belongs each
            % voxel, one determinate the right layer

            % Gray matter
            V1 = spm_vol(fullfile(dir,['processed_c1',name,'.nii']));
            Y1 = spm_read_vols(V1);
            clear V1;
            spm_progress_bar('Set',6);

            % White matter
            V2 = spm_vol(fullfile(dir,['processed_c2',name,'.nii']));
            Y2 = spm_read_vols(V2);
            clear V2;
            spm_progress_bar('Set',7);

            % CSF
            V3 = spm_vol(fullfile(dir,['processed_c3',name,'.nii']));
            Y3 = spm_read_vols(V3);
            clear V3;
            spm_progress_bar('Set',8);

            % Skull
            V4 = spm_vol(fullfile(dir,['processed_c4',name,'.nii']));
            Y4 = spm_read_vols(V4);
            clear V4;
            spm_progress_bar('Set',9);

            % Skin
            V5 = spm_vol(fullfile(dir,['processed_c5',name,'.nii']));
            Y5 = spm_read_vols(V5);
            clear V5;
            spm_progress_bar('Set',10);

            % Gray matter
            V1m = spm_vol(fullfile(dir,['c1',name,'.nii']));
            Y1m = spm_read_vols(V1m);
            clear V1m;
            spm_progress_bar('Set',11);

            % White matter
            V2m = spm_vol(fullfile(dir,['c2',name,'.nii']));
            Y2m = spm_read_vols(V2m);
            clear V2m;
            spm_progress_bar('Set',12);

            % CSF
            V3m = spm_vol(fullfile(dir,['c3',name,'.nii']));
            Y3m = spm_read_vols(V3m);
            clear V3m;
            spm_progress_bar('Set',13);

            % Skull
            V4m = spm_vol(fullfile(dir,['c4',name,'.nii']));
            Y4m = spm_read_vols(V4m);
            clear V4m;
            spm_progress_bar('Set',14);

            % Skin
            V5m = spm_vol(fullfile(dir,['c5',name,'.nii']));
            Y5m = spm_read_vols(V5m);
            %clear V5m; This image is required below!
            spm_progress_bar('Set',15);

            Ysegmented = zeros(size(Y5m));

            % Give the number of layers each voxel belongs to (between 0 and 2)
            [dummy,belongs2n_V] = spm_imcalc_ui([fullfile(dir,['processed_c1',name,'.nii']);...
                fullfile(dir,['processed_c2',name,'.nii']);...
                fullfile(dir,['processed_c3',name,'.nii']);...
                fullfile(dir,['processed_c4',name,'.nii']);
                fullfile(dir,['processed_c5',name,'.nii'])],...
                fullfile(dir,'belongs2n_processed_ci.nii'),...
                'i1+i2+i3+i4+i5');
            belongs2n_Y = spm_read_vols(belongs2n_V);

            spm_progress_bar('Set',17);

            % Head shadow is determined to know if a voxel is situated in the
            % head or not
            se_size_hs = job.head_shadow.se_size_hs;
            thresh_hs  = job.head_shadow.thresh_hs;

            head_shadow_fname = nirs_MCsegment_head_shadow(V.fname, se_size_hs, thresh_hs);

            head_shadow_V = spm_vol(head_shadow_fname);
            head_shadow_Y = spm_read_vols(head_shadow_V);

            spm_progress_bar('Set',19);

            %--> belongs_2n==1 : voxel is directly attributed
            Ysegmented( belongs2n_Y>=0.5 & belongs2n_Y<1.51 & Y5>=thresh_as ) = 5;%Y5m
            Ysegmented( belongs2n_Y>=0.5 & belongs2n_Y<1.51 & Y4>=thresh_as ) = 4;%Y4m

            Ysegmented( belongs2n_Y>=0.5 & belongs2n_Y<1.51 & Y2>=thresh_as ) = 2;
            Ysegmented( belongs2n_Y>=0.5 & belongs2n_Y<1.51 & Y1>=thresh_as ) = 1;
            Ysegmented( belongs2n_Y>=0.5 & belongs2n_Y<1.51 & Y3>=thresh_as ) = 3;
            spm_progress_bar('Set',20);
            %--> belongs_2n==2 : voxel is attributed according to the biggest
            % probability (of the random variable) [Y4m and Y5m are used because
            % Y4 and Y5 are binary images]
            Ysegmented( belongs2n_Y>=1.51 & Y4m>=Y1 & Y4m>=Y2 & Y4m>=Y3 & Y4m>=Y5m ) = 4;
            Ysegmented( belongs2n_Y>=1.51 & Y5m>=Y1 & Y5m>=Y2 & Y5m>=Y3 & Y5m>=Y4m ) = 5;

            Ysegmented( belongs2n_Y>=1.51 & Y2>=Y1 & Y2>=Y3 & Y2>=Y4m & Y2>=Y5m ) = 2;
            Ysegmented( belongs2n_Y>=1.51 & Y1>=Y2 & Y1>=Y3 & Y1>=Y4m & Y1>=Y5m ) = 1;
            Ysegmented( belongs2n_Y>=1.51 & Y3>=Y1 & Y3>=Y2 & Y3>=Y4m & Y3>=Y5m ) = 3;
            spm_progress_bar('Set',25);
            %--> belongs_2n==0 : voxel has'nt been attributed to any
            %layer...
            % All those belonging to the head shadow are attributed depending
            % on the values of random variables
            Ysegmented( belongs2n_Y<0.5 & head_shadow_Y==1 & Y4m>=Y1 & Y4m>=Y2 & Y4m>=Y3 & Y4m>=Y5m & Y4m~=0 ) = 4;
            Ysegmented( belongs2n_Y<0.5 & head_shadow_Y==1 & Y5m>=Y1 & Y5m>=Y2 & Y5m>=Y3 & Y5m>=Y4m & Y5m~=0 ) = 5;
            Ysegmented( belongs2n_Y<0.5 & head_shadow_Y==1 & Y2>=Y1 & Y2>=Y3 & Y2>=Y4m & Y2>=Y5m & Y2~=0 ) = 2;
            Ysegmented( belongs2n_Y<0.5 & head_shadow_Y==1 & Y1>=Y2 & Y1>=Y3 & Y1>=Y4m & Y1>=Y5m & Y1~=0 ) = 1;
            Ysegmented( belongs2n_Y<0.5 & head_shadow_Y==1 & Y3>=Y1 & Y3>=Y2 & Y3>=Y4m & Y3>=Y5m & Y3~=0 ) = 3;
            spm_progress_bar('Set',30);
            %--> for rebel voxels
            rebel_count = 0;
            % rebel_count_processed = 0;

            for i=1+rebel_surrounding:size(Ysegmented,1)-rebel_surrounding
                % non attributed voxels are detected
                [posAlpha,posAleph]=find(squeeze(head_shadow_Y(i,:,:))>=rebel_thresh_hs & squeeze(Ysegmented(i,:,:))==0);
                posAlpha = floor(posAlpha);
                posAleph = floor(posAleph);

                % processing require a surrounding (hence to reduce the size of the image processed)
                for i_alpha=1:size(posAlpha,1)
                    if posAlpha(i_alpha,1)<=rebel_surrounding
                        posAlpha(i_alpha,1)=1+rebel_surrounding;
                    elseif posAlpha(i_alpha,1)>=V5m.dim(2)-rebel_surrounding+1
                        posAlpha(i_alpha,1)=V5m.dim(2)-rebel_surrounding;
                    end
                end
                for i_aleph=1:size(posAleph,1)
                    if posAleph(i_aleph,1)<=rebel_surrounding
                        posAleph(i_aleph,1)=1+rebel_surrounding;
                    elseif posAleph(i_aleph,1)>=V5m.dim(3)-rebel_surrounding+1
                        posAleph(i_aleph,1)=V5m.dim(3)-rebel_surrounding;
                    end
                end
                spm_progress_bar('Set',40);
                % on each surrounding, one compute the number of voxels of
                % each layer
                for j=1:size(posAlpha,1)
                    nbrhood = Ysegmented(i-rebel_surrounding:i+rebel_surrounding,posAlpha(j)-rebel_surrounding:posAlpha(j)+rebel_surrounding,posAleph(j)-rebel_surrounding:posAleph(j)+rebel_surrounding);
                    % the most widespread value is the one given to the
                    % central voxel
                    rebel_value_count = zeros(6,1);
                    for i_value =1:6
                        rebel_value_count(i_value) = size(find(nbrhood==(i_value-1)),1);
                    end
                    [dummy,layer] = max(rebel_value_count);
                    layer = layer-1;
                    % if layer ==0 hence air, then values from the
                    % NewSegment layers are used
                    if layer==0
                        v_c1 = Y1(i,posAlpha(j),posAleph(j));
                        v_c2 = Y2(i,posAlpha(j),posAleph(j));
                        v_c3 = Y3(i,posAlpha(j),posAleph(j));
                        v_c4 = Y4m(i,posAlpha(j),posAleph(j));
                        v_c5 = Y5m(i,posAlpha(j),posAleph(j));
                        [dummy,layer] = max([v_c1,v_c2,v_c3,v_c4,v_c5]);
                        Ysegmented(i,posAlpha(j),posAleph(j)) = layer;
                    else
                        Ysegmented(i,posAlpha(j),posAleph(j)) = layer;
                    end
                end
                rebel_count = rebel_count +  size(posAlpha,1);
            end
            disp(['Remaining rebel voxels before processing : ',int2str(rebel_count)]);
            rebel_count_processed = 0;
            for i=1+rebel_surrounding:size(Ysegmented,1)-rebel_surrounding
                [posAlpha,dummy]=find(squeeze(head_shadow_Y(i,:,:))>=rebel_thresh_hs & squeeze(Ysegmented(i,:,:))==0);
                rebel_count_processed = rebel_count_processed +  size(posAlpha,1);
            end

            % Last try at filling holes in c4 (cranium... usually holes around
            % the sinuses)
            cranium = (Ysegmented==4);
            cranium = imfill(cranium,'holes'); % fill holes
            Ysegmented = Ysegmented .* abs(1-cranium); % 0-out layer c4
            Ysegmented = Ysegmented + 10*cranium; % then add cranium with holes filled
            % then reset original layers where they were
            Ysegmented(Ysegmented==15) = 5;
            Ysegmented(Ysegmented==10) = 4;
            Ysegmented(Ysegmented==13) = 3;
            Ysegmented(Ysegmented==12) = 2;
            Ysegmented(Ysegmented==11) = 1;

            disp(['Remaining rebel voxels after processing : ',int2str(rebel_count_processed)]);
            
            % calcul des covariances des voxels (voir Diffuse optical tomography with a priori anatomical information, Guven, 2005)
            % on calcule l ecart type pour l'ense;ble des voxels
            % appartenant a une couche
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% A CODER %

            % save of nifti image :
            Vsegmented = struct('fname',fullfile(dir,[output_prefix,'_segmented_',name,'.nii']),...
                'dim',  V.dim,...
                'dt',   V.dt,...
                'pinfo',V.pinfo,...
                'mat',  V.mat);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%% pinfo(1,1) doit etre egal a 1...
            Vsegmented.pinfo(1,1)=1;
            spm_write_vol(Vsegmented, Ysegmented);

            spm_progress_bar('Set',50);
            spm_progress_bar('Clear');

            % Save name of segmented image in NIRS matrix
            NIRS.Dt.ana.T1seg = fullfile(dir,[output_prefix,'_segmented_',name,'.nii']);
            if NIRSok
                %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 NIRS.Cs.mcs.seg = Vsegmented.fname;
                save(outNIRSmat{Idx},'NIRS');
            else
                %what to do if we don't have a NIRS structure already?
            end

        catch
            disp(['MC Segment failed to run for subject ' int2str(Idx)]);
        end
        
    end
    %Returns the location of the processed image
    %outVsegmented = [outVsegmented; Vsegmented.fname];
end
%out.Vsegmented = outVsegmented;
out.NIRSmat = outNIRSmat;
