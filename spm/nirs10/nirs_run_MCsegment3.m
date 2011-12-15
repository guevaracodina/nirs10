function out = nirs_run_MCsegment3(job)
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
try 
    force_reprocess = job.force_reprocess;
catch
    force_reprocess = 0;
end
for Idx=1:nsubj
    if NIRSok
        try
            load(outNIRSmat{Idx});
        catch
            if ~isempty(outNIRSmat{Idx})
                disp(['Could not load NIRS.mat for subject ' int2str(Idx)]);
            end
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
            disp(['Could not find anatomical image for subject ' int2str(Idx) ]);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROVOQUE DES MILLIARDS DE BUGS !!!!
    %Try to use field corrected image if possible
%     [dirA filA extA] = fileparts(V.fname);
%     tmp_file = fullfile(dirA,['m' filA extA]);
%     if exist(tmp_file,'file')
%         V.fname = tmp_file;
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
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
    se = strel('disk',se_size_pi);
    gaussfilt_size = job.process_image.gaussfilt_size;
    gaussfilt_sdev = job.process_image.gaussfilt_sdev;
    
    % Check if already run
    [dir1, file1, dummy] = fileparts(V.fname);
    %tmpf = spm_select('List',dir1,['_segmented_' file1]);
    %if ~isempty(tmpf)
    
    if spm_existfile(fullfile(dir1,[output_prefix,'_segmented_',file1,'.nii'])) && ~force_reprocess
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
            
            nirs_MCsegment_process_image(fullfile(dir,['c5',name,'.nii']),...
                c5method, se_size_pi, gaussfilt_size, gaussfilt_sdev);
            
            spm_progress_bar('Set',5);
            
            % gathering of layers
            % NewSegment is considered to provide a perfect segmentation
            % Control (sum of processed_ci) : according to the number of layers
            % to which belongs each
            % voxel, one determinate the right layer
            
            V5 = spm_vol(fullfile(dir,['processed_c5',name,'.nii']));Y5m = spm_read_vols(V5);clear V5;
            
            V1m = spm_vol(fullfile(dir,['c1',name,'.nii']));Y1 = spm_read_vols(V1m);clear V1m;
            V2m = spm_vol(fullfile(dir,['c2',name,'.nii']));Y2 = spm_read_vols(V2m);clear V2m;
            V3m = spm_vol(fullfile(dir,['c3',name,'.nii']));Y3 = spm_read_vols(V3m);clear V3m;
            V4m = spm_vol(fullfile(dir,['c4',name,'.nii']));Y4 = spm_read_vols(V4m);clear V4m;
            V6m = spm_vol(fullfile(dir,['c6',name,'.nii']));Y6 = spm_read_vols(V6m);clear V6m;
            
            spm_progress_bar('Set',15);
            Ysegmented = zeros(size(Y5m));
            
            % Give the number of layers each voxel belongs to (between 0 and 2)
            %%% belong2
            % Give the number of layers each voxel belongs to (between 0 and 2)
            if ~spm_existfile(fullfile(dir,filesep,'belongs2n_ci.nii'))
                
                matlabbatch{1}.spm.util.imcalc.input = {
                    fullfile(dir,['c1',name,'.nii'])
                    fullfile(dir,['c2',name,'.nii'])
                    fullfile(dir,['c3',name,'.nii'])
                    fullfile(dir,['c4',name,'.nii'])
                    fullfile(dir,['processed_c5',name,'.nii'])};
                matlabbatch{1}.spm.util.imcalc.output = 'belongs2n_ci.nii';
                matlabbatch{1}.spm.util.imcalc.outdir = {dir};
                matlabbatch{1}.spm.util.imcalc.expression = 'i1+i2+i3+i4+i5';
                matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
                matlabbatch{1}.spm.util.imcalc.options.mask = 0;
                matlabbatch{1}.spm.util.imcalc.options.interp = 1;
                matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
                spm_jobman('run',matlabbatch);
                
                Vb = spm_vol(fullfile(dir,'belongs2n_ci.nii'));Yb = spm_read_vols(Vb);
                
                 %%% Head shadow
                % is determined to know whether a voxel is situated in the head or not
                % once the sum is calculated an erosion is carried out to smoothen the
                % image and to uniformise the mask
                Yh = Yb;
%                 Y6nus = zeros(size(Yb));
                
                for i=1:size(Yh,1)
                    Yh(i,:,:) = medfilt2(squeeze(Yh(i,:,:)));
                end
                for i=1:size(Yh,2)
                    Yh(:,i,:) = medfilt2(squeeze(Yh(:,i,:)));
                    level = graythresh(squeeze(Yh(:,i,:)));
                    Yh(:,i,:) = im2bw(squeeze(Yh(:,i,:)),level);
                    Yh(:,i,:) = imfill(squeeze(Yh(:,i,:)),'holes');
                end
                nirs_create_vol(fullfile(dir,'head_shadow.nii'),...
                    Vb.dim, Vb.dt, Vb.pinfo, Vb.mat, Yh);
                
                Y6h = Y6.*Yh;

                % renormalisation
                Ysom = Y1+Y2+Y3+Y4+Y5m+Y6h;
                nirs_create_vol(fullfile(dir,'som.nii'),...
                    Vb.dim, Vb.dt, Vb.pinfo, Vb.mat, Ysom);
                
                Y1 = Y1./Ysom;
                Y2 = Y2./Ysom;
                Y3 = Y3./Ysom;
                Y4 = Y4./Ysom;
                Y5m = Y5m./Ysom;
                Y6h = Y6h./Ysom;
                nirs_create_vol(fullfile(dir,['c6h' name '.nii']),...
                    Vb.dim, Vb.dt, Vb.pinfo, Vb.mat, Y6h);

            else
                Vb = spm_vol(fullfile(dir,'belongs2n_ci.nii'));Yb = spm_read_vols(Vb);
                Vh = spm_vol(fullfile(dir,'head_shadow.nii')); Yh = spm_read_vols(Vh);
                V6h = spm_vol(fullfile(dir,['c6h',name,'.nii']));Y6h = spm_read_vols(V6h);clear V6h;
            end
            
%             Y6nus = zeros(size(Yb));
% 
%             se = strel('disk',3);
%             % pour les sinus
%             for i=1:size(Yb,1)
%                 Y6nus1(i,:,:) = imdilate(squeeze(Yb(i,:,:)),se);
%                 Y6nus1(i,:,:) = imerode(squeeze(Y6nus1(i,:,:)),se);
%                 Y6nus1(i,:,:) = medfilt2(squeeze(Y6nus1(i,:,:)));
%             end
%             for i=1:size(Yb,2)
%                 Y6nus2(:,i,:) = imdilate(squeeze(Yb(:,i,:)),se);
%                 Y6nus2(:,i,:) = imerode(squeeze(Y6nus2(:,i,:)),se);
%                 Y6nus2(:,i,:) = medfilt2(squeeze(Y6nus2(:,i,:)));
%             end
%             for i=1:size(Yb,3)
%                 Y6nus3(:,:,i) = imdilate(squeeze(Yb(:,:,i)),se);
%                 Y6nus3(:,:,i) = imerode(squeeze(Y6nus3(:,:,i)),se);
%                 Y6nus3(:,:,i) = medfilt2(squeeze(Y6nus3(:,:,i)));
%             end
%             
%             Y6 = Y6nus1+Y6nus2+Y6nus3/3;
%             Y6n = ones(size(Y6)) - Y6/max(Y6(:));
            % on a renormalise les Yi pour conserver la valeur de 1
            % belong n a plus de sens...
% % % % %             Ysegmented( Yh>0 & Y6nus<0.1) =6;
            % the voxel is attributed to the layer for which it has the
            % biggest probability with a threshold to ensure a certain
            % probability
            Ysegmented( Yh>0 & Y1>Y2  & Y1>Y3  & Y1>Y4  & Y1>Y5m & Y1>Y6h  & Y1>=thresh_as ) = 1;
            Ysegmented( Yh>0 & Y2>Y1  & Y2>Y3  & Y2>Y4  & Y2>Y5m & Y2>Y6h  & Y2>=thresh_as ) = 2;
            Ysegmented( Yh>0 & Y3>Y1  & Y3>Y2  & Y3>Y4  & Y3>Y5m & Y3>Y6h  & Y3>=thresh_as ) = 3;
            Ysegmented( Yh>0 & Y4>Y1  & Y4>Y2  & Y4>Y3  & Y4>Y5m & Y4>Y6h  & Y4>=thresh_as ) = 4;
            Ysegmented( Yh>0 & Y5m>Y1 & Y5m>Y2 & Y5m>Y3 & Y5m>Y4 & Y5m>Y6h & Y5m>=thresh_as ) = 5;
            Ysegmented( Yh>0 & Y6h>Y1 & Y6h>Y2 & Y6h>Y3 & Y6h>Y4 & Y6h>Y5m & Y6h>=thresh_as ) = 6;

            % for the remainers one use the code of rebel voxels (layer is attributed depending on the surrounding)
            jobr.rebel_surrounding = job.rebel_surrounding;
            jobr.output_prefix = output_prefix;
            jobr.name = name;
            jobr.dir = dir;
            jobr.Yseg = Ysegmented;
            jobr.Yb = fullfile(dir,'belongs2n_ci.nii');
            jobr.Yh = fullfile(dir,'head_shadow.nii');
            outr = nirs_MCsegment_rebels2(jobr);
            
            % Save name of segmented image in NIRS matrix
            NIRS.Dt.ana.T1seg = outr;
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
