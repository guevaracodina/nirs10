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
            disp(['Could not find anatomical image for ' int2str(Idx) 'th subject.']);
        end
    end
    
    
    % USER OPTIONS %
    thresh_as = job.thresh_as;
    
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
%         try
            [dir,name] = fileparts(V.fname);
            
            %%% on traite les couches
            nirs_MCsegment_process_image(fullfile(dir,['c5',name,'.nii']),...
                c5method, se_size_pi, gaussfilt_size, gaussfilt_sdev);
            
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
                spm_jobman('run_nogui',matlabbatch);
                
%                 [dummy,Vb] = spm_imcalc_ui([fullfile(dir,['c5',name,'.nii']);...
%                     fullfile(dir,['c2',name,'.nii']);...
%                     fullfile(dir,['c3',name,'.nii']);...
%                     fullfile(dir,['c4',name,'.nii']);...
%                     fullfile(dir,['c1',name,'.nii'])],...
%                     fullfile(dir,'belongs2n_ci.nii'),...
%                     'i1+i2+i3+i4+i5');
                
                %%% Head shadow
                % is determined to know whether a voxel is situated in the head or not
                % once the sum is calculated an erosion is carried out to smoothen the
                % image and to uniformise the mask
                Vb = spm_vol(fullfile(dir,'belongs2n_ci.nii'));
                Yb = spm_read_vols(Vb);
                Yh = Yb;
                
                for i=1:size(Yh,1)
                    Yh(i,:,:) = medfilt2(squeeze(Yh(i,:,:)));
                end
                for i=1:size(Yh,2)
                    Yh(:,i,:) = medfilt2(squeeze(Yh(:,i,:)));
                    level = graythresh(squeeze(Yh(:,i,:)));
                    Yh(:,i,:) = im2bw(squeeze(Yh(:,i,:)),level);
                    Yh(:,i,:) = imfill(squeeze(Yh(:,i,:)),'holes');
                end
                
                Vh = struct('fname',fullfile(dir,'head_shadow.nii'),...
                    'dim',  Vb.dim,...
                    'dt',   Vb.dt,...
                    'pinfo',Vb.pinfo,...
                    'mat',Vb.mat);
                Vh = spm_create_vol(Vh);
                spm_write_vol(Vh, Yh);
            else
                Vb = spm_vol(fullfile(dir,'belongs2n_ci.nii'));Yb = spm_read_vols(Vb);
                Vh = spm_vol(fullfile(dir,'head_shadow.nii')); Yh = spm_read_vols(Vh);
            end
            
            Ysegmented = zeros(size(Yb));
            
            V1 = spm_vol(fullfile(dir,['c1',name,'.nii']));Y1 = spm_read_vols(V1);clear V1;
            V2 = spm_vol(fullfile(dir,['c2',name,'.nii']));Y2 = spm_read_vols(V2);clear V2;
            V3 = spm_vol(fullfile(dir,['c3',name,'.nii']));Y3 = spm_read_vols(V3);clear V3;
            V4 = spm_vol(fullfile(dir,['c4',name,'.nii']));Y4 = spm_read_vols(V4);clear V4;
            V5 = spm_vol(fullfile(dir,['processed_c5',name,'.nii']));Y5 = spm_read_vols(V5);
            
            %%% choosing layer...
            %--> belongs_2n==1 : voxel is directly attributed
            Ysegmented( Yb>=0.5 & Yb<1.51 & Y5>=thresh_as ) = 5;
            Ysegmented( Yb>=0.5 & Yb<1.51 & Y4>=thresh_as ) = 4;
            
            Ysegmented( Yb>=0.5 & Yb<1.51 & Y2>=thresh_as ) = 2;
            Ysegmented( Yb>=0.5 & Yb<1.51 & Y1>=thresh_as ) = 1;
            Ysegmented( Yb>=0.5 & Yb<1.51 & Y3>=thresh_as ) = 3;
            
            spm_progress_bar('Set',20);
            %--> belongs_2n==2 : voxel is attributed according to the biggest
            % probability (of the random variable)
            Ysegmented( Yb>=1.51 & Y4>=Y1 & Y4>=Y2 & Y4>=Y3 & Y4>=Y5 ) = 4;
            Ysegmented( Yb>=1.51 & Y5>=Y1 & Y5>=Y2 & Y5>=Y3 & Y5>=Y4 ) = 5;
            
            Ysegmented( Yb>=1.51 & Y2>=Y1 & Y2>=Y3 & Y2>=Y4 & Y2>=Y5 ) = 2;
            Ysegmented( Yb>=1.51 & Y1>=Y2 & Y1>=Y3 & Y1>=Y4 & Y1>=Y5 ) = 1;
            Ysegmented( Yb>=1.51 & Y3>=Y1 & Y3>=Y2 & Y3>=Y4 & Y3>=Y5 ) = 3;
            spm_progress_bar('Set',25);
            %--> belongs_2n==0 : voxel has'nt been attributed to any
            %layer...
            % All those belonging to the head shadow are attributed depending
            % on the values of random variables
            Ysegmented( Yb<0.5 & Yh==1 & Y4>=Y1 & Y4>=Y2 & Y4>=Y3 & Y4>=Y5 & Y4~=0 ) = 4;
            Ysegmented( Yb<0.5 & Yh==1 & Y5>=Y1 & Y5>=Y2 & Y5>=Y3 & Y5>=Y4 & Y5~=0 ) = 5;
            Ysegmented( Yb<0.5 & Yh==1 & Y2>=Y1 & Y2>=Y3 & Y2>=Y4 & Y2>=Y5 & Y2~=0 ) = 2;
            Ysegmented( Yb<0.5 & Yh==1 & Y1>=Y2 & Y1>=Y3 & Y1>=Y4 & Y1>=Y5 & Y1~=0 ) = 1;
            Ysegmented( Yb<0.5 & Yh==1 & Y3>=Y1 & Y3>=Y2 & Y3>=Y4 & Y3>=Y5 & Y3~=0 ) = 3;
            
            jobr.rebel_surrounding = job.rebel_surrounding;
            jobr.output_prefix = output_prefix;
            jobr.name = name;
            jobr.dir = dir;
            jobr.Yseg = Ysegmented;
            jobr.Yb = fullfile(dir,'belongs2n_ci.nii');
            jobr.Yh = fullfile(dir,'head_shadow.nii');
            
            nirs_MCsegment_process_image(fullfile(dir,['c1',name,'.nii']),...
                c1method, se_size_pi, gaussfilt_size, gaussfilt_sdev);
            nirs_MCsegment_process_image(fullfile(dir,['c2',name,'.nii']),...
                c2method, se_size_pi, gaussfilt_size, gaussfilt_sdev);
            nirs_MCsegment_process_image(fullfile(dir,['c3',name,'.nii']),...
                c3method, se_size_pi, gaussfilt_size, gaussfilt_sdev);
            nirs_MCsegment_process_image(fullfile(dir,['c4',name,'.nii']),...
                c4method, se_size_pi, gaussfilt_size, gaussfilt_sdev);
            
            jobr.Y1 = fullfile(dir,['processed_c1',name,'.nii']);
            jobr.Y2 = fullfile(dir,['processed_c2',name,'.nii']);
            jobr.Y3 = fullfile(dir,['processed_c3',name,'.nii']);
            jobr.Y4 = fullfile(dir,['processed_c4',name,'.nii']);
            jobr.Y5 = fullfile(dir,['processed_c5',name,'.nii']);
            outr = nirs_MCsegment_rebels(jobr);
                        
            % Save name of segmented image in NIRS matrix
            NIRS.Dt.ana.T1seg = outr;
            if NIRSok
                save(outNIRSmat{Idx},'NIRS');
            else
                %what to do if we don't have a NIRS structure already?
            end
            
%         catch
%             disp(['MC Segment failed to run for subject ' int2str(Idx)]);
%         end   
    end
end
out.NIRSmat = outNIRSmat;
