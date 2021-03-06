function out = nirs_create_roi_run(job)
% Batch function to create ROIs.
%_______________________________________________________________________________
% Copyright (C) 2013 LIOM Laboratoire d'Imagerie Optique et Moleculaire
%                    Ecole Polytechnique de Montreal
%_______________________________________________________________________________

% ------------------------------------------------------------------------------
% REMOVE AFTER FINISHING THE FUNCTION //EGC
% ------------------------------------------------------------------------------
% fprintf('Work in progress...\nEGC\n')
% out.NIRSmat = job.NIRSmat;
% return
% ------------------------------------------------------------------------------
try
    SelectPreviousROI = job.SelectPreviousROI;
catch
    SelectPreviousROI = false;
end
% Manual/automatic selection of ROIs/seeds (pointNclickROI ManualROIspline)
fNamesROIchoice         = fieldnames(job.AutoROIchoice);
autoROI                 = false;
graphicalROI            = false;
ManualROIspline         = false;
pointNclickROI          = false;
pointNclickSquare       = false;
switch(fNamesROIchoice{1})
    case 'AutoROI'
        autoROI         = 1;
        ROIsize         = job.AutoROIchoice.AutoROI.ArrayROI;
    case 'ManualROI'
        graphicalROI    = 1;
    case 'ManualROIspline'
        ManualROIspline = 1;
    case 'pointNclickROI'
        pointNclickROI  = 1;
        radiusX         = job.AutoROIchoice.pointNclickROI.ManualROIradius;
        radiusY         = job.AutoROIchoice.pointNclickROI.ManualROIradius;
    case 'pointNclickSquare'
        pointNclickSquare  = 1;
%         ManualROIwidth  = job.AutoROIchoice.pointNclickSquare.ManualROIwidth;
%         ManualROIheight = job.AutoROIchoice.pointNclickSquare.ManualROIheight;
        radiusX         = job.AutoROIchoice.pointNclickSquare.ManualROIwidth;
        radiusY         = job.AutoROIchoice.pointNclickSquare.ManualROIheight;
    case 'ManualEnterROI'
        % Do nothing
    otherwise
        % Do nothing
end

for scanIdx=1:length(job.NIRSmat)
    try
        %Load NIRS.mat information
        [NIRS NIRSmat dir_nirsmat]= pat_get_PATmat(job,scanIdx);
        [~, ~, ~, ~, ~, ~, splitStr] = regexp(NIRS.input_dir,'\\');
        scanName = splitStr{end-1};
        if pointNclickROI
            % radius in pixels
            radiusX = job.AutoROIchoice.pointNclickROI.ManualROIradius/NIRS.PAparam.pixWidth;
            radiusY = job.AutoROIchoice.pointNclickROI.ManualROIradius/NIRS.PAparam.pixDepth;
        end
        if pointNclickSquare
            % radius in pixels
            radiusX = job.AutoROIchoice.pointNclickSquare.ManualROIwidth/NIRS.PAparam.pixWidth;
            radiusY = job.AutoROIchoice.pointNclickSquare.ManualROIheight/NIRS.PAparam.pixDepth;
        end
        if ~isfield(NIRS.flags,'ROIOK') || job.force_redo
            if job.RemovePreviousROI
                try
                    NIRS.res = rmfield(NIRS.flags,'ROIOK');
                end
                try
                    for i1=1:length(NIRS.res.NIRS)
                        %clean up: delete ROI mask files
                        delete(NIRS.res.ROI{i1}.fname);
                    end
                end
                try
                    NIRS.res = rmfield(NIRS.res,'ROI');
                end
                index = 0;
            else
                if isfield(NIRS.res,'ROI')
                    index = length(NIRS.res.ROI);
                else
                    index = 0;
                end
            end
            
            % Display anatomical image
            try
                vol_anat = spm_vol(NIRS.res.file_anat);
            catch
                disp('Could not find anatomical image');
                [t sts] = spm_select(1,'image','Select anatomical image','',dir_nirsmat,'.*',1);
                NIRS.res.file_anat = t;
                vol_anat = spm_vol(NIRS.res.file_anat);
            end
            [dir1 fil1] = fileparts(vol_anat.fname);
            im_anat = spm_read_vols(vol_anat);

            
            spm_figure('GetWin', 'Graphics');
            spm_figure('Clear', 'Graphics');
            
            if isfield(job,'displayBrainmask')
                if job.displayBrainmask == 1 && isfield(NIRS,'fcPAT') && isfield(NIRS.fcPAT,'mask') && isfield(NIRS.fcNIRS.mask,'fname')
                    % Display only brain pixels mask
                    vol = spm_vol(NIRS.fcNIRS.mask.fname);
                    full_mask = logical(spm_read_vols(vol));
                else
                    % Display all the image
                    full_mask = ones(size(im_anat));
                    if job.displayBrainmask
                        disp('Could not find brain mask')
                    end
                end
            else
                full_mask = ones(size(im_anat));
            end
            if SelectPreviousROI
                % Goto interactive window
                h2 = spm_figure('GetWin', 'Interactive');
                % Ypos = 1, always first line
                spm_input(['Subject ' int2str(scanIdx) ' ' scanName], 1, 'd');
                SelPrevROI = spm_input('Select a previous list of ROIs?','+1','y/n');
                if SelPrevROI == 'y'
                    [tPrevROI stsPrevROI] = spm_select(1,'mat','Select NIRS.mat structure containing information on desired ROIs','',dir_nirsmat,'NIRS.mat',1);
                    if stsPrevROI
                        PAT0 = NIRS; %Store current NIRS
                        try
                            load(tPrevROI);
                            PAT_withROIs = NIRS;
                            NIRS = PAT0;
                            try
                                if ~isfield(NIRS.res,'ROI')
                                    NIRS.ROI.ROIname = PAT_withROIs.ROI.ROIname;
                                    NIRS.res.ROI = PAT_withROIs.res.ROI;
                                else
                                    NIRS.ROI.ROIname = [NIRS.ROI.ROIname; PAT_withROIs.ROI.ROIname];
                                    NIRS.res.ROI = [NIRS.res.ROI PAT_withROIs.res.ROI];
                                end
                                index = index+length(NIRS.res.ROI);
                                clear PAT_withROIs
                            catch
                                NIRS = PAT0;
                                disp('Specified NIRS.mat structure does not contain valid ROI information')
                            end
                        catch
                            disp('Could not load NIRS.mat structure containing desired ROIs')
                        end
                    end
                end
            end
            
            if ~autoROI
                %Display images of changes from 10th to 90th percentile for all sessions
                try % <--- Needs a catch statement //EGC
                    nCols = ceil(sqrt(length(NIRS.sess_res)));
                    nRows = ceil(length(NIRS.sess_res) / nRows);
                    for i0=1:length(NIRS.sess_res)
                        % Only open 1 figure for all the sessions
                        if i0 == 1,
                            hs = figure;
                            set(hs, 'Name', '10-90 percentile changes')
                        end
                        V = spm_vol(NIRS.sess_res{i0}.fname_change_90_10{1}); %color green
                        tmp_image = spm_read_vols(V);
                        figure(hs);
                        subplot(nRows, nCols, i0);
                        imagesc(tmp_image); 
                        % axis image
                        set(gca,'DataAspectRatio',[1 NIRS.PAparam.pixWidth/NIRS.PAparam.pixDepth 1])
                        title(['Session ' int2str(i0) ': ratio of 90th to 10th percentile']);
                    end
                end
                oneMoreROI = 1;
                % Goto interactive window
                h2 = spm_figure('GetWin', 'Interactive');
                
                % Ypos = 1, always first line
                spm_input(['Subject ' int2str(scanIdx) ' ' scanName], 1, 'd');
                linecount = 0;
                while oneMoreROI
                    figure(h2);
                    oneMoreROI = spm_input('Add an ROI?',2*linecount+2,'y/n');
                    if oneMoreROI == 'y', oneMoreROI = 1; else oneMoreROI = 0; end
                    if oneMoreROI
                        % h1 = figure('Position',[20 50 3*size(im_anat,1) 3*size(im_anat,2)]);
                        % Display anatomical image on SPM graphics window
                        minVal = min(im_anat(:));
                        maxVal = max(im_anat(:));
                        spm_figure('GetWin', 'Graphics');
                        spm_figure('Clear', 'Graphics');
                        imagesc(im_anat .* full_mask, [minVal, maxVal]);
                                    % Gray-scale colormap to enhance contrast
                        if useGrayContrast
                            cmap = contrast(im_anat);
                            colormap(cmap);
                        else
                            colormap(gray);
                        end
                        % axis image;
                        set(gca,'DataAspectRatio',[1 NIRS.PAparam.pixWidth/NIRS.PAparam.pixDepth 1])
                        if graphicalROI
                            % Specify polygonal region of interest (ROI)
                            title(sprintf('Make ROI polygon, then double click in it to create ROI (scan %d of %d).', scanIdx, length(job.NIRSmat)));
                            mask = roipoly;
                        elseif ManualROIspline
                            % Start interactive ROI tool to choose spline
                            % ROI/seed
                            mask = pat_roi_spline(im_anat,[],[],sprintf('Mark spline points, then right-click in it to create ROI/seed (scan %d of %d)',scanIdx, length(job.NIRSmat)));
                        else
                            if pointNclickROI
                                % Specify center of circular ROI/seed with mouse
                                % point & click on the anatomical image
                                title(sprintf('Click the center of circular ROI/seed (scan %d of %d)', scanIdx, length(job.NIRSmat)))
                                % Circular seed setup
                                t = 0:pi/100:2*pi;
                                % Prompt user to point & click
                                p = ginput(1);
                                % Center of the seed coordinates
                                x0 = p(1); y0 = p(2);
                                
                                % Parametric function for a circle
                                xi = radiusX * cos(t) + x0;
                                yi = radiusY * sin(t) + y0;
                                LineHandler = line(xi,yi,'LineWidth',3,'Color',[.8 0 0]);
                                % Create ROI/seed mask
                                mask = poly2mask(xi, yi, size(im_anat,1), size(im_anat,2));
                                % Save coordinates of seed for later display.
                                % NOTE: Row is 1st coordinate, Column is 2nd
                                NIRS.res.ROI{index+1}.center = [y0 x0];
                                NIRS.res.ROI{index+1}.radius = job.AutoROIchoice.pointNclickROI.ManualROIradius;
                            else
                                if pointNclickSquare
                                    % Specify center of a square ROI/seed with mouse
                                    % point & click on the anatomical image
                                    title(sprintf('Click the center of rectangular ROI/seed (scan %d of %d)', scanIdx, length(job.NIRSmat)))
                                    % Square setup
%                                     hR = imrect(gca,[0 0 ManualROIwidth ManualROIheight]);
                                    hR = imrect(gca,[0 0 radiusX radiusY]);
                                    pR = wait(hR);
                                    xi = [pR(1) pR(1)+pR(3) pR(1)+pR(3) pR(1)];
                                    yi = [pR(2) pR(2)       pR(2)+pR(4) pR(2)+pR(4)];
                                    pRr = round(pR);
                                    disp(['Width = ' int2str(pRr(3)) ', Height = ' int2str(pRr(4)) ...
                                        ', Center(x) = ' int2str(round(pRr(1)+pRr(3)/2)) ...
                                        ', Center(y) = ' int2str(round(pRr(2)+pRr(4)/2))])
                                    % Create the mask
                                    mask = poly2mask(xi, yi, size(im_anat,1), size(im_anat,2));
                                  
                                    % Save coordinates of seed for later display.
                                    % NOTE:  [xmin ymin width height]
                                    NIRS.res.ROI{index+1}.center = [round(pRr(1)+pRr(3)/2) round(pRr(2)+pRr(4)/2)];
                                    NIRS.res.ROI{index+1}.width_height = [pRr(3) pRr(4)];
                                    NIRS.res.ROI{index+1}.width_height_mm = [job.AutoROIchoice.pointNclickSquare.ManualROIwidth job.AutoROIchoice.pointNclickSquare.ManualROIheight];
                                else
                                    % Manual ROI coordinate entry
                                    linecount = linecount + 1;
                                    rc = spm_input('Enter [row,column] of center',2*linecount+1,'e',[],2);
                                    linecount = linecount + 1;
                                    radius = spm_input('Enter radius in pixels',2*linecount+1,'e',[],1);
                                    radius = round(radius);
                                    if radius < 0, radius = 0; end
                                    mask = zeros(size(im_anat));
                                    for x1=-radius:radius
                                        for y1=-radius:radius
                                            if x1^2+y1^2 <= radius^2
                                                try %will skip pixels outside the image
                                                    mask(rc(1)+y1,rc(2)+x1) = 1;
                                                end
                                            end
                                        end
                                    end
                                    % Save coordinates of seed for later display
                                    % NOTE: Row is 1st coordinate, Column is 2nd
                                    NIRS.res.ROI{index+1}.center = rc;
                                    NIRS.res.ROI{index+1}.radius = radius;
                                end
                            end
                        end
                        
                        mask = single(mask);
                        % Update ROI's display
                        full_mask = full_mask - mask;
                        index = index + 1; linecount = linecount + 1;
                        if job.select_names
                            figure(h2);
                            name = spm_input(['Enter name of ROI' int2str(index)],2*linecount+1,'s');
                        else
                            name = int2str(index);
                        end
                        NIRS.res.ROI{index}.name = name;
                        if index < 10, str0 = '0'; else str0 = ''; end
                        str = [str0 int2str(index)];
                        % Save nifti files in ROI sub-folder //EGC
                        fname_mask = fullfile(dir_nirsmat,[fil1 '_ROI_' str '_' name '.nii']);
                        NIRS.res.ROI{index}.fname = fname_mask;
                        
                        % Correct data type for logical masks
                        switch (vol_anat(1).dt(1))
                            case 512
                                % unsigned integer 16-bit
                                vol_anat(1).dt = [64 0];
                            case 2
                                % unsigned integer 8-bit
                                vol_anat(1).dt = [64 0];
                            case 64
                                % Float 64 poses no problem
                            otherwise
                                % Convert to unsigned int 8-bit
                                vol_anat(1).dt = [64 0];
                        end

                        % Create and write a NIFTI file with ROI mask
                        pat_create_vol(fname_mask, vol_anat(1).dim, vol_anat(1).dt,...
                            vol_anat(1).pinfo, vol_anat(1).mat, 1, mask);
                    end
                end
                try close(hs); end
                
            else
                %automatic ROIs
                %h1 = figure('Position',[20 50 3*size(im_anat,1) 3*size(im_anat,2)]);
                sz = size(im_anat);
                sz = sz(1:2); %remove 3rd component, which is one
                sz = floor(sz./ROIsize);
                for i1 = 1:ROIsize(1) %N
                    for i2 = 1:ROIsize(2) %M
                        %ROI created in following order
                        % 1      2       3 ...   ROIsize(2)=M
                        % M+1   M+2     M+3 ...  2*M
                        %  .    .        .       .
                        %  .    .        .       N*M
                        index = index+1;
                        %vertices in the order: UL, UR, LR, LL
                        mask = roipoly(im_anat,...
                            [1+sz(2)*(i2-1) 1+sz(2)*i2 1+sz(2)*i2 1+sz(2)*(i2-1)],...
                            [1+sz(1)*(i1-1) 1+sz(1)*(i1-1) 1+sz(1)*i1 1+sz(1)*i1]);
                        mask = single(mask);
                        if index < 10, str0 = '0'; else str0 = ''; end
                        str = [str0 int2str(index)];
                        % Save nifti files in ROI sub-folder //EGC
                        fname_mask = fullfile(dir_nirsmat,[fil1 '_ROI_' str '_' int2str(i1) 'x' int2str(i2) '.nii']);
                        NIRS.res.ROI{index}.fname = fname_mask;
                        
                        % Correct data type for logical masks
                        switch (vol_anat(1).dt(1))
                            case 512
                                % unsigned integer 16-bit
                                vol_anat(1).dt = [64 0];
                            case 2
                                % unsigned integer 8-bit
                                vol_anat(1).dt = [64 0];
                            case 64
                                % Float 64 poses no problem
                            otherwise
                                % Convert to unsigned int 8-bit
                                vol_anat(1).dt = [64 0];
                        end
                        
                        % Create and write a NIFTI file with ROI mask
                        pat_create_vol(fname_mask, vol_anat(1).dim, vol_anat(1).dt,...
                            vol_anat(1).pinfo, vol_anat(1).mat, 1, mask);
                    end
                end
            end
            NIRS.ROI.ROIname = {};
            for i0=1:length(NIRS.res.ROI)
                if isfield(NIRS.res.ROI{i0},'name')
                    NIRS.ROI.ROIname = [NIRS.ROI.ROIname; NIRS.res.ROI{i0}.name];
                else
                    NIRS.ROI.ROIname = [NIRS.ROI.ROIname; ['ROI' gen_num_str(i0,3)]];
                end
            end
            % ROI creation succesful
            NIRS.flags.ROIOK = true;
            save(NIRSmat,'NIRS');
        end
        fprintf('Scan %s, %d of %d complete\n', scanName, scanIdx, length(job.NIRSmat));
        out.NIRSmat{scanIdx} = NIRSmat;
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.NIRSmat{scanIdx} = job.NIRSmat{scanIdx};
    end % End try
end % Scans loop
end % End function
% EOF
