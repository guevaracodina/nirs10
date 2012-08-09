function out = nirs_run_detectVitamins(job)
% This module is used for semi-automatically detecting fiducial markers 
% (e.g. marking the positions of the optodes) on an MRI anatomical image.
% The positions are detected and saved in the NIRS.mat matrix, and the 
% the markers are erased from the anatomical image (a copy of the original
% image is saved). The name of this new "cleaned" anatomical image is also
% saved in the nirs matrix in place of the original image.
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
%
% Michèle Desjardins
% 2011-03

% Boolean for "semi-manual"/debug mode
manualMode = 1;
% Boolean for when the semi-automatic coregistration failed and one just
% wants to manually define the order of the optodes
wrongWhenAutomatic = 1;
% Boolean for an exception... always leave =0 except for P2S24 & P2S27.
% Sorry about that.
needFixForMissingVit = 0;

% INPUTS %
%%%%%%%%%%

% NIRS matrices for all subjects
NIRSmat = job.NIRSmat;
%load(NIRSmat{:});

% (Optional) Path to anatomical image for all subjects
%job.anatT1; -> is read for each subject in the loop below

% Prefix to be added to the name of anatT1 file when saving image with
% markers removed
output_prefix_woVit = job.output_prefix_woVit;


% OUTPUTS %
%%%%%%%%%%
outNIRSmat = {};



% Number of subjects
nSubj = size(job.NIRSmat,1);

% Loop over subjects
for iSubj = 1:nSubj
              
    [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{iSubj,1},job.NIRSmatCopyChoice,job.force_redo);
    job.NIRSmat{iSubj,1} = newNIRSlocation;
    if ~isempty(NIRS) && (~isfield(NIRS,'flags') || ~isfield(NIRS.flags,'detVit_OK') || job.force_redo)

    % Exeception.......
    temp = needFixForMissingVit;
    theNIRSmat = load(NIRSmat{iSubj});
    if strcmp(theNIRSmat.NIRS.Dt.s.p(end-3:end),'S024') || ...
            strcmp(theNIRSmat.NIRS.Dt.s.p(end-3:end),'S027')
            %iSubj==24 || iSubj==27
        needFixForMissingVit = 1;
    else
        needFixForMissingVit = temp;
    end
    
    if needFixForMissingVit % Found in Michèle's stuff... do not use except to fix
        % the bug for 2 particular subjects of study MDEIEP2 (P2S24 &
        % P2S27)
        fix_coreg_missing_vitamin(fileparts(NIRSmat{iSubj}),theNIRSmat.NIRS.Dt.s.p(end-3:end));
        
    % Use the normal code    
    else
    
    temp = wrongWhenAutomatic;
    if strcmp(theNIRSmat.NIRS.Dt.s.p(end-3:end),'S021') || ...
            strcmp(theNIRSmat.NIRS.Dt.s.p(end-3:end),'S040') || ...
            strcmp(theNIRSmat.NIRS.Dt.s.p(end-3:end),'S041') || ...
            strcmp(theNIRSmat.NIRS.Dt.s.p(end-3:end),'S043') || ...
            strcmp(theNIRSmat.NIRS.Dt.s.p(end-3:end),'S044')
    %if iSubj==21 || ...% iSubj==24 || iSubj==27 || ...
    %        iSubj==40 || iSubj==41 || iSubj==43 || iSubj==44
        wrongWhenAutomatic = 1;

    else
        wrongWhenAutomatic = temp;
    end
    
    clear NIRS Vanat Yanat
    
    % Read anatomical image and NIRS matrix
    try    
        load(NIRSmat{iSubj});
        if isempty(job.anatT1{1,1})
            try
                image_in = NIRS.Dt.ana.T1;
            catch
                disp('Could not find an anatomical image');
            end
        else
            % Store T1 file location
            NIRS.Dt.ana.T1 = job.anatT1{iSubj,:};
            image_in = job.anatT1{iSubj,:};
        end
        Vanat = spm_vol(image_in);
        Yanat = spm_read_vols(Vanat);
    catch
        disp(['Could not read NIRS.mat or anatomical image for ' int2str(iSubj) 'th subject.' ]);
    end
    
    subjPath = NIRS.Dt.s.p;
    if ~exist(fullfile(subjPath,'spots.mat'),'file')
    
        % Define structuring elements
        se =  strel('ball',5,5);
        se_larger = strel('ball',7,7);

        % Fill head image (paint-bucket)
        filledYanat = [];
        for i=1:size(Yanat,3)
            filledYanat(:,:,i) = imfill(Yanat(:,:,i),'holes');
        end

        % Erode, then dilate with ball structuring element. This operation will
        % preserve objects larger than the structuring element (with unmodyfied
        % contours) but will have erased objects smaller than it.
        erodedBW = imerode(filledYanat,se);
        dilatedBW = imdilate(erodedBW,se);

        % Close with a larger strucring element to fill in some gaps (is this
        % step necessary?)
        closedBW = imclose(dilatedBW,se_larger);
        % Fill remaining holes in the head
        filledBW = [];
        for i=1:size(closedBW,3)
            filledBW(:,:,i) = imfill(closedBW(:,:,i),'holes');
        end

        % We are now left with an image where the head is intact, but the
        % smaller (fiducials & noise) features of the image have been erased.

        % Transform to a binary (mask) image, using simple thresholding
        % (improve this method?? Global thresholding using Otsu's method
        % (greythresh) did not work here.

        thresh=0.05;
        mask = filledBW>(thresh.*max(filledBW(:)));
        % This should be a 0/1 mask of the head with the contours unchanged 
        % from those of the original iamge.

        % Original image without fiducials - will be saved
        Yanat_WOspots = Yanat.*mask;
        % Image of the fiducials (with background noise) - will be used to
        % detect the fiducials' positions
        spots = Yanat - Yanat_WOspots;

        save(fullfile(subjPath,'spots.mat'),'spots');
    
    else
        load(fullfile(subjPath,'spots.mat'));
        Yanat_WOspots = Yanat - spots;
    end
    
    if ~exist(fullfile(subjPath,'spots.nii'),'file')
        % Save also as .nii
        Vspots = Vanat;
        Vspots.fname = fullfile(subjPath,'spots.nii');
        spm_write_vol(Vspots,spots);
    end
    
    % In order to detect the fiducials, we want to isolate them from the
    % background noise. This operation will select all regional maxima
    % higher than a certain threshold (defined in % of the max intensity)
    % and 0-out everything else in the image.
    thresh=0.5;
    vitaminsBW = imextendedmax(spots,thresh*max(spots(:)),26);
    % Connectivity 26, we are working in 3D.
    % The resulting image is a binary mask where the fiducials are 1s and
    % the background is 0. All that is left is to detect their positions.

    % Identify centroid of the objects in the image (fiducials) - thank you
    % image processing toolbox!
    stats = regionprops(vitaminsBW, 'Centroid', 'Area');
    coord_fid = [];
    size_fid = [];
    for iFiducial=1:size(stats,1)
        coord_fid(iFiducial,:) = stats(iFiducial).Centroid(:);
        size_fid(iFiducial,:) = stats(iFiducial).Area(:);
    end
    
    if manualMode
        % Display for manual corrections...
        figure(101), plot3(coord_fid(:,1),coord_fid(:,2),coord_fid(:,3),'or')
        %addpath(genpath('G:\MesProgrammes\tMCimg_scripts\Scripts\PourPreparerSimuMC\New_SPMnewSegment'));
        %optodesCoord_click(Yanat,'Verify correspondance of the optodes',...
        %   round(coord_fid(:,[2 1 3])))
    end
    
    % Read helmet information
    
    % Either use helemt info (src & det positions) in NIRS matrix or use
    % template BrainSight-coregistered info
    % For now, we will use nirs file info.
     nSrc = NIRS.Cf.H.S.N;
     nDet = NIRS.Cf.H.D.N;
     coord_helmet = [NIRS.Cf.H.S.r.o.mm.p; NIRS.Cf.H.D.r.o.mm.p]; % nOptodes x 3, mm
                
    
    % We now have the coordinates of the centroïds of the fiducials. In
    % order to identify each source and detector and their names, as well
    % as correct for potential errors (fake markers due to noise mistaken
    % for a fiducial), we will compare with the information given in the
    % NIRS matrix regarding the helmet used for the acquisition.
    % Thus we need to "match" coord_helmet and coord_fid by reordering the
    % rows of coord_fid. (Both are of size nOptodesx3).
   
    % If there is a "fake" marker (remaining noise mistaken for a
    % fiducial), chances are it will be smaller or larger than the rest
    nFakes = size(coord_fid,1) - size(coord_helmet,1);
    toKeep = 1:size(coord_fid,1);
    if nFakes > 0 % Arterfactual fiducials
        dSize = abs(size_fid - mean(size_fid));
        [dSize_sort,idx] = sort(dSize,'descend');
        toKeep = setdiff(idx,idx(1:nFakes));
        % Remove the n fiducials that have the largest size
        % difference relative to the mean size
            
    elseif nFakes < 0 % Missing fiducials...?
        % What to do??
    end
    
    
    try % optodes previously manaully corrected and saved
        load(fullfile(subjPath,'outCorrected.mat'),'outCorrected');
        coord_fid = outCorrected(:,[2 1 3]);
    catch 
        if manualMode % manaually correct positions that were detected automatically
           idxSubjID = strfind(subjPath,'P2S');
           subjID = subjPath(idxSubjID:idxSubjID+4);
           outCorrected = optodesCoord_click(Yanat,[subjID '. Left: Keep, Right: Reject'],...
               round(coord_fid(toKeep,[2 1 3])),round(coord_fid(setdiff(idx,toKeep),[2 1 3])));
           if ~isempty(outCorrected) % "OK"
              coord_fid = outCorrected(:,[2 1 3]);
              save(fullfile(subjPath,'outCorrected.mat'),'outCorrected');
           else % "cancel"
               coord_fid = coord_fid(toKeep,:);
           end
           
        else % optodes detected automatically
            coord_fid = coord_fid(toKeep,:);
        end
    end
    
    nOptodes = size(coord_fid,1);
    
    % Project coorde_fid on a (well-chosen) 2D plane
    
    % Find plane z = a1*x+a2*y+a0 that is closest to the points
    Z = coord_fid(:,3);
    X = [coord_fid(:,1) coord_fid(:,2) ones(nOptodes,1)];
    aHat = inv(X'*X)*X'*Z; % [a1 a2 a0]
    a1 = aHat(1);
    a2 = aHat(2);
    a0 = aHat(3);

    % Project orthogonally on that plane
    % 2 vectors that form an ortogonal basis for the plane are
    u1 = [1; 0; a1]; % 3x1
    u2 = [-a1*a2; 1+a1^2; a2];
    u1 = u1./norm(u1);
    u2 = -u2./norm(u2);
    % (It can be shown that these 2 vectors belong to the plane and that they
    % are not only linearly independent, but also orthogonal)
    % The projection of the coordinates on this plane is given by the (vector)
    % sum of the projection on each of the vectors of the orthogonal basis,
    % i.e. Proj = (coord*basisVector/(basisVector'*basisVector)) .* basis vector
    coord_fid_proj2D = [];
    for iOptode = 1:nOptodes
        coord_fid_proj2D(iOptode,:) =  [ (coord_fid(iOptode,:)*u1)./(u1'*u1) ...
                                       (coord_fid(iOptode,:)*u2)./(u2'*u2) ];
    end

    % Center on 0 in u1-u2 plane
    ctr_mass_fid = [mean(coord_fid_proj2D(:,1)); mean(coord_fid_proj2D(:,2))];
    coord_fid_proj2D = coord_fid_proj2D - (ctr_mass_fid*ones(1,size(coord_fid_proj2D,1)))'; 
    
    if manualMode
        % Display for manual corrections...
        figure(102), plot(coord_fid_proj2D(:,1),coord_fid_proj2D(:,2),'or')
    end
    
    % Center helmet on 0 and use 2D coordinates
    ctr_mass_helmet = [mean(coord_helmet(:,1)); mean(coord_helmet(:,2)); mean(coord_helmet(:,3))];
    coord_helmet = coord_helmet - (ctr_mass_helmet*ones(1,size(coord_helmet,1)))'; 
    coord_helmet2D = coord_helmet(:,1:2);
    
    % We must now match 2 sets of 2D coordinates, which will require to
    % first find the right orientation, using some manipulations of the
    % images of the coordinates
    
    % Fix for when semi-automatic coregistration failed...
    if wrongWhenAutomatic
        if strcmp(theNIRSmat.NIRS.Dt.s.p(end-3:end),'S040') 
            optOrder = [13 9 6 2 11 10 7 1 12 8 5 4 3];
        elseif strcmp(theNIRSmat.NIRS.Dt.s.p(end-3:end),'S021') 
            optOrder = [12 8 5 1 10 9 6 2 11 7 4 3];
        elseif strcmp(theNIRSmat.NIRS.Dt.s.p(end-3:end),'S041') 
            optOrder = [13 9 7 4 11 10 6 1 12 8 5 2 3];
        elseif strcmp(theNIRSmat.NIRS.Dt.s.p(end-3:end),'S043') 
            optOrder = [13 9 6 2 11 10 7 1 12 8 5 4 3];
        elseif strcmp(theNIRSmat.NIRS.Dt.s.p(end-3:end),'S044') 
            optOrder = [13 9 6 2 12 10 7 3 11 8 5 4 1];
        else 
        
            % Then we will plot the projected 2D geometry and prompt the user
            % to enter the optimal order in the command window... (lazier and
            % less clean than a GUI but this situation should not arise
            % too often...)
            fff=figure(103);
            set(fff,'Units','normalized','Position',[0.05 0.2 0.8 0.7]),
            plot(coord_fid_proj2D(:,1),coord_fid_proj2D(:,2),'ok','MarkerSize',50,'MarkerFaceColor','k');
            hold on;
            for iOptode=1:nOptodes
                    hold on, 
                    text(coord_fid_proj2D(iOptode,1),...
                        coord_fid_proj2D(iOptode,2),0,['F' int2str(iOptode)],...
                        'Color','r','FontWeight','bold','FontSize',14)
            end
            %set(gca,'Visible','off')
            xlim(1.5*[min(coord_helmet2D(:,1)) max(coord_helmet2D(:,1))])
            ylim(2.5*[min(coord_helmet2D(:,2)) max(coord_helmet2D(:,2))])
            idxSubjID = strfind(subjPath,'P2S');
            subjID = subjPath(idxSubjID:idxSubjID+4);
            title({[subjID ': Manual coregistration'];...
                'Go to Matlab Command Window and type in the numbers of the fiducials';...
                'in the image in the order that corresponds to S1-4, then D1-8 or 9, ';...
                'in square brackets.      Example: [13 9 6 2 11 10 7 1 12 8 5 4 3]'})
            optOrder = [];
            while (length(optOrder)<12) || (length(optOrder)>13) ||...
                    (sum(sort(optOrder,'ascend')-[1:nOptodes])~=0)
                optOrder = input(['Type in the numbers of the fiducials\n' ...
                    'in the image in the order that corresponds to S1-4,\n then D1-8 or 9, ' ...
                    'in square brackets. \n     Example: [13 9 6 2 11 10 7 1 12 8 5 4 3] :\n' ...
                    ' (or enter "0" to abort) \n\n']);
                if length(optOrder)==1 && optOrder==0
                    optOrder = 1:nOptodes;
                    break;
                end
            end
        end
        % Now let's go back to 3D coordinates, rearranging the order of the
        % markers to match the order that was to found to match the setup
        % helmet
        I = eye(nOptodes);
        coord_fid_opt = I(optOrder,:) * coord_fid;

        %close(gcf);    
        
    % Normal semi-automatic coregistration    
    else


        % Plot 2 geometries, then save as images and load them
        fff=figure;
        set(fff,'Units','normalized','Position',[0.05 0.2 0.8 0.7]),
        hfid=plot(coord_fid_proj2D(:,1),coord_fid_proj2D(:,2),'ok','MarkerSize',70,'MarkerFaceColor','k');
        set(gca,'Visible','off')
        xlim(1.5*[min(coord_helmet2D(:,1)) max(coord_helmet2D(:,1))])
        ylim(2.5*[min(coord_helmet2D(:,2)) max(coord_helmet2D(:,2))])
        saveas(hfid,fullfile(subjPath,'im_fid.png'),'png');
        close(gcf);
        figure('Units','normalized','Position',[0.05 0.2 0.8 0.7]),
        hhelmet=plot(coord_helmet2D(:,1),coord_helmet2D(:,2),'ob','MarkerSize',70,'MarkerFaceColor','b');
        xlim(1.5*[min(coord_helmet2D(:,1)) max(coord_helmet2D(:,1))])
        ylim(2.5*[min(coord_helmet2D(:,2)) max(coord_helmet2D(:,2))])
        set(gca,'Visible','off')
        saveas(hhelmet,fullfile(subjPath,'im_helmet.png'),'png');
        close(gcf);

        im_fid = double(imread(fullfile(subjPath,'im_fid.png'),'png'));
        im_helmet = double(imread(fullfile(subjPath,'im_helmet.png'),'png'));
        % No need for RGB info, just want a mask where background=0 and optodes=1
        im_fid = squeeze(im_fid(:,:,1));
        im_fid2 = im_fid;
        im_fid2(im_fid<max(im_fid(:)/2)) = 1; % optodes
        im_fid2(im_fid>=max(im_fid(:)/2)) = 0; % background
        im_helmet = squeeze(im_helmet(:,:,1));
        im_helmet2 = im_helmet;
        im_helmet2(im_helmet<max(im_helmet(:)/2)) = 1; % optodes
        im_helmet2(im_helmet>=max(im_helmet(:)/2)) = 0; % background
        % Crop borders (where was the figure background around the axes)
        lim1 = round(0.1*size(im_helmet,1));
        lim2 = round(0.1*size(im_helmet,2));
        im_helmet = im_helmet2(lim1:end-lim1,lim2:end-lim2);
        im_fid = im_fid2(lim1:end-lim1,lim2:end-lim2);
        clear im_fid2 im_helmet2

        % We now have 2 images of 2D probes and want to match their orientation
        % by testing all possible rotations and choosing the one that yields
        % the best match between both images

        angles = 1:2:360;
        correlation = [];
        for iTheta=1:length(angles)
            rotation = imrotate(im_fid,angles(iTheta),'bilinear','crop');
            correlation(iTheta)  =  sum(rotation(:).*im_helmet(:));
        end
        [cmax imax] = max(correlation);
        theta = angles(imax);
        %figure,imagesc(im_helmet)
        %figure,imagesc(imrotate(im_fid,theta,'bilinear','crop'))

        % The orientation of the coordinates will be either the angle found or
        % 180+it (it is almost symmetrical along the main axis). 
        theta = theta*pi/180;
        mat_rot1 = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        mat_rot2 = [cos(theta+pi) -sin(theta+pi); sin(theta+pi) cos(theta+pi)];
        coord_fid_proj2D_rot1 = [mat_rot1 * coord_fid_proj2D']';
        coord_fid_proj2D_rot2 = [mat_rot2 * coord_fid_proj2D']';
        % Also test for reflection around y axis
        coord_fid_proj2D_rot3 = coord_fid_proj2D_rot1;
        coord_fid_proj2D_rot3(:,2) = -1*coord_fid_proj2D_rot3(:,2);
        coord_fid_proj2D_rot4 = coord_fid_proj2D_rot2;
        coord_fid_proj2D_rot4(:,2) = -1*coord_fid_proj2D_rot4(:,2);

    %     if manualMode
    %         names_helmet = ['S1';'S2';'S3';'S4';'D1';'D2';'D3';'D4';'D5';'D6';'D7';'D8';'D9']; 
    %         figure, plot(coord_fid_proj2D_rot3(:,1),coord_fid_proj2D_rot3(:,2),'or');
    %         hold on, plot(coord_helmet2D(:,1),coord_helmet2D(:,2),'ob')
    %         for iOptode=1:nOptodes
    %             hold on, 
    %             text(coord_fid_proj2D_rot3(iOptode,1),...
    %                 coord_fid_proj2D_rot3(iOptode,2),0,['F' int2str(iOptode)],...
    %                 'Color','k')
    %             text(coord_helmet2D(iOptode,1),...
    %                 coord_helmet2D(iOptode,2),0,[names_helmet(iOptode,:)],...
    %                 'Color','b')
    %         end
    %     end

        % At this point, the 2D coordinates of the "theoretical" (setup) and
        % "experimental" (vitamins) helmets should be in relatively good
        % correspondance, and what remains to be done is to match each vitamin with
        % the right optode (of which we know the "names" from the setup
        % information). 

        % Two cases must be examined, because of the near-symmetry of the probe, 
        % its 180 degree rotation is very similar. In the end the best
        % correspondance will reveal the right choice. Idem for its reflection
        % about the main axis.

        % Distance between each optode and the closest fiducial
        % - this is not so robust in case the probe is much deformed so that the
        % 2D-projected coordinates do not closely match the setup helmet (or if the
        % SD distances in it are very small, similar in scale to that of ). Eventually find a fix...

        % Case rotation theta
        optOrder1 = zeros(1,nOptodes);
        for iOptode = 1:nOptodes-mod(nOptodes,2)
            dd = sqrt(sum((coord_fid_proj2D_rot1 - ones(nOptodes,1)*coord_helmet2D(iOptode,:)).^2,2));
            [err iFid] = min(dd);
            iFid2 = [];
            while ~isempty(find(optOrder1==iFid)) % if already assigned
                iFid2 = [iFid2 iFid];
                iFidDiff = setdiff(1:nOptodes-mod(nOptodes,2),iFid2);
                ddDiff = dd(iFidDiff);
                [err idxFid] = min(ddDiff); 
                iFid = iFidDiff(idxFid);
            end
           optOrder1(iOptode) = iFid;
        end
        % Pick the last one for the closest SD pair
        if mod(nOptodes,2); optOrder1(nOptodes)=setdiff(1:nOptodes,optOrder1); end
        % Case rotation theta+pi
        optOrder2 = zeros(1,nOptodes);
        for iOptode = 1:nOptodes-mod(nOptodes,2)
            dd = sqrt(sum((coord_fid_proj2D_rot2 - ones(nOptodes,1)*coord_helmet2D(iOptode,:)).^2,2));
            [err iFid] = min(dd);
            iFid2 = [];
            while ~isempty(find(optOrder2==iFid)) % if already assigned
                iFid2 = [iFid2 iFid];
                iFidDiff = setdiff(1:nOptodes-mod(nOptodes,2),iFid2);
                ddDiff = dd(iFidDiff);
                [err idxFid] = min(ddDiff); 
                iFid = iFidDiff(idxFid);
            end
            optOrder2(iOptode) = iFid;
        end
        % Pick the last one for the closest SD pair
        if mod(nOptodes,2); optOrder2(nOptodes)=setdiff(1:nOptodes,optOrder2); end
        % Case rotation theta & reflection
        optOrder3 = zeros(1,nOptodes);
        for iOptode = 1:nOptodes-mod(nOptodes,2)
            dd = sqrt(sum((coord_fid_proj2D_rot3 - ones(nOptodes,1)*coord_helmet2D(iOptode,:)).^2,2));
            [err iFid] = min(dd);
            iFid2 = [];
            while ~isempty(find(optOrder3==iFid)) % if already assigned
                iFid2 = [iFid2 iFid];
                iFidDiff = setdiff(1:nOptodes-mod(nOptodes,2),iFid2);
                ddDiff = dd(iFidDiff);
                [err idxFid] = min(ddDiff); 
                iFid = iFidDiff(idxFid);
            end
            optOrder3(iOptode) = iFid;
        end
        % Pick the last one for the closest SD pair
        if mod(nOptodes,2); optOrder3(nOptodes)=setdiff(1:nOptodes,optOrder3); end
        % Case rotation theta+pi & reflection
        optOrder4 = zeros(1,nOptodes);
        for iOptode = 1:nOptodes-mod(nOptodes,2)
            dd = sqrt(sum((coord_fid_proj2D_rot4 - ones(nOptodes,1)*coord_helmet2D(iOptode,:)).^2,2));
            [err iFid] = min(dd);
            iFid2 = [];
            while ~isempty(find(optOrder4==iFid)) % if already assigned
                iFid2 = [iFid2 iFid];
                iFidDiff = setdiff(1:nOptodes-mod(nOptodes,2),iFid2);
                ddDiff = dd(iFidDiff);
                [err idxFid] = min(ddDiff); 
                iFid = iFidDiff(idxFid);
            end
            optOrder4(iOptode) = iFid;
        end
        % Pick the last one for the closest SD pair
        if mod(nOptodes,2); optOrder4(nOptodes)=setdiff(1:nOptodes,optOrder4); end

        % Find out which of the 2 orientations yields the best match with the
        % setup helmet. Check using 2D distance.
        I = eye(nOptodes);
        coord_fid_perm1 = I(optOrder1,:)*coord_fid_proj2D_rot1;
        coord_fid_perm2 = I(optOrder2,:)*coord_fid_proj2D_rot2;
        coord_fid_perm3 = I(optOrder3,:)*coord_fid_proj2D_rot3;
        coord_fid_perm4 = I(optOrder4,:)*coord_fid_proj2D_rot4;
        %... find the one that corresponds best to the setup helmet
        dist1 = sum((coord_fid_perm1(:)-coord_helmet2D(:)).^2);
        dist2 = sum((coord_fid_perm2(:)-coord_helmet2D(:)).^2);
        dist3 = sum((coord_fid_perm3(:)-coord_helmet2D(:)).^2);
        dist4 = sum((coord_fid_perm4(:)-coord_helmet2D(:)).^2);

        if manualMode
            % For displaying the 2 geometries in 2D
            names_helmet = ['S1';'S2';'S3';'S4';'D1';'D2';'D3';'D4';'D5';'D6';'D7';'D8';'D9']; 
            figure(104), plot(coord_fid_proj2D_rot1(:,1),coord_fid_proj2D_rot1(:,2),'or');
            hold on, plot(coord_helmet2D(:,1),coord_helmet2D(:,2),'ob')
            for iOptode=1:nOptodes
                hold on, 
                text(coord_fid_proj2D_rot1(iOptode,1),...
                    coord_fid_proj2D_rot1(iOptode,2),0,['F' int2str(iOptode)],...
                    'Color','k')
                text(coord_helmet2D(iOptode,1),...
                    coord_helmet2D(iOptode,2),0,[names_helmet(iOptode,:)],...
                    'Color','b')
            end
            figure(105), plot(coord_fid_proj2D_rot2(:,1),coord_fid_proj2D_rot2(:,2),'or');
            hold on, plot(coord_helmet2D(:,1),coord_helmet2D(:,2),'ob')
            for iOptode=1:nOptodes
                hold on, 
                text(coord_fid_proj2D_rot2(iOptode,1),...
                    coord_fid_proj2D_rot2(iOptode,2),0,['F' int2str(iOptode)],...
                    'Color','k')
                text(coord_helmet2D(iOptode,1),...
                    coord_helmet2D(iOptode,2),0,[names_helmet(iOptode,:)],...
                    'Color','b')
            end
            figure(106), plot(coord_fid_proj2D_rot3(:,1),coord_fid_proj2D_rot3(:,2),'or');
            hold on, plot(coord_helmet2D(:,1),coord_helmet2D(:,2),'ob')
            for iOptode=1:nOptodes
                hold on, 
                text(coord_fid_proj2D_rot3(iOptode,1),...
                    coord_fid_proj2D_rot3(iOptode,2),0,['F' int2str(iOptode)],...
                    'Color','k')
                text(coord_helmet2D(iOptode,1),...
                    coord_helmet2D(iOptode,2),0,[names_helmet(iOptode,:)],...
                    'Color','b')
            end
            figure(107), plot(coord_fid_proj2D_rot4(:,1),coord_fid_proj2D_rot4(:,2),'or');
            hold on, plot(coord_helmet2D(:,1),coord_helmet2D(:,2),'ob')
            for iOptode=1:nOptodes
                hold on, 
                text(coord_fid_proj2D_rot4(iOptode,1),...
                    coord_fid_proj2D_rot4(iOptode,2),0,['F' int2str(iOptode)],...
                    'Color','k')
                text(coord_helmet2D(iOptode,1),...
                    coord_helmet2D(iOptode,2),0,[names_helmet(iOptode,:)],...
                    'Color','b')
            end
        end

        optOrders{1}=optOrder1;
        optOrders{2}=optOrder2;
        optOrders{3}=optOrder3;
        optOrders{4}=optOrder4;
        [dummy,optIdx] = min([dist1 dist2 dist3 dist4]);
        optOrder = optOrders{optIdx};
    %     if dist1<dist2
    %         optOrder = optOrder1;
    %     else
    %         optOrder = optOrder2;
    %     end
        % Now let's go back to 3D coordinates, rearranging the order of the
        % markers to match the order that was to found to match the setup
        % helmet
        I = eye(nOptodes);
        coord_fid_opt = I(optOrder,:) * coord_fid;
        
        
        % Last "security" fixes
        % This should be improved, because it is completely specific to one
        % geometry (2 rows of detectors surrounding on row of sources 
        % and in some cases + 1 detector (#13))...

        % 1) For 12 optodes
        % Finally, since the probe can be completely symmetrical around its main
        % axis (relative to a rotation around it), we must use this little fix
        % to ensure that we match the geometry correctly.
        if mod(nOptodes,2)==0
            %dets_height = NIRS.Cf.H.D.r.o.mm.p(:,2); % nOptodes x 3, mm
            % Convention: the first row of detectors (1:ndets/2) is on top, the
            % second one is below
            % So if 1st row seems lower...
            if sum(coord_fid_opt(nSrc+1:nSrc+floor(nDet/2),3)) < ...
                    sum(coord_fid_opt(nSrc+floor(nDet/2)+1:nSrc+2*floor(nDet/2),3))
                % then swap the 2 rows
                row1 = coord_fid_opt(nSrc+1:nSrc+floor(nDet/2),:);
                coord_fid_opt(nSrc+1:nSrc+floor(nDet/2),:) = ...
                    coord_fid_opt(nSrc+floor(nDet/2)+1:nSrc+2*floor(nDet/2),:);
                coord_fid_opt(nSrc+floor(nDet/2)+1:nSrc+2*floor(nDet/2),:) = row1;
            end

        % For 13 optodes : optode #13 must be closest to 3,4,11 in 
        % optiomized gemotry. If not, go back and choose another one of the 4
        % possible rotations...
        else
            % If the three optodes #3,4,11 are not the
            % closest to #13, we will try another configuration
            check = 0;
            nOrder = 0;
            temp = coord_fid_opt;
            while check < 3
                nOrder = nOrder + 1;
                if nOrder > 4
                    % Return to original order... but warning!
                    coord_fid_opt = temp;
                    idxSubjID = strfind(subjPath,'P2S');
                    subjID = subjPath(idxSubjID:idxSubjID+4);
                    disp(['Warning: potential error in coregistration for subject ' subjID]);
                    break; % no infinite loops!!
                end              
                coord_fid_opt = I(optOrders{nOrder},:) * coord_fid;
                dist_to_last_opt = sum(abs(coord_fid_opt -...
                    ones(size(coord_fid_opt),1)*coord_fid_opt(end,:)),2); % nOptodes x 1
                [nothing,order_close_to_13] = sort(dist_to_last_opt,'ascend');
                three_closest = order_close_to_13(2:4); % closest one is #13 itself
                check = length(find(three_closest==3))+length(find(three_closest==4))+...
                    length(find(three_closest==11));                      
            end 

        end
        
    end
    
    
    
    

    if manualMode
        % For displaying the optimized geometry in 3D
        figure(108), plot3(coord_fid_opt(1:nSrc,1),coord_fid_opt(1:nSrc,2),coord_fid_opt(1:nSrc,3),'or')
        hold on, plot3(coord_fid_opt(nSrc+1:nSrc+nDet,1),...
            coord_fid_opt(nSrc+1:nSrc+nDet,2),coord_fid_opt(nSrc+1:nSrc+nDet,3),'xb')
        for iF=1:nOptodes
            text(coord_fid_opt(iF,1),coord_fid_opt(iF,2),...
                coord_fid_opt(iF,3),['F' int2str(iF)])
        end
    end


    % Save the results of all these computations %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Create nifti object for anatomical without markers
    [dir,name,ext] = fileparts(Vanat.fname);    
    Vanat_WOspots = Vanat; 
    Vanat_WOspots.fname = fullfile(dir,[output_prefix_woVit,name,ext]);
    
    % Update NIRS matrix with new info:
    NIRS.Dt.ana.T1 = Vanat_WOspots.fname; % New path to anatomical
    % Positions of all points in mm : [pos_mm(x;y;z); 1] = V.mat * [pos_vox(x;y;z); 1]
    Pp_rmv = coord_fid_opt(:,[2 1 3])'; % column vectors of coordinates
    Pp_rmm = Vanat.mat * [Pp_rmv; ones(1,size(Pp_rmv,2))];
    NIRS.Cf.H.P.r.m.mm.p = Pp_rmm(1:3,:);
    NIRS.Cf.H.P.r.m.v.p = Pp_rmv;
    
    % Save outputs: anatomical without markers and NIRS matrix updated
    % with the positions of the markers
    NIRS.flags.detVit_OK = 1;
    save(NIRSmat{iSubj},'NIRS');
    spm_write_vol(Vanat_WOspots, Yanat_WOspots);
    
    % Avoid multiple windows when working in a loop
    for ii=1:8
        try close(100+ii); end
    end
    
    
    end  % end if needFixForMissingVit
    end

end



% Updated NIRS matrix is made available as a dependency
out.NIRSmat = NIRSmat;

end



function varargout = optodesCoord_click( volume,subjectID,initSrcPos,initDetPos )
% CLICKOPTODESPOS  Allow the user to manually select the positions of
%                 optodes relative to the volume previous to Monte-Carlo
%                 simulations. This GUI is called by the main GUI
%                 "prepMCsim". 
%
%       Comments displayed at the command line in response 
%       to the help command. 
% (Leave a blank line following the help.)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Initialization tasks (before components creation)  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Data initializations
    
    % Input arguments :
    % volume : 3D matrix of head MRI volume
    % subjectID : name used for window title
    % initSrc/DetPos : initial opotdes positions (those already entered
    % before, the newly entered ones will be added to the list)
    if ~exist('initSrcPos','var') || isempty(initSrcPos)
        initSrcPos = zeros(0,3);
    end
    if ~exist('initDetPos','var') || isempty(initDetPos)
        initDetPos = zeros(0,3);
    end

    % Output arguments : will contain user-defined optodes positions, or an
    % empty matrix if no position is entered or figure is closed
    outputPos = {};    
    
    % Optodes positions (to be updated to match with user-defined positions
    % (from clicks or text editing)
    optodes.srcPos = initSrcPos;
    optodes.detPos = initDetPos;
    
    % Volume min and max values (for color mapping for display)
    clims = [min(volume(:)) max(volume(:))];


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Construct the components  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%  Create and then hide the GUI as it is being constructed.
   hfig = figure('Visible','off','Position',[100 200 1100 700],...
          'Color',[.925 .914 .847],...
          'DeleteFcn',{@hfig_DeleteFcn});

% --- Project title and instructions for the user
   htitle = uicontrol('Style','text','String',char(subjectID),...
       'Position',[0 670 250 30],'FontSize',14,...
       'Min',0,'Max',2,...
       'BackgroundColor',[0 0 0],'ForegroundColor',get(hfig,'Color'));
   hinstructions = uicontrol('Style','text',...
       'Position',[50 270 1000 20],'FontSize',12,...
       'Min',0,'Max',2,...
       'BackgroundColor',[.8 .8 .75],'ForegroundColor',[0 0 0]);
   set(hinstructions,'String',['Repérez les optodes',...
       ' en cliquant sur un des trois axes,',...
       ' puis à la position de la source (clic gauche)',...
       ' ou du détecteur (clic droit).']);

% --- Axes and sliders for volume display and graphical inputs
   haxes_1axial = axes('Parent',hfig,'NextPlot','replacechildren',...
                  'Units','pixels','Position',[50 370 300 250],...
                  'Tag','axes_1axial',...
                  'ButtonDownFcn',{@axes_ButtonDownFcn});
                  title({'Bas -> haut'; '  -  Axiale  -  '},...
                      'HorizontalAlignment','Right',...
                      'Units','normalized',...
                      'Position',[1.25 1.0]);
   haxes_2coronal = axes('Parent',hfig,'NextPlot','replacechildren',...
                    'Units','pixels','Position',[400 370 300 250],...
                    'Tag','axes_2coronal',...
                    'ButtonDownFcn',{@axes_ButtonDownFcn});
                    title({'Derrière -> Devant'; ' - Coronale - '},...
                     'HorizontalAlignment','Right',...
                      'Units','normalized',...
                      'Position',[1.25 1.0]);
   haxes_3sagittal = axes('Parent',hfig,'NextPlot','replacechildren',...
                     'Units','pixels','Position',[750 370 300 250],...
                     'Tag','axes_3sagittal',...
                     'ButtonDownFcn',{@axes_ButtonDownFcn});
                     title({'Gauche -> droite'; ' - Sagittale - '},...
                      'HorizontalAlignment','Right',...
                      'Units','normalized',...
                      'Position',[1.25 1.0]);
   set([haxes_1axial haxes_2coronal haxes_3sagittal],...
       'DataAspectRatio',[1 1 1],...
       'PlotBoxAspectRatio',[1 1 1]);
   
   hslider_1axial = uicontrol(hfig,'Style','slider',...
                    'Position',[50 320 300 20],...
                    'Value',1,...
                    'Callback',{@slider_1axial_Callback});
   hslider_2coronal = uicontrol(hfig,'Style','slider',...
                    'Position',[400 320 300 20],...
                    'Value',1,...
                    'Callback',{@slider_2coronal_Callback});
   hslider_3sagittal = uicontrol(hfig,'Style','slider',...
                    'Position',[750 320 300 20],...
                    'Value',1,...
                    'Callback',{@slider_3sagittal_Callback});  
                
   % Volume dimension for initializing sliders
   dim = size(volume);
          set(hslider_1axial,'Enable','on',...
              'Min',1,'Max',dim(3),...
              'SliderStep',[1/dim(3) 15/dim(3)]);
          set(hslider_2coronal,'Enable','on',...
              'Min',1,'Max',dim(2),...
              'SliderStep',[1/dim(2) 15/dim(2)]);
          set(hslider_3sagittal,'Enable','on',...
              'Min',1,'Max',dim(1),...
              'SliderStep',[1/dim(1) 15/dim(1)]);

% --- Edit texts for text inputs
   hpanel_edits = uipanel('Parent',hfig,'Title','Positions entrées',...
                  'BackgroundColor',[.8 .8 .75],...
                  'Units','pixels','Position',[50 25 500 225]);   
     hedit_srcPos = uicontrol(hpanel_edits,'Style','edit',...
                       'Units','normalized','Position',[.1 .05 .35 .9],...
                       'Min',0,'Max',10,'String',num2str(optodes.srcPos),...
                       'Tag','s',...
                       'Callback',{@edit_Pos_Callback});
     hedit_detPos = uicontrol(hpanel_edits,'Style','edit',...
                       'Units','normalized','Position',[.55 .05 .35 .9],...
                       'Min',0,'Max',10,'String',num2str(optodes.detPos),...
                       'Tag','d',...
                       'Callback',{@edit_Pos_Callback});


% --- Pushbuttons for user interactions
   hpushbutton_OK = uicontrol(hfig,'Style','pushbutton',...
                  'String','Terminé',...
                  'Units','pixels','Position',[800 110 150 130],...
                  'Callback',{@pushbutton_OK_Callback});
   hpushbutton_cancel = uicontrol(hfig,'Style','pushbutton',...
                  'String','Annuler',...
                  'Units','pixels','Position',[800 25 150 50],...
                  'Callback',{@pushbutton_cancel_Callback});

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Initialization tasks (after components creation)  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% --- GUI initializations
   set([hfig htitle haxes_1axial...
       haxes_2coronal haxes_3sagittal hslider_1axial...
       hslider_2coronal hslider_3sagittal],...
   'Units','normalized'); % So components resize automatically

   % General figure appearance 
   set(hfig,'Name','Repérage des optodes',... % Window title
       'NumberTitle','off'); % Do not display figure number
   % Move the GUI to the center of the screen, then up a little
   movegui(hfig,'center');
   set(hfig,'Position',get(hfig,'Position')+[0 .1 0 0]);
   
   % Initialize each plot and invert vertical axes
   % (see explanation below in "draw" function
    updateDisplay(1:3);
    initializeAxesTicks;
    
   
   % Make the GUI visible.
   set(hfig,'Visible','on');
   
   % Block execution until user has interacted with or closed the figure
    uiwait(hfig);

    % Return the coordinates selected by the user if it is requested
    outputPos{1} = optodes.srcPos;
    outputPos{2} = optodes.detPos;
    if nargout>0
        [varargout{1:nargout}] = outputPos{:};
    end

   
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Callbacks for CLICKOPTODESPOS  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % These callbacks automatically
   %  have access to component handles and initialized data 
   %  because they are nested at a lower level.
 

 
% --- GRAPHICAL INPUTS

   % --- Extracting optode coordinates upon user click on a location
    function axes_ButtonDownFcn(hObject,eventdata)
        % When the user double-clicks, the first click will call this
        % ButtonDown function and the second one will respond to the
        % "ginput" command:
        [horizClicked vertClicked SorD] = ginput(1);
        %get(hObject,'ButtonDownFcn')
        % The newly selected optode will be added to the list
        updateClick(hObject,round(horizClicked),round(vertClicked),SorD);
%         % Re-define the ButtonDown function which has been reset (why?)
%         set(hObject,'ButtonDownFcn',{@axes_ButtonDownFcn})
    end

    % --- Slider callbacks
     % --- Navigate thourgh 3D current volume with sliders
   function slider_1axial_Callback(hObject,eventdata)
       % Update display on 1st plot
       updateDisplay(1);
   end

   function slider_2coronal_Callback(hObject,eventdata)
       % Update display on second plot
       updateDisplay(2);
   end

   function slider_3sagittal_Callback(hObject,eventdata)
       % Update display on 3rd plot
       updateDisplay(3);
   end


% --- EDIT TEXT INPUTS
      % Edit texts are updated whenever a graphical input is given
      % Graphics are updated to match new entered position

    function edit_Pos_Callback(hObject,eventdata)
        SorD = get(hObject,'Tag'); % source or detector edit text
        updateEdit(hObject,SorD);
    end


% --- FINISHING

   % --- Terminate selections
    function pushbutton_OK_Callback(hObject,eventdata)
        % Values in optodes.srcPos and optodes.detPos will be output when
        % execution resumes
       uiresume;
        % End
        delete(hfig);
    end

     function pushbutton_cancel_Callback(hObject,eventdata)
        % Cancel : return to initial values (possibly empty)
        optodes.srcPos = initSrcPos;
        optodes.detPos = initDetPos;
        uiresume;
        % End
        delete(hfig);
      end

    % --- If figure is closed prematurely
    function hfig_DeleteFcn(hObject,eventdata)
        % Return the positions anyway
        uiresume;
    end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Utility functions for PREPMCSIM2  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % --- Initialize vertical axes ticks (inverted because of imrotate)
   % (see exmplanation in "draw" function below)
   % Must have initialized display with updateDisplay prior to this!
    function initializeAxesTicks
        axesHandles = [haxes_1axial; haxes_2coronal; haxes_3sagittal];
        for view = 1:3 % for each axes
            hax = axesHandles(view);
            if ~isempty(get(hax,'Children')) % make sure display was initialized
                ticksYdef = get(hax,'YTick'); % get current ticks (ex: 0-256; i=5,255)
                %ticksYnew = size(imrotate(image2D,90),1) - ticksYdef;
                ticksYnew = max(get(hax,'ylim')) - ticksYdef; % calculate new ticks : y -> -y + max(y) (invert) (i=251,1)
                set(hax,'YTick',sort(ticksYnew,'ascend')) % set new ticks (at i=251,1)
                set(hax,'YTickLabel',num2str(sort(ticksYdef','descend'))) % but label them (label i=5,255)
            end
        end
    end

   
   % --- Update display to match current matrix choice and display options
    function updateDisplay(viewNumbers)
        % viewNumbers indicate which plot(s) to update : numbers 1, 2 or 3
        
        % Display each view (axial, coronal, sagittal) in the corresponding axes
        sliderHandles = [hslider_1axial; hslider_2coronal; hslider_3sagittal];
        axesHandles = [haxes_1axial; haxes_2coronal; haxes_3sagittal];
        if ~isempty(volume)
            for view = viewNumbers
                % Slice number is indicated by slider
                sliceNo = ceil(get(sliderHandles(view),'Value'));
                switch view
                    case 1
                        theSlice = squeeze(volume(:,:,sliceNo));
                        coordHoriz = 1;
                        coordVert = 2;
                    case 2
                        theSlice = squeeze(volume(:,sliceNo,:));
                        coordHoriz = 1;
                        coordVert = 3;
                    case 3
                        theSlice = squeeze(volume(sliceNo,:,:));
                        coordHoriz = 2;
                        coordVert = 3;
                    otherwise
                        % Plots are not updated
                end
                
                % Clear axes of previous children (image, rectangle...)
                cla(axesHandles(view));
                
                % Draw the image of the slice on the plot, indicate slice
                % number on top and label colorbar on the left
                draw(theSlice,axesHandles(view),clims,...
                    sliceNo);

                if any(optodes.srcPos) % inialized to [0 0 0]
                % Display sources and detectors which are in this plane
                    for nsrc = 1:size(optodes.srcPos,1); % For each source of known position
                        if optodes.srcPos(nsrc,4-view) == sliceNo;
                        % 4 - view : axial/z , coronal/y, sagittal/x
                             axes(axesHandles(view))
                             hold on, plot(optodes.srcPos(nsrc,coordHoriz),...
                                 size(imrotate(theSlice,90),1)-optodes.srcPos(nsrc,coordVert),... % vertical axis is inverted
                                 '*y','MarkerSize',10);
                        end
                    end
                end
                if any(optodes.detPos)
                    for ndet = 1:size(optodes.detPos,1); % For each source of known position
                        if optodes.detPos(ndet,4-view) == sliceNo;
                        % 4-view : axial/z , coronal/y, sagittal/x
                             axes(axesHandles(view))
                             hold on, plot(optodes.detPos(ndet,coordHoriz),...
                                 size(imrotate(theSlice,90),1)-optodes.detPos(ndet,coordVert),... % vertical axis is inverted
                                 'og','MarkerSize',10);
                        end
                    end
                end
                
            end
            
        else % clear (but do not reset) axes
            axes(haxes_1axial); cla;
            axes(haxes_2coronal); cla;
            axes(haxes_3sagittal); cla;
        end
        
    end
   

    % --- Draw an image on specified axes
    function draw(image2D,h_axes,clims,sliceNumber,ticks)

    % inputs :
    % plan : 2D mage to display (could be a slice of a 3D matrix)
    % h_axes : handles to (existing) axes where the image will be drawn
    % optional inputs :
    % clims : [min max] limits for color mapping (see imagesc)
    % sliceNumber : will be displayed on top of axes (ID the slice within
    % the 3D volume)
    % ticks : labels for the colorbar in a cell array of strings
    % ex: {'Air'; 'Mgr'; 'Mbl'; 'CSF'; 'Skull'; 'Scalp'}

    % Remove singleton dimensions
    if size(size(image2D))>2
        image2D = squeeze(image2D);
    end
    
    % Default values for arguments
        switch nargin
            case 2
                clims = [min(image2D(:)) max(image2D(:))];
                sliceNumber = [];
                ticks = [];
            case 3
                sliceNumber = 0;
                ticks = [];
            case 4
                ticks = [];
            case 5
            otherwise
                return
        end

       
        axes(h_axes)
        imagesc(imrotate(image2D,90),clims);
        axis ij; axis tight; colormap gray;
        % Rotate because imagesc plots the first coordinate vertical
        % and the second horizontal. However this rotation causes the
        % vertical axis to be inverted relative to the convention, so reset
        % the graduations so that they reflect the real axis direction.     
       
     
        % Colorbar labelled or not
        if ~isempty(ticks)
            colorbar('Location','EastOutside','Ytick',0:length(ticks)-1,'YTickLabel',ticks)
            % maximum 5 tissues + air
        else
            colorbar('Location','EastOutside')
        end

        % Display slice number within the 3D matrix
        if sliceNumber~=0
            txt = text(0.01,1.04,'','Units','normalized');
            set(txt,'String',strcat('Tranche #',int2str(sliceNumber)),...
                'BackgroundColor',[.8 .8 .75],'Color','k','FontSize',10,'FontWeight','Normal');
        end

        
    end
        

    % --- Update the list and display of optodes positions to include the
    % newly clicked optode (source or detector)
    function updateClick(haxes,horizPos,vertPos,SorD)
        % Optode type : source (1) or detector (3) (if 2, do nothing)
        
        % Coordinates of clicked point
        switch get(haxes,'Tag')
            case 'axes_1axial'
                % Position sélectionnée :
                x = round(horizPos);
                y = round(max(get(haxes,'ylim')) - vertPos); % inverted graduation on vertical axis
                z = ceil(get(hslider_1axial,'Value'));

            case 'axes_2coronal'
                % Position sélectionnée :
                x = round(horizPos);
                z = round(max(get(haxes,'ylim')) - vertPos); % inverted graduation on vertical axis
                y = ceil(get(hslider_2coronal,'Value'));
               
            case 'axes_3sagittal'
                % Position sélectionnée :
                y = round(horizPos);
                z = round(max(get(haxes,'ylim')) - vertPos); % inverted graduation on vertical axis
                x = ceil(get(hslider_3sagittal,'Value'));
        end
        
        % Protection : abort if user clicked outside the axes
        dim = size(volume);
        if ( (1<=x && x<=dim(1)) && (1<=y && y<=dim(2)) && (1<=z && z<=dim(3)) )
            % List of positions
            if SorD == 1 % left click : source
                optodes.srcPos = [optodes.srcPos; x y z];
            elseif SorD == 3 % right click : detector
                optodes.detPos = [optodes.detPos; x y z];
            end % middle click : do nothing

            % Display optode on graphics
            updateDisplay(1:3);

            % Update list of positions in edit text
            set(hedit_srcPos,'String',num2str(optodes.srcPos));
            set(hedit_detPos,'String',num2str(optodes.detPos));
            
        end
        
    end




  % --- Update the list and display of optodes positions to include the
    % newly enetered optode (source or detector)
    function updateEdit(hObject,SorD)
        % hObject : handles to the edit text object that called the
        % function
        % SorD : 's' (source) or 'd' (detector)
        
        % Get the coordinates entered
        coords = str2num(get(hObject,'String'));
        if ~isempty(coords)
            allx = round(coords(:,1));
            ally = round(coords(:,2));
            allz = round(coords(:,3));
            dim = size(volume); % Make sure the new position is inside the volume
            % (and while we are at it, this will also check that the input is
            % really a number (char, Inf or NaN will not respect this condition)
            if isreal(sum(allx(:))+sum(ally(:))+sum(allz)) && ...
               isfinite(sum(allx(:))+sum(ally(:))+sum(allz)) && ...
                (1 <= min(allx)) && (max(allx) <= dim(1)) && ...
                (1 <= min(ally)) && (max(ally) <= dim(2)) && ...
                (1 <= min(allz)) && (max(allz) <= dim(3))
                    % Update list of positions
                    if strcmp(SorD,'s')
                        optodes.srcPos = [allx ally allz];
                    elseif strcmp(SorD,'d')
                        optodes.detPos = [allx ally allz];
                    end
                    % Display optodes on graphics
                    updateDisplay(1:3);
            else  % Do not consider the last change
                if strcmp(SorD,'s')
                    set(hObject,'String',num2str(optodes.srcPos));
                elseif strcmp(SorD,'d')
                    set(hObject,'String',num2str(optodes.detPos));
                end
                    
            end

        elseif ~isempty(get(hObject,'String')) % Unless user did delete all positions,
            % do not consider the last change
            if strcmp(SorD,'s')
                set(hObject,'String',num2str(optodes.srcPos));
            elseif strcmp(SorD,'d')
                set(hObject,'String',num2str(optodes.detPos));
            end
        end
            
    end

    
    


 
end




