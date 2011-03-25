function out = nirs_run_detectVitamins(job)
% This module is used for semi-automatically detecting fiducial markers 
% (e.g. marking the positions of the optodes) on an MRI anatomical image.
% The positions are detected and saved in the NIRS.mat matrix, and the 
% the markers are erased from the anatomical image (a copy of the original
% image is saved).
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
%
% Michèle Desjardins
% 2011-03

% INPUTS %
%%%%%%%%%%

% NIRS matrices for all subjects
NIRSmat = job.NIRSmat;
%load(NIRSmat{:});

% Path to anatomical image for all subjects
image_in = job.image_in;
%V = spm_vol(anatT1{:});
%Y = spm_read_vols(V);

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
    
    % Read anatomical image and NIRS matrix
    try    
        load(NIRSmat{iSubj});
        Vanat = spm_vol(image_in{iSubj});
        Yanat = spm_read_vols(Vanat);
    catch
        disp(['Could not read NIRS.mat or anatomical image for ' int2str(iSubj) 'th subject.' ]);
    end
    
    % Define structuring elements
    se =  strel('ball',5,5);
    se_larger = strel('ball',7,7);
    
    % Fill head image (paint-bucket)
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
    for iFiducial=1:size(stats,1)
        coord_fid(iFiducial,:) = stats(iFiducial).Centroid(:);
        size_fid(iFiducial,:) = stats(iFiducial).Area(:);
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
    
    coord_fid = coord_fid(toKeep,:);
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
    u2 = u2./norm(u2);
    % (It can be shown that these 2 vectors belong to the plane and that they
    % are not only linearly independent, but also orthogonal)
    % The projection of the coordinates on this plane is given by the (vector)
    % sum of the projection on each of the vectors of the orthogonal basis,
    % i.e. Proj = (coord*basisVector/(basisVector'*basisVector)) .* basis vector
    for iOptode = 1:nOptodes
        coord_fid_proj2D(iOptode,:) =  [ (coord_fid(iOptode,:)*u1)./(u1'*u1) ...
                                       (coord_fid(iOptode,:)*u2)./(u2'*u2) ];
    end

    % Center on 0 in u1-u2 plane
    ctr_mass_fid = [mean(coord_fid_proj2D(:,1)); mean(coord_fid_proj2D(:,2))];
    coord_fid_proj2D = coord_fid_proj2D - (ctr_mass_fid*ones(1,size(coord_fid_proj2D,1)))'; 
    
    % Center helmet on 0 and use 2D coordinates
    ctr_mass_helmet = [mean(coord_helmet(:,1)); mean(coord_helmet(:,2)); mean(coord_helmet(:,3))];
    coord_helmet = coord_helmet - (ctr_mass_helmet*ones(1,size(coord_helmet,1)))'; 
    coord_helmet2D = coord_helmet(:,1:2);
    
    % We must now match 2 sets of 2D coordinates, which will require to
    % first find the right orientation, using some manipulations of the
    % images of the coordinates
    
    % Plot 2 geometries, then save as images and load them
    figure('Units','normalized','Position',[0.05 0.2 0.8 0.7]),
    hfid=plot(coord_fid_proj2D(:,1),coord_fid_proj2D(:,2),'ok','MarkerSize',70,'MarkerFaceColor','k');
    set(gca,'Visible','off')
    xlim(1.5*[min(coord_helmet2D(:,1)) max(coord_helmet2D(:,1))])
    ylim(2.5*[min(coord_helmet2D(:,2)) max(coord_helmet2D(:,2))])
    saveas(hfid,'im_fid.png','png');
    figure('Units','normalized','Position',[0.05 0.2 0.8 0.7]),
    hhelmet=plot(coord_helmet2D(:,1),coord_helmet2D(:,2),'ob','MarkerSize',70,'MarkerFaceColor','b');
    xlim(1.5*[min(coord_helmet2D(:,1)) max(coord_helmet2D(:,1))])
    ylim(2.5*[min(coord_helmet2D(:,2)) max(coord_helmet2D(:,2))])
    set(gca,'Visible','off')
    saveas(hhelmet,'im_helmet.png','png');

    im_fid = double(imread('im_fid.png','png'));
    im_helmet = double(imread('im_helmet.png','png'));
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

%     figure, plot(coord_fid_projP_rot2(:,1),coord_fid_projP_rot2(:,2),'or');
%     hold on, plot(coord_helmet2D(:,1),coord_helmet2D(:,2),'ob')
%     for iOptode=1:nOptodes
%         hold on, 
%         text(coord_fid2_proj_ctr_rot2(iOptode,1),...
%             coord_fid2_proj_ctr_rot2(iOptode,2),0,['F' int2str(iOptode)],...
%             'Color','k')
%         text(coord_helmet2D(iOptode,1),...
%             coord_helmet2D(iOptode,2),0,[names_helmet(iOptode,:)],...
%             'Color','b')
%     end

    % At this point, the 2D coordinates of the "theoretical" (setup) and
    % "experimental" (vitamins) helmets should be in relatively good
    % correspondance, and what remains to be done is to match each vitamin with
    % the right optode (of which we know the "names" from the setup
    % information). 

    % Two cases must be examined, because of the near-symmetry of the probe, 
    % its 180 degree rotation is very similar. In the end the best
    % correspondance will reveal the right choice.

    % Distance between each optode and the closest fiducial
    % - this is not so robust in case the probe is much deformed so that the
    % 2D-projected coordinates do not closely match the setup helmet (or if the
    % SD distances in it are very small, similar in scale to that of ). Eventually find a fix...
    
    % Case rotation theta
    optOrder1 = zeros(1,nOptodes);
    for iOptode = 1:nOptodes
        dd = sqrt(sum((coord_fid_proj2D_rot1 - ones(nOptodes,1)*coord_helmet2D(iOptode,:)).^2,2));
        [err iFid] = min(dd);
        optOrder1(iOptode) = iFid;
    end
    % Case rotation theta+pi
    optOrder2 = zeros(1,nOptodes);
    for iOptode = 1:nOptodes
        dd = sqrt(sum((coord_fid_proj2D_rot2 - ones(nOptodes,1)*coord_helmet2D(iOptode,:)).^2,2));
        [err iFid] = min(dd);
        optOrder2(iOptode) = iFid;
    end
    % Find out which of the 2 orientations yields the best match with the
    % setup helmet. Check using 2D distance.
    I = eye(nOptodes);
    coord_fid_perm1 = I(optOrder1,:)*coord_fid_proj2D_rot1;
    coord_fid_perm2 = I(optOrder2,:)*coord_fid_proj2D_rot2;
    %... find the one that corresponds best to the setup helmet
    dist1 = sum((coord_fid_perm1(:)-coord_helmet2D(:)).^2);
    dist2 = sum((coord_fid_perm2(:)-coord_helmet2D(:)).^2);

    if dist1<dist2
        optOrder = optOrder1;
    else
        optOrder = optOrder2;
    end

    % Now let's go back to 3D coordinates, rearranging the order of the
    % markers to match the order that was to found to match the setup
    % helmet
    I = eye(nOptodes);
    coord_fid_opt = I(optOrder,:) * coord_fid;

%     figure, plot3(coord_fid_opt(1:nSrc,1),coord_fid_opt(1:nSrc,2),coord_fid_opt(1:nSrc,3),'or')
%     hold on, plot3(coord_fid_opt(nSrc+1:nSrc+nDet,1),...
%         coord_fid_opt(nSrc+1:nSrc+nDet,2),coord_fid_opt(nSrc+1:nSrc+nDet,3),'xb')
%     for iF=1:nOptodes2
%         text(coord_fid_opt(iF,1),coord_fid_opt(iF,2),...
%             coord_fid_opt(iF,3),['F' int2str(iF)])
%     end


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
    save(NIRSmat{iSubj},'NIRS');
    spm_write_vol(Vanat_WOspots, Yanat_WOspots);
    
end

% Updated NIRS matrix is made available as a dependency
out.NIRSmat = NIRSmat;

return


