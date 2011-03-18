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
    
%     addpath(('Z:\Programmes\spm8_300710_r4010'));
%     addpath(genpath('G:\MesProgrammes\SPM\pour_nirs10'))
%     addpath(genpath('G:\MesProgrammes\tMCimg_scripts\Scripts\PourPreparerSimuMC\New_SPMnewSegment'));
%     cd('V:\PreProc_Phase2\P2S1\test_nirs10');
%     V=spm_vol('s201009150900-0004-00001-000176-01.nii');
%     Y=spm_read_vols(V);

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
    %optodesCoord_click(test3,'bob',round(coord([2 1 3])),[]);
    
    % Read helmet information
    
    % Either keep helemt info (src & det positions) in NIRS matrix or use
    % template BrainSight-coregistered info
    %nn=load('V:\PreProc_Phase2\P2S1\IOD\P2S1_scan3.nirs','-mat')
    coord_helmet = [nn.SD.SrcPos; nn.SD.DetPos];
    names_helmet = ['S1';'S2';'S3';'S4';'D1';'D2';'D3';'D4';'D5';'D6';'D7';'D8';];
    
   
    
    % We now have the coordinates of the centroïds of the fiducials. In
    % order to identify each source and detector and their names, as well
    % as correct for potential errors (fake markers due to noise mistaken
    % for a fiducial), we will compare with the information given in the
    % NIRS matrix regarding the helmet used for the acquisition.
   
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
    
    %load('spots.mat');
    %Vanat=spm_vol('s201009150900-0004-00001-000176-01.nii');
    Pp_rmv = coord_fid'; % column vectors of coordinates
    Pp_rmm = Vanat.mat * [Pp_rmv; ones(1,size(Pp_rmv,2))];
    %Pp_rmm = Pp_rmm(1:3,:);
    
    ctr_mass = [mean(Pp_rmm(1,:)); mean(Pp_rmm(2,:)); mean(Pp_rmm(3,:))];
    Pp_rmm_ctr = Pp_rmm - ctr_mass*ones(1,size(Pp_rmm,2));  
    
    figure;
    for idet=1:size(Pp_rmm,2)
        xx(idet)=Pp_rmm(1,idet);
        yy(idet)=Pp_rmm(2,idet);
        zz(idet)=Pp_rmm(3,idet);
    end
    plot3(xx,yy,zz,'or','MarkerSize',14,'MarkerFaceColor','r')
    figure, plot(yy,zz,'or','MarkerSize',14,'MarkerFaceColor','r')
%     pca = princomp(Pp_rmm(2:3,:)');
%     for idet=1:size(pca,1)
%         xx1(idet)=pca(idet,1);
%         yy1(idet)=pca(idet,2);
%         %zz1(idet)=pca(idet,3);
%     end
%     %hold on, plot3(xx1,yy1,zz1,'xb','MarkerSize',14,'MarkerFaceColor','b')
%     figure, plot(yy,zz,'or','MarkerSize',14,'MarkerFaceColor','r')
%     hold on, plot(xx1,yy1,'xb','MarkerSize',14,'MarkerFaceColor','b')

    % Projeté sur un des plans (je pense coronal - le plan yz (enlever coord x)),
    % on voit 2 directions principales du "quadrillé"... peut-etre des "composantes principales"
    % ou qqch du genre pourrait nous donner la bonne direction où regarder
    % les coord extrêmes?
    

% Nice idea, maybe in 2050...
% ----------------------------------------------------------------------- % 
    % Minimize the Euclidian distance (L-2 norm) between the "theoretical"
    % helmet (from the setup info in the NIRS matrix) and the ensemble of
    % positions identified using the fiducials:
    
%    I = eye(nOptodes);
%    % Of all possible rearragements of the optodes...
%    thePermutations = perms(1:nOptodes);
%    for iPerm=1:factorial(nOptodes)
%        I = I(thePermutations(iPerm,:),:);
%        coord_fid_perm = I*coord_fid;
%        %... find the one that corresponds best to the setup helmet
%        dist(iPerm) = (sum(coord_fid_perm(:)-coord_helmet(:)).^2);
%    end
%    
%    [theMinDist,theIdx] = min(dist);
%    theOrder = thePermutations(theIdx,:);
%    coord_fid_opt = I(theOrder,:) * coord_fid;
%    names_fid_opt = inv(I(theOrder),:)*names_helmet;
%    
%    % Check
%    nSrc = nn.SD.nSrc;
%    nDet = nn.SD.nDet;
%    optodesCoord_click(Yanat,'Verify correspondance of the optodes',...
%        coord_fid_opt([2 1 3],1:nSrc),coord_fid_opt([2 1 3],nSrc+1:nSrc+nDet))
% ----------------------------------------------------------------------- % 
   
    
    % Create nifti object for anatomical without markers
    [dir,name,ext] = fileparts(Vanat.fname);    
    Vanat_WOspots = Vanat; 
    Vanat_WOspots.fname = fullfile(dir,[output_prefix_woVit,name,ext]);
    
     % Update NIRS matrix with new info:
    NIRS.Dt.ana.T1 = Vanat_WOspots.fname; % New path to anatomical
    % Positions of all points in mm : [pos_mm(x;y;z); 1] = V.mat * [pos_vox(x;y;z); 1]
    Pp_rmv = coord_fid'; % column vectors of coordinates
    Pp_rmm = Vanat.mat * [Pp_rmv; ones(1,size(Pp_rmv,2))];
    NIRS.Cf.H.P.r.m.mm.p = Pp_rmm(1:3,:);
    NIRS.Cf.H.P.r.m.v.p = Pp_rmv;
    
    % Save outputs: anatomical without markers and NIRS matrix updated
    % with the positions of the markers
    save(outNIRSmat{iSubj},'NIRS');
    spm_write_vol(Vanat_WOspots, Yanat_WOspots);
    
end

% clear NIRS
% load(NIRSmat{:});
% NIRS.Cf.H.P.r.m.mm.fp =Pfp_rmm;
% NIRS.Cf.H.P.r.m.vx.fp = Pfp_rmv(1:3,:);
% save(NIRSmat{:},'NIRS');

%out = 1;
% Updated NIRS matrix is made available as a dependency
out.NIRSmat = outNIRSmat;

return






% OLD TRIES
% 
% im2d=squeeze(Y(:,:,:));
%         t3 = 35;
%         t2 = 168;
%         % Tranche où on voit 4 sources!: Y(:,:,14)
%         figure(4), subplot(1,2,1), imagesc(squeeze(im2d(:,:,t3)));
%         subplot(1,2,2), imagesc(squeeze(im2d(:,t2,:)));
%         %se = strel('disk',9,0);    
%         se =  strel('ball',5,5);
%         se2 = strel('ball',7,7);
%         
%         for i=1:size(im2d,3)
%             filledim2d(:,:,i) = imfill(im2d(:,:,i),'holes');
%         end
%         figure(44), subplot(1,2,1), imagesc(squeeze(filledim2d(:,:,t3)));
%         subplot(1,2,2), imagesc(squeeze(filledim2d(:,t2,:)));
%         
%         erodedBW = imerode(filledim2d,se);
%         figure(5), subplot(1,2,1), imagesc(squeeze(erodedBW(:,:,t3)));
%         subplot(1,2,2), imagesc(squeeze(erodedBW(:,t2,:)));
%         
%         dilatedBW = imdilate(erodedBW,se);
%         figure(6), subplot(1,2,1), imagesc(squeeze(dilatedBW(:,:,t3)));
%         subplot(1,2,2), imagesc(squeeze(dilatedBW(:,t2,:)));
%         
%         closedBW = imclose(dilatedBW,se2);
%         figure(7), subplot(1,2,1), imagesc(squeeze(closedBW(:,:,t3)));
%         subplot(1,2,2), imagesc(squeeze(closedBW(:,t2,:)));
%         
% %         filledBW = imfill(closedBW,'holes');
% %         figure(8), subplot(1,2,1), imagesc(squeeze(filledBW(:,:,t3)));
% %         subplot(1,2,2), imagesc(squeeze(filledBW(:,t2,:)));
% %         
%         for i=1:size(closedBW,3)
%             filledBW(:,:,i) = imfill(closedBW(:,:,i),'holes');
%         end
%         figure(8), subplot(1,2,1), imagesc(squeeze(filledBW(:,:,t3)));
%         subplot(1,2,2), imagesc(squeeze(filledBW(:,t2,:)));
%         
%         %thresh = 0.05;
%         for i=1:size(filledBW,2)
%             filtBW(:,i,:) = medfilt2(squeeze(filledBW(:,i,:)));
%         end
%         for i=1:size(filtBW,1)
%             filtBW(i,:,:) = medfilt2(squeeze(filtBW(i,:,:)));
%             level = graythresh(squeeze(filtBW(i,:,:)));
%             otsuBW(i,:,:) = im2bw(squeeze(filtBW(i,:,:)),level);
%         end
%         figure(99), subplot(1,2,1), imagesc(squeeze(otsuBW(:,:,t3)));
%         subplot(1,2,2), imagesc(squeeze(otsuBW(:,t2,:)));
%         %thresh=0.05;
%         %mask = filledBW>(thresh.*max(filledBW(:)));
%         %figure(9), subplot(1,2,1), imagesc(squeeze(mask(:,:,t3)));
%         %subplot(1,2,2), imagesc(squeeze(mask(:,t2,:)));
%         
% %         filledMask = imfill(mask,'holes');
% %         figure(10), subplot(1,2,1), imagesc(squeeze(filledMask(:,:,t3)));
% %         subplot(1,2,2), imagesc(squeeze(filledMask(:,t2,:)));
% %         closedMask = imclose(double(mask),se2);
% %         figure(10), subplot(1,2,1), imagesc(squeeze(closedMask(:,:,t3)));
% %         subplot(1,2,2), imagesc(squeeze(closedMask(:,t2,:)));
%         
%         
%         Yanat_WOspots = im2d.*mask;
%         spots = im2d-Yanat_WOspots;
%         figure(10), subplot(1,2,1), imagesc(squeeze(Yanat_WOspots(:,:,t3)));
%         subplot(1,2,2), imagesc(squeeze(Yanat_WOspots(:,t2,:)));
%         figure(11), subplot(1,2,1), imagesc(squeeze(spots(:,:,t3)));
%         subplot(1,2,2), imagesc(squeeze(spots(:,t2,:)));
%         % optodesCoord_click(spots,'Spots');
%         
% 
%         
%         
% %         % Median filter + otsu
% %         for i=1:size(spots,2)
% %             spots2(:,i,:) = medfilt2(squeeze(spots(:,i,:)));
% %         end
% %         
% %         for i=1:size(spots,1)
% %             spots2(i,:,:) = medfilt2(squeeze(spots2(i,:,:)));
% %             level = graythresh(squeeze(spots2(i,:,:)));
% %             spots3(i,:,:) = im2bw(squeeze(spots2(i,:,:)),level);
% %         end
% %         
% %         figure(12), subplot(1,2,1), imagesc(squeeze(spots2(:,:,t3)));
% %         subplot(1,2,2), imagesc(squeeze(spots2(:,t2,:)));
% %         figure(13), subplot(1,2,1), imagesc(squeeze(spots3(:,:,t3)));
% %         subplot(1,2,2), imagesc(squeeze(spots3(:,t2,:)));
% %         
% %         % Histograms
% %         [n,xout] = hist(spots(spots>0));
% %         figure(14), bar(xout,n);
% %         
% %         [n2,xout2] = hist(spots2(spots2>0));
% %         figure(15), bar(xout2,n2);
%         
% %         test = imregionalmax(spots,26);
% %         figure(16), subplot(1,2,1), imagesc(squeeze(test(:,:,t3)));
% %         subplot(1,2,2), imagesc(squeeze(test(:,t2,:)));
% %         % Je voudrais les regional max mais avec un seuil sur la taille...
% %         % Ou seuiller l'image avant de
% %         % faire le regional maxima stuff. C'est pareil que imextendedmax...
%         
%         % REGIONAL MAXIMA
%         thresh=0.5;
%         vitaminsBW = imextendedmax(spots,thresh*max(spots(:)),26);
%         figure(17), subplot(1,2,1), imagesc(squeeze(vitaminsBW(:,:,t3)));
%         subplot(1,2,2), imagesc(squeeze(vitaminsBW(:,t2,:)));
% 
%         % Identify centroid
%         stats = regionprops(test3, 'Centroid');
%         %Stats(1).Centroid;
%         for ii=1:size(stats,1)
%             coord(ii,:) = round(stats(ii).Centroid([2 1 3]));
%         end
%         optodesCoord_click(test3,'bob',coord,[]);
%         
% %         thresh=0.01;
% %         test2 = imhmax(spots,thresh*max(spots(:)));
% %         figure(14), subplot(1,2,1), imagesc(squeeze(test2(:,:,t3)));
% %         subplot(1,2,2), imagesc(squeeze(test2(:,t2,:)));
% %         

% Reflexions obsoletes.... et tata
 %  min sqrt((order*coord_fid-coord_helmet).^2) with respect to order of optodes
    % will this converge?? (numOptodes variables optimized... not really%?).
    % where order=a permutation of 1:nOptodes... in fact a rearrangement of
    % the lines of eye(nopt,nopt). (Possible rearrangements of 12 lines)
    % I([order],:)*pos_vit = pos_helmet. where order=randperm(1:12).
    % (Means vitNames = inv(I)*helmet.names)).    
    % x = lsqlin(C,d,A,b) solves the linear system C*x = d in the least-squares sense subject to A*x ? b,
     % optimize wrt order, but not any 12 numbers... its 12 numbers taken
    % out of 1:12 (1 of each). order=randperm(1:12);
    % Fonction to minimize:
    % f = bob(which argument??)
    
    
