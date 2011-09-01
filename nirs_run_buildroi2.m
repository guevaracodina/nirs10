function out = nirs_run_buildroi2(job)
% Builds the smallest ROI containing all the selected source-detector pairs
% FORMAT NIRSmat = nirs_run_buildroi2(NIRS, image_in, keepAllChannels, outputprefix)
% NIRS            - NIRS matrix
% image_in        - segmented image from which the ROI will be extracted
% keepAllChannels - Channels to keep in the ROI
% outputprefix    - Prefix of the ROI image name
%_______________________________________________________________________
%
% The image will be cropped with respect to the selected source-detector
% pairs (channels) to be kept for the MonteCarlo simulation. You
% only have to enter the channel numbers for the first wavelength.
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire

% Clément Bonnéry
% 2011-03

try
    NewNIRSdir = job.NewDirCopyNIRS.CreateNIRSCopy.NewNIRSdir;
    NewDirCopyNIRS = 1;
catch
    NewDirCopyNIRS = 0;
end

for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        NIRS = [];
        load(job.NIRSmat{Idx,1});
        [dir0,dummy,dummy2] = fileparts(job.NIRSmat{Idx,1});
        
        if ~isempty(job.image_in{:})
            V = spm_vol(job.image_in{:});            
        else
            V = spm_vol(NIRS.Dt.ana.T1seg);
        end
        Y = spm_read_vols(V);
        [dir,name] = fileparts(V.fname);
        
        NS = NIRS.Cf.H.S.N;
        Pfp_rmm = NIRS.Cf.H.P.r.m.mm.fp;
        Pp_rmm = NIRS.Cf.H.P.r.m.mm.p;
        Pp_c1_rmm = NIRS.Cf.H.P.r.m.mm.c1.p;
        Cid = NIRS.Cf.H.C.id;
        NC = NIRS.Cf.H.C.N;
        wl = NIRS.Cf.dev.wl;
        nc = NC/length(wl);
        
        try
            Ckpt = job.keepAllChannels.keepChannels;
            NewNIRSdir = ['roi_channels-' strrep(int2str(job.keepAllChannels.keepChannels),'  ','-')];
        catch
            Ckpt = Cid(1,:);
            NewNIRSdir = 'roi_all-channels';
        end
        
        Ckpt_owl = [];
        Skpt =[];
        Dkpt=[];
        %Loop over channels chosen by user
        for i=1:length(Ckpt)
            %keep track of channel location
            Cbloup = (1:length(Cid(1,:))).*(Cid(1,:) == Ckpt(i));
            Skpt = [Skpt Cid(2,sum(Cbloup))];
            Dkpt = [Dkpt Cid(3,sum(Cbloup))];
            Ckpt_owl = [Ckpt_owl Cid(1,Cid(2,:)==Cid(2,sum(Cbloup)) & Cid(3,:)==Cid(3,sum(Cbloup)))];%other wavelengths
        end
        
        Ckpt = unique([Ckpt Ckpt_owl]);
        Skpt = unique(Skpt);
        Dkpt = unique(Dkpt);
        Pkpt=[Skpt Dkpt+NS];
        
        if NewDirCopyNIRS
            dir2 = [dir0 filesep NewNIRSdir];
            if ~exist(dir2,'dir'), mkdir(dir2); end;
        else
            dir2 = dir0;
        end
        
        %%%  vérifier que ça a de l'intérêt de garder tout ça
        Pfp_roi_rmm     = zeros(3,size(Pkpt,2));
        Pp_roi_rmm      = zeros(3,size(Pkpt,2));
        Pp_roi_c1_rmm   = zeros(3,size(Pkpt,2));
        %4th row in order to multiply by V.mat
        PfpR_roi_rmv    = zeros(4,size(Pkpt,2));
        Pfp_roi_rmvtemp = zeros(4,size(Pkpt,2));
        
        for i=1:size(Pkpt,2)
            Pfp_roi_rmm(:,i) = Pfp_rmm(:,Pkpt(i));
            Pfp_roi_rmvtemp(:,i) = V.mat\[Pfp_roi_rmm(:,i);1];
            
            Pp_roi_rmm(:,i) = Pp_rmm(:,Pkpt(i));
            Pp_roi_c1_rmm(:,i) = Pp_c1_rmm(:,Pkpt(i));
        end
        Pfp_roi_rmv = Pfp_roi_rmvtemp(1:3,:);
        
        method = 'S_cubecenter';
        switch method
            case 'ancienne'
%                 bbv(1,1) = min(Pfp_roi_rmv(1,:));
%                 bbv(1,2) = max(Pfp_roi_rmv(1,:));
%                 
%                 bbv(2,1) = min(Pfp_roi_rmv(2,:));
%                 bbv(2,2) = max(Pfp_roi_rmv(2,:));
%                 
%                 bbv(3,1) = min(Pfp_roi_rmv(3,:));
%                 bbv(3,2) = max(Pfp_roi_rmv(3,:));
                
                bbv(:,1) = min(Pfp_roi_rmv,[],2);
                bbv(:,2) = max(Pfp_roi_rmv,[],2);
                bbv = round(bbv); %?
                marge = 40;%5
                
            case 'S_cubecenter'
                %additional 4 cm (40 mm) added to source 
                ray = 40;%max(NIRS.Cf.H.C.gp(Ckpt)+3);
                bbm = zeros(3,length(Skpt));
                bbM = zeros(3,length(Skpt));
                for iS=1:length(Skpt);
                    % cube a 4cm du centre S
                    S_cc = Pfp_roi_rmm(:,iS);
                    bbm(:,iS) = S_cc-ray*ones(3,1);
                    bbM(:,iS) = S_cc+ray*ones(3,1);
                    % on le coupe la ou il depasse trop
                    %A FAIRE
                    bbmv1 = V.mat\[bbm;1];
                    bbMv1 = V.mat\[bbM;1];
                    %min and max of 
                    bbmv(1:3,iS) = round(bbmv1(1:3));
                    bbMv(1:3,iS) = round(bbMv1(1:3));
                end
                bbv(1,1) = min(bbmv(1,:),bbMv(1,:));
                bbv(1,2) = max(bbmv(1,:),bbMv(1,:));
                
                bbv(2,1) = min(bbmv(2,:),bbMv(2,:));
                bbv(2,2) = max(bbmv(2,:),bbMv(2,:));
                
                bbv(3,1) = min(bbmv(3,:),bbMv(3,:));
                bbv(3,2) = max(bbmv(3,:),bbMv(3,:));
                marge =0;
        end
        % the size of the plotted image can be bigger than the size read in
        % the header, in such case the value kept is the one of header
        for i =1:3
            bbv(i,1) = max(1,bbv(i,1)-marge); 
            bbv(i,2) = min(V.dim(i),bbv(i,2)+marge);
        end
        
        % With resizing
        %  mat_roi2raw is the transformation matrix from the voxel space
        %  of the ROI to the voxel space of the raw image
        % Translation : de l'image brute a la roi
        mat_roi2raw = eye(4) + [zeros(4,3) [bbv(:,1);0]];
        
        V_roi.dim = [bbv(1,2)-bbv(1,1)+1 ...
                     bbv(2,2)-bbv(2,1)+1 ...
                     bbv(3,2)-bbv(3,1)+1];
        % V_roi.mat is the mat from voxel ROI to the raw image space in mm
        V_roi.mat = V.mat*mat_roi2raw; %just a translation
        %select data in ROI
        Y_roi = Y(bbv(1,1):bbv(1,2),bbv(2,1):bbv(2,2),bbv(3,1):bbv(3,2));
       
        V_roi = nirs_create_vol(fullfile(dir2,[job.output_prefix,name,'.nii']),...
                    V_roi.dim, V.dt, V.pinfo, V_roi.mat, Y_roi);               
        
        % Sources and detectors positions must be updated also
        %voxel positions only are updated (real positions are unchanged)
        for i=1:size(Pkpt,2)
            PfpR_roi_rmv(:,i) = [Pfp_roi_rmv(:,i);1] - mat_roi2raw(:,4);
        end
        
        %then saved
        if ~isfield(NIRS,'Cs'), NIRS.Cs={}; end
        if isfield(NIRS.Cs,'temp'), clear NIRS.Cs.temp; end
        NIRS.Cs.temp.Ckpt = Ckpt; %kept channels
        NIRS.Cs.temp.Pkpt = Pkpt; %kept points (excluding channels)
        NIRS.Cs.temp.NSkpt = size(Skpt,2);
        NIRS.Cs.temp.NDkpt = size(Dkpt,2);
        %Fitted position of points in MNI voxel space
        NIRS.Cs.temp.Pfp_roi_rmv = PfpR_roi_rmv(1:3,:); 
        %Same but in MNI real space
        NIRS.Cs.temp.Pfp_roi_rmm = Pfp_roi_rmm;
        NIRS.Cs.temp.Pp_roi_rmm = Pp_roi_rmm;
        NIRS.Cs.temp.Pp_roi_c1_rmm = Pp_roi_c1_rmm;
        NIRS.Cs.temp.segR = fullfile(dir2,[job.output_prefix,name,'.nii']);
        
        if NewDirCopyNIRS
            newNIRSlocation = fullfile(dir2,'NIRS.mat');
            save(newNIRSlocation,'NIRS');
            job.NIRSmat{Idx,1} = newNIRSlocation;
        else
            save(job.NIRSmat{Idx,1},'NIRS');
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not run build roi for subject' int2str(Idx)]);
    end
end
out.NIRSmat = job.NIRSmat;