function out = nirs_run_view3d(job)
% View optodes en 3D
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
%
% based on Michèle Desjardins version of prepMCsim2
%
% Clément Bonnéry
% 2010-06

% Loop over subjects
for iSubj=1:size(job.NIRSmat,1)
    
    % Load NIRS.mat
    NIRS = [];
    try
        load(job.NIRSmat{iSubj,1});
    catch
        disp(['Could not load NIRS.mat for ' int2str(Idx) 'th subject in nirs_run_view3d.']);
    end

    if isempty(job.segT1_4fit{1,1})
        try
            if ~job.on_cortex
                segT1_4fit = NIRS.Dt.ana.T1seg;
            else
                % To implement: option to display optodes on cortex
                %try segT1_4fit = ['c3' NIRS.Dt.ana.T1];...
                disp('Need to specify ''c3'' image in order to display optodes on cortex in nirs_run_view3d.');
                disp('Optodes will be displayed on scalp instead.');
                segT1_4fit = NIRS.Dt.ana.T1seg;
            end
        catch
            disp(['Could not find a segmented image to display positions on for ' int2str(iSubj) 'th subject.']);
        end
    else
        try
            segT1_4fit = job.segT1_4fit{iSubj,:};
        catch
            segT1_4fit = NIRS.Dt.ana.T1seg;
        end
    end
    
    V = spm_vol(segT1_4fit);
    Y = spm_read_vols(V);

    if isfield(NIRS.Cf.H.P.r.m.mm,'p')
        Pp_rmm = NIRS.Cf.H.P.r.m.mm.p;% positions on head surface
        Pp_c1_rmm = NIRS.Cf.H.P.r.m.mm.c1.p;% positions on cortex
        Pfp_rmm = NIRS.Cf.H.P.r.m.mm.fp;% positions fitted on skin
    end

    NP = NIRS.Cf.H.P.N;
    NS = NIRS.Cf.H.S.N;
    ND = NIRS.Cf.H.D.N;

    % % Fitted positions (positions are fitted with respect to the scalp)
    % jobe.NIRS = NIRS;
    % jobe.image_in = segT1_4fit;
    % out = nirs_fit_probe(jobe);
    % Pfp_rmm = out{1};

    % from MNI real space (mm) to MNI voxel space
    Pp_rmv = [];
    Pp_c1_rmv = [];
    Pfp_rmv = [];
    for i=1:NP
        Pp_rmv(:,i) = V.mat\[Pp_rmm(:,i);1];
        Pp_c1_rmv(:,i) = V.mat\[Pp_c1_rmm(:,i);1];
        Pfp_rmv(:,i) = V.mat\[Pfp_rmm(:,i);1];
    end
    % return to 3D space (V.mat is (:,4))
    Pp_rmv = Pp_rmv(1:3,:);
    Pp_c1_rmv = Pp_c1_rmv(1:3,:);

    %% Display 3D image of the ROI with optodes pointing towards their
    % respective directions, allowing a final check of the setup
    % Surface de l'IRM
    set(0,'defaultfigurevisible','on');
    hfig_all3D = figure;

    patch(isosurface(smooth3(Y),4),...
        'FaceColor',[1,.75,.65],...
        'EdgeColor','none',...
        'FaceAlpha',0.5);

    % Views adjustments
    view(0,0)
    daspect([1,1,1])
    lightangle(70,-70);
    set(hfig_all3D,'Renderer','zbuffer');
    lighting phong

    % Display optodes and their directions
    hold on

    % adding tags
    for Pi = 1:NP
        if Pi<=NS
            list{1,1} = 'r';
            list{2,1} = 10;
            list{3,1} = 'S#';
            Pinum = int2str(Pi);
        elseif (Pi>NS && Pi<=NS+ND)
            list{1,1} = 'b';
            list{2,1} = 2;
            list{3,1} = 'D#';
            Pinum = int2str(Pi-NS);
        elseif Pi>NS+ND
            list{1,1} = 'g';
            list{2,1} = 2;
            list{3,1} = 'Q#';
            Pinum = int2str(Pi-(NS+ND));
        end

        xp = Pp_rmv(2,Pi);
        yp = Pp_rmv(1,Pi);
        zp = Pp_rmv(3,Pi);
        text(xp,yp,zp,[list{3,1} Pinum],'FontWeight','bold','Color',list{1,1});

        xfp = Pfp_rmv(2,Pi);
        yfp = Pfp_rmv(1,Pi);
        zfp = Pfp_rmv(3,Pi);
        text(xfp,yfp,zfp,'X','FontWeight','bold','Color','black');%list{1,1}

        plot3([Pp_c1_rmv(2,Pi),Pp_rmv(2,Pi)],...
            [Pp_c1_rmv(1,Pi),Pp_rmv(1,Pi)],...
            [Pp_c1_rmv(3,Pi),Pp_rmv(3,Pi)], 'Linewidth',list{2,1},'Color',list{1,1});
    end
    hold off

    clear NIRS
    load(job.NIRSmat{iSubj,:});
    NIRS.Cf.H.P.r.m.mm.fp = Pfp_rmm;
    NIRS.Cf.H.P.r.m.vx.fp = Pfp_rmv(1:3,:);
    save(job.NIRSmat{iSubj,:},'NIRS');

    % Save and close figure, or not, according to user-interface defined option
    if job.save_figure
        saveas(gcf,fullfile(NIRS.Dt.s.p,'coreg_3Dview.fig'));
        close(gcf);
    end
    
end

out =1;
return