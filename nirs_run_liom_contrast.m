function out = nirs_run_liom_contrast(job)
%Contrasts are generated in one of 2 ways:
%1- If no contrast is user-specified, automated contrasts are generated
%2- If contrasts are user-speficied, automated contrasts are not generated
%Contrasts need to be put in the SPM xCon structures
%If there is more than one session, contrasts can be run by session or
%over all sessions - the latter is the SPM standard way
%Fill Z Structure, for passing the most generic data
Z = get_contrast_group_common_options(job);
%Views:
Z.views_to_run = job.view;
% 1: 'ventral'
% 2: 'dorsal'
% 3: 'right_lateral'
% 4: 'left_lateral'
% 5: 'frontal'
% 6: 'occipital'
Z.DoStats = job.GenerateStats; %Boolean; if set to false, will skip LKC, statistics and maps 
Z.LKC = job.StatMethod; %Use Lipschitz-Killing curvature
Z.UseCorrelRes = job.UseCorrelRes;
if Z.LKC
    Z.StatStr = 'EC';
else
    Z.StatStr = 'Tube';
    Z.StatStr2 = 'Bonf';
end

%SPM contrasts or automatic contrasts
if isfield(job.ContrastChoice,'user_contrasts')
    if isempty(job.ContrastChoice.user_contrasts.consess)
        disp('No contrast specified. Aborting. Either specify contrasts or use automated-contrasts option');
        return
    else
        Z.automated_contrasts = 0;
        TF.consess = job.ContrastChoice.user_contrasts.consess;
        Z.GroupMultiSession = job.ContrastChoice.user_contrasts.GroupMultiSession;
    end
else
    TF.consess = [];
    %automated contrasts if no user-specified contrast
    Z.automated_contrasts = 1;
    Z.GroupMultiSession = 0;
    %Settings for epilepsy study on nonlinearities only
    Z.NonlinearEpilepsyOn = job.ContrastChoice.automated_contrasts.NonlinearEpilepsyOn;
end
%Gaussian spatial LPF
if isfield(job.spatial_LPF,'spatial_LPF_On')
    Z.radius = job.spatial_LPF.spatial_LPF_On.spatial_LPF_radius;
    Z.spatial_LPF = 1;
else
    Z.spatial_LPF = 0;
end
if ~isempty(job.Sessions)
    Z.sessions = job.Sessions;
else
    Z.sessions = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loop over all subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        run_contrast_OK = 1;
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isfield(NIRS,'flags'), NIRS.flags = []; end
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'con_OK') || job.force_redo)
            NC = NIRS.Cf.H.C.N;
            [rendered_MNI run_contrast_OK NIRS] = nirs_load_TopoData(job,NIRS,run_contrast_OK);
            %load SPM - first GLM - might want to generalize
            try
                [dir1 dummy] = fileparts(NIRS.SPM{1});
            catch exception
                disp(exception.identifier);
                disp(exception.stack(1));
                disp('Could not find SPM file as specified in NIRS.mat')
            end
            try 
                load(NIRS.SPM{end});
            catch
                load(fullfile(NIRS.SPM{end},'SPM.mat')); %old format
            end
            if ~exist('SPM','var')
                run_contrast_OK = 0;
                disp('SPM not found');
            end
            if run_contrast_OK
                %load TOPO (topographic maps) if already (partially) generated
                [TOPO newDir fTOPO] = nirs_load_TOPO(newNIRSlocation);
                %Careful, this is one aspect of Z that is subject specific
                Z.dir1 = newDir;
                Z.Idx = Idx;
                try
                if isfield(NIRS.Dt.s,'subj_id')
                    Z.subj_id = NIRS.Dt.s.subj_id;
                end
                end
                %Calculate xCon and put it in TOPO -- also may need SSxCon
                TOPO = nirs_get_contrasts(SPM,Z,TF,TOPO);
                %Calculate required information on F-stats, which is
                %independent of the view
                if Z.LKC || Z.UseCorrelRes
                    TOPO = precalculate_F(SPM,Z,TOPO);
                end
                %Big loop over views
                for v1=1:size(Z.views_to_run,2)
                    try
                        %Fill W structure, which is view and subject specific
                        %Structure for passing GLM and interpolation data
                        clear W
                        brain_view = Z.views_to_run(v1);
                        [W.side_hemi W.spec_hemi] = nirs_get_brain_view(brain_view);
                        Z.spec_hemi = W.spec_hemi; %needed to initialize assembled figures
                        % channel information
                        rchn = rendered_MNI{W.side_hemi}.rchn;
                        cchn = rendered_MNI{W.side_hemi}.cchn;
                        W.AllowExtrapolation = job.AllowExtrapolation;
                        W.no_interpolation = job.no_interpolation;
                        %find channels which are visible from this projection view
                        W.index_ch = find(rchn ~= -1);
                        if isfield(rendered_MNI{W.side_hemi},'view_mask_2d') % for back-compatibility
                            W.brain_view_mask_2d = rendered_MNI{W.side_hemi}.view_mask_2d;                            
                        end
                        
                        if isempty(W.index_ch)
                            TOPO.v{W.side_hemi}.Warning = 'No channel found for this view';
                            if v1 > 1 && v1 < 6
                                disp(['No channel for view ' int2str(v1) ': Probable coregistration problem. Skipping this view']);
                            end
                        else
                            %rendering surface
                            brain = rendered_MNI{W.side_hemi}.ren;                                                        
                            W.brain = brain * 0.5;
                            W = nirs_get_boundary(W,job);
                            W.s1 = size(brain, 1);
                            W.s2 = size(brain, 2);
                            %split into HbO and HbR interpolations
                            W.ch_HbO = W.index_ch;
                            W.ch_HbR = NC/2 + W.index_ch;
                            W.ch_HbT = NC+W.index_ch;
                            W.rchn = rchn(W.index_ch);
                            W.cchn = cchn(W.index_ch);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %Generate the contrasts
                            TOPO = prepare_constrast_core_call(Z,W,SPM,TOPO);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            TOPO.v{W.side_hemi}.s1 = W.s1; %sizes of topographic projection
                            TOPO.v{W.side_hemi}.s2 = W.s2;
                            TOPO.v{W.side_hemi}.view = W.spec_hemi; %%% view of the brain
                        end
                    catch exception
                        disp(exception.identifier);
                        disp(exception.stack(1));
                        disp(['Could not create contrasts for view ' W.spec_hemi ' for subject ' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
                    end
                end %end for v1
                save(fTOPO,'TOPO','-v7.3'); %file can be large - because it stores all the contrast data
            end
            NIRS.flags.con_OK = 1;
            NIRS.TOPO = fTOPO;
            save(job.NIRSmat{Idx,1},'NIRS');
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not create contrasts for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
    end
end
out.NIRSmat = job.NIRSmat;