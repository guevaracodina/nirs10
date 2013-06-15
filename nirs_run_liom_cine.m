function out = nirs_run_liom_cine(job)
%Cine
Z = get_contrast_group_common_options(job);
%Views:
Z.views_to_run = job.view;
% 1: 'ventral'
% 2: 'dorsal'
% 3: 'right_lateral'
% 4: 'left_lateral'
% 5: 'frontal'
% 6: 'occipital'
Z.DoStats = 1; 
Z.erdf_default = 100; %will give thresholds in the limit of large # of degrees of freedom
Z.LKC = 0;
Z.UseCorrelRes = 0;
if Z.LKC
    Z.StatStr = 'EC';
else
    Z.StatStr = 'Tube';
    Z.StatStr2 = 'Bonf';
end

Z.automated_contrasts = 1;
Z.GroupMultiSession = 0;
Z.NonlinearEpilepsyOn = 0;
%Z.Avg = 1;
%W.Avg = Z.Avg;
Z.output_unc = 1;

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
Z.InterpolationKernel = job.InterpolationKernel;
Z.select_chromophore = job.select_chromophore;
Z.onsetInfo = job.onsetInfo;
Z.CineFilters = job.CineFilters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loop over all subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Idx=1:size(job.NIRSmat,1)
    %Load NIRS.mat information
    try
        run_cine_OK = 1;
        [NIRS newNIRSlocation]= nirs_load(job.NIRSmat{Idx,1},job.NIRSmatCopyChoice,job.force_redo);
        job.NIRSmat{Idx,1} = newNIRSlocation;
        if ~isfield(NIRS,'flags'), NIRS.flags = []; end
        if ~isempty(NIRS) && (~isfield(NIRS.flags,'cine_OK') || job.force_redo)
            NC = NIRS.Cf.H.C.N;
            [rendered_MNI run_cine_OK NIRS] = nirs_load_TopoData(job,NIRS,run_cine_OK);           
            if run_cine_OK
                %load CINE (topographic maps) if already (partially) generated
                [CINE newDir fCINE] = nirs_load_TOPO(newNIRSlocation);
                %Careful, this is one aspect of Z that is subject specific
                Z.dir1 = newDir;
                Z.Idx = Idx;                
                if isfield(NIRS.Dt.s,'subj_id')
                    Z.subj_id = NIRS.Dt.s.subj_id;
                end
                %Big loop over views
                cv = 0; %view counter
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
                            CINE.v{W.side_hemi}.Warning = 'No channel found for this view';
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
                            cv = cv + 1;
                            if ~isempty(Z.InterpolationKernel{1})
                                load(Z.InterpolationKernel{cv});
                            else
                                %Q = interpolation_kernel(q0,W,Z.LKC,Z.UseCorrelRes,Z.DoStats);
                                Q = interpolation_kernel_cine(W);
                                Qname = fullfile(newDir,['Q_' W.spec_hemi '.mat']);
                                save(Qname,'Q');
                            end
                            W.Q = Q;                           
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            CINE = prepare_cine_core_call(Z,W,NIRS,CINE);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            CINE.v{W.side_hemi}.s1 = W.s1; %sizes of topographic projection
                            CINE.v{W.side_hemi}.s2 = W.s2;
                            CINE.v{W.side_hemi}.view = W.spec_hemi; %%% view of the brain
                        end
                    catch exception
                        disp(exception.identifier);
                        disp(exception.stack(1));
                        disp(['Could not create cine for view ' W.spec_hemi ' for subject ' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
                    end
                end %end for v1
                save(fCINE,'CINE','-v7.3'); %file can be large - because it stores all the contrast data
            end
            NIRS.flags.con_OK = 1;
            NIRS.CINE = fCINE;
            save(job.NIRSmat{Idx,1},'NIRS');
        end
    catch exception
        disp(exception.identifier);
        disp(exception.stack(1));
        disp(['Could not create cine for subject' int2str(Idx) ' for ' job.NIRSmat{Idx,1}]);
    end
end
out.NIRSmat = job.NIRSmat;