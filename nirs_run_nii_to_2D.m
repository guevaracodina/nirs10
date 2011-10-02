function out = nirs_run_nii_to_2D(job)
%loop over specified files and folders
ff = job.convert_files;
%Define W
W = [];
W.gen_fig = 0;
W.gen_tiff = 1;
for f1=1:length(ff)
    file = ff{f1};
    if spm_existfile(file) 
        W = do_convert(file,W);        
    else
        file = ff{f1}(1:end-2);
        if spm_existfile(file) 
            W = do_convert(file,W);
        else
            disp(['Could not convert: ' ff{f1}]);
        end
    end
end
ff = job.convert_folders;
for f1=1:length(ff)
    if exist(ff{f1},'dir')
        [files dirs] = spm_select('FPList',ff{f1},'.*');
    end
    for f2=1:size(files,1)
        W = do_convert(files(f2,:),W);
    end
end
out = [];
end

function W = do_convert(file,W)
[dir fil ext] = fileparts(file);
%loop through the 6 views - inverse function compared to nirs_get_brain_view
views = {'ventral','dorsal','right','left','frontal','occipital'};
side_hemi = 0;
for v1=1:length(views)
    if strfind(fil,views{v1})
        side_hemi = v1;
        spec_hemi = views{v1};
    end
end
W.askDOF = 1;
if side_hemi

    if askDOF
        [null,YPos]=spm_input('Enter one # for T-stat, 2 for F-stat',1,'d!','Degrees of freedom');
        [null,YPos]=spm_input('enter R-C (# of rows less columns in design matrix)','+1','d!','T');
        [null,YPos]=spm_input('write R-C, a space, then write number of contrasts','+1','d!','F');        
        p = spm_input('Enter degrees of freedom:','+1');
        p2 = spm_input('Apply to all images','+1','y/n');
        p3 = spm_input('Enter p value','+1');
        if strcmp(p2,'y')
            W.askDOF = 0;
        end
        W.p_value = p3;
        W.erdf = p(1);
        if length(p) == 2
            W.eidf = p(2);
        end
    end

    %Define W
    
    %Define F
    load Split
    F.split = split;
    F.ftiff = fullfile(dir,[fil '.tiff']);
    F.ffig = fullfile(dir,[fil '.fig']);
    F.erdf;
    F.eidf;
    F.tstr;
    F.T_map;
    F.hb;
    nirs_simple_figure(F,W);       
else
    disp(['View not recognized for: ' file]);
end
end
