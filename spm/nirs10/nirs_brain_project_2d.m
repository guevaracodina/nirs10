function nirs_brain_project_2d(NIRS,dir_coreg,rendered_MNI,rendered_MNI2,color1,color2,suffix,Pvoid)
N1 = size(rendered_MNI{1}.rchn,1);
if isempty(rendered_MNI2)
    sdmode = 0;
else
    sdmode = 1;
    N2 = size(rendered_MNI2{1}.rchn,1);
end
load Split
if ~isempty(suffix)
    suffix = ['_' suffix];
end
for kk=1:6
    fh0(kk) = figure;
    brain = rendered_MNI{kk}.ren;
    brain = brain.* 0.5;
    sbar = linspace(0, 1, 128);
    sbrain = ((-sbar(1) + sbar(64))/(0.5)).* brain + sbar(1);
    sbrain(1,1) = 1;
    rchn = rendered_MNI{kk}.rchn;
    cchn = rendered_MNI{kk}.cchn;
    for jj = 1:N1
        if isempty(Pvoid) || ~Pvoid(jj)
            if rchn(jj) ~= -1 && cchn(jj) ~= -1
                if rchn(jj) < 6 || cchn(jj) < 6
                    sbrain(rchn(jj), cchn(jj)) = 0.9; % 0.67
                else
                    sbrain(rchn(jj)-5:rchn(jj)+5, cchn(jj)-5:cchn(jj)+5) = 0.9;
                end
            end
        end
    end
    if sdmode
        %copy above
        rchn2 = rendered_MNI2{kk}.rchn;
        cchn2 = rendered_MNI2{kk}.cchn;
        for jj = 1:N2
            if isempty(Pvoid) || ~Pvoid(jj+N1)
                if rchn2(jj) ~= -1 && cchn2(jj) ~= -1 %% updated 2009-02-25
                    if rchn2(jj) < 6 || cchn2(jj) < 6
                        sbrain(rchn2(jj), cchn2(jj)) = 0.9; % 0.67
                    else
                        sbrain(rchn2(jj)-5:rchn2(jj)+5, cchn2(jj)-5:cchn2(jj)+5) = 0.9;
                    end
                end
            end
        end
    end
    imagesc(sbrain);
    colormap(split);
    axis image;
    axis off;
    for jj = 1:N1
        if isempty(Pvoid) || ~Pvoid(jj)
            if rchn(jj) ~= -1 && cchn(jj) ~= -1
                text(cchn(jj)-5, rchn(jj), num2str(jj), 'color', color1);
            end
        end
    end
    if sdmode
        for jj = 1:N2
            if isempty(Pvoid) || ~Pvoid(jj+N1)
                if rchn2(jj) ~= -1 && cchn2(jj) ~= -1
                    text(cchn2(jj)-5, rchn2(jj), num2str(jj), 'color', color2);
                end
            end
        end
    end
    try
        [side_hemi spec_hemi] = nirs_get_brain_view(kk);
        subj_str = '';
        try
            if isfield(NIRS.Dt.s,'subj_id')
                subj_str = [NIRS.Dt.s.subj_id '_'];
            end
        end
        filen = fullfile(dir_coreg,[subj_str spec_hemi suffix]);
        saveas(fh0(kk),[filen '.fig'],'fig');
        print(fh0(kk), '-dpng', [filen '.png'],'-r300');
        close(fh0(kk));
    end
end