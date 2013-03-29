function W = nirs_get_common_brain_mask(W,big_TOPO,v1)
if ~isfield(W,'ignore_brain_mask')
    W.ignore_brain_mask = 0;
end
if ~isempty(big_TOPO) && isfield(W,'brain_view_mask_2d') && ~W.ignore_brain_mask
    for s0 = 1:length(big_TOPO)
        W.brain_view_mask_2d = W.brain_view_mask_2d .* big_TOPO{s0}.rendered_MNI{v1}.view_mask_2d; 
    end
end