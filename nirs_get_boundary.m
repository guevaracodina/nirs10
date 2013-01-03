function W = nirs_get_boundary(W,job)
if isfield(job.display_options,'show_boundary')
    show_boundary = job.display_options.show_boundary;
else
    show_boundary = 1;
end
try
    if show_boundary
        [bx by] = size(W.brain_view_mask_2d);
        W.boundary = [diff(W.brain_view_mask_2d,1,1); zeros(1,by)];
        W.boundary = abs(W.boundary) + abs([diff(W.brain_view_mask_2d,1,2) zeros(bx,1)]);
        W.boundary_mask = logical(W.boundary);
        W.brain(W.boundary_mask) = 0;
    end
catch
    disp('Cannot add boundary mask -- check brain_view_mask_2d')
end