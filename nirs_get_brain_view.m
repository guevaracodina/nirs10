function [side_hemi spec_hemi] = nirs_get_brain_view(brain_view)
switch brain_view
    case 1 % 'ventral'
        spec_hemi = 'ventral';
        side_hemi = 1;
    case 2 % 'dorsal'
        spec_hemi = 'dorsal';
        side_hemi = 2;
    case 3 %'right_lateral'
        spec_hemi = 'right';
        side_hemi = 3;
    case 4 %'left_lateral'
        spec_hemi = 'left';
        side_hemi = 4;
    case 5 %'frontal'
        spec_hemi = 'frontal';
        side_hemi = 5;
    case 6 %'occipital'
        spec_hemi = 'occipital';
        side_hemi = 6;
end