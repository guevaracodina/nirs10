function [bool, energie_quantifiee] = nirs_runMOB_walktest(fft_slab,fft_size)
% Cl�ment Bonn�ry

bool=0;

energie = fft_slab.^2;
energie_roi = energie(30:fft_size/2);

energie_quantifiee = sum(energie_roi);


% tester le seuil
if energie_quantifiee>1.4*10^8;
    bool =1;
end