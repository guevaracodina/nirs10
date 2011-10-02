function fY = nirs_resize_no_save(Y,dimf)
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Clément Bonnéry 11/2010

%Problem with orientation of images

dimi = size(Y);
% Create linear spaces for initial and final volume
for i=1:3
x{i} = linspace(0,dimf(i)-1,dimi(i));
f{i} = 0:1:dimf(i)-1;
end
% Make grids, to use with 3-dimensional interpolation
[xi1,yi1,zi1] = meshgrid(x{1},x{2},x{3});
[xf1,yf1,zf1] = meshgrid(f{1},f{2},f{3});
Y = permute(Y ,[2,1,3]); %This is a peculiarity due to interp3
fY = interp3(xi1,yi1,zi1,Y,xf1,yf1,zf1,'nearest');
fY= permute(fY,[2,1,3]); %Permute it back
end