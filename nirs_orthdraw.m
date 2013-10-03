function nirs_orthdraw(xSPM, render)
%Ke Peng, 2013-08-28, called in nirs_run_liom_orth_coreg

V_render = spm_vol(render);

%Transfer mm -> voxel
tSPM.tXYZmm = [xSPM.XYZmm;ones(1,size(xSPM.XYZmm,2))];
tSPM.tXYZ = V_render.mat\tSPM.tXYZmm;

xSPM.XYZ = tSPM.tXYZ;
tSPM.Z = xSPM.Z;
tSPM.label = xSPM.label;
%xSPM.XYZ = tXYZ(1:3,:);

%Interpolation, as the points are too small

for x0 = -3 : 3
    for y0 = -3 : 3
        for z0 = -3 : 3
            if x0 == 0 && y0 == 0 && z0 == 0
                continue;
            else
                new_points = [tSPM.tXYZ(1,:)+x0; tSPM.tXYZ(2,:)+y0; tSPM.tXYZ(3,:)+z0; tSPM.tXYZ(4,:)];
                xSPM.XYZ = [xSPM.XYZ new_points];
                xSPM.Z = [xSPM.Z tSPM.Z];
                xSPM.label = [xSPM.label tSPM.label];
            end
        end
    end
end

%Transfer back to mm

xSPM.XYZmm = [];
xSPM.XYZmm = V_render.mat * xSPM.XYZ;
xSPM.XYZ = xSPM.XYZ(1:3,:);
xSPM.XYZmm = xSPM.XYZmm(1:3,:);

xSPM.M = V_render.mat;
xSPM.DIM = V_render.dim;

%Draw
nirs_orthviews('Reset');
h = nirs_orthviews('Image', render, [0.05 0.05 0.9 0.7]);
nirs_orthviews('AddContext', h);
nirs_orthviews('MaxBB');
%if ~isempty(hReg), spm_orthviews('Register', hReg); end
nirs_orthviews('AddBlobs', h, xSPM.XYZ, xSPM.Z, xSPM.M);
nirs_orthviews('Redraw');