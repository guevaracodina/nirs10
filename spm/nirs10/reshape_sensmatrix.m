function out = reshape_sensmatrix(sens,Volbase)

Vbase = spm_vol(Volbase);
sens_reshaped = zeros(Vbase.dim);
for i=1:size(sens,1)
    sens_reshaped = sens_reshaped + reshape(sens(i,:),Vbase.dim);
end

V.dim = [99 95 80];%V.dim = [ny nx nz];
V.mat = Vbase.mat;
V.dt = [4,0];
V.pinfo = [1;0;352];

V = struct('fname',fullfile('D:\Users\Clément\Projet_ReML\donnees\test_roi\MCconfigHopingpongPINGPOUNG',['sensReconstruct' '.nii']),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',V.mat);

V = spm_create_vol(V);
V = spm_write_vol(V, sens_reshaped);

out =1;
end
