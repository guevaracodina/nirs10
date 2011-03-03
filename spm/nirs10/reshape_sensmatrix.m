function out = reshape_sensmatrix(sens,Volbase)

% lambda = [830 690];
% 
% t =sens;
% %Load sensitivity matrix for one wavelength
% % [t,sts] = spm_select([1 inf],'.mat','Select sensitivity matrix for one wavelength (sens - .mat file in T1\dataMC\)',[],0);
% % if ~sts, return; end
% load(t); %very long
% 
% %See which wavelength the user has selected
% [dir1 file1 ext1] = fileparts(t);
% for Idx=1:2
%     if strfind(file1,int2str(lambda(Idx))), wlength = Idx; end
% end
% 
% meas = reshape(meas,2*n_chn,[]);
% npt = size(meas,2);
% %select first or second half of measures based on the chosen wavelength
% meas = meas(1+(wlength-1)*n_chn:wlength*n_chn,:);
% %normalize measures
% for Idx=1:n_chn
%     meas(Idx,:) = meas(Idx,:)/meas(Idx,1) -ones(1,npt);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%
% sens_reshaped = zeros(99,95,80);
% for i=1:size(sens,1)
%     sens_reshaped = sens_reshaped + reshape(sens(i,:),[99 95 80]);
% end
sens_reshaped = reshape(sens(1,:),[95 99 80]);
Vbase = spm_vol(Volbase);
V.dim = [95 99 80];%V.dim = [ny nx nz];
V.mat = Vbase.mat;
V.dt = [4,0];
V.pinfo = [1;0;352];

V = struct('fname',fullfile('D:\Users\Clément\Projet_ReML\donnees\test_roi\MCconfigHopingpongPINGPOUNG',['sens1' '.nii']),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',V.mat);

V = spm_create_vol(V);
V = spm_write_vol(V, sens_reshaped);


out =1;
end
