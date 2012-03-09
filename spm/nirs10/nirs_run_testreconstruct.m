function out = nirs_run_testreconstruct(job)
% fait avec le sujet 003 et les 9 premiers canaux
sDtp = 'D:\Users\Clément\Projet_ReML\donnees\fantome_ptitet\S002';
% if exist([sDtp '\phantom'],'dir')==7
% else
% mkdir(sDtp,'phantom');
% end

% homogeneous phantom
V = spm_vol(fullfile(sDtp,'roi_00044_segmented_s201007051500-0002-00001-000160-01.nii'));
Y = spm_read_vols(V);
Yphantom = Y;

job.layers_opt =0;
job.inclusion = 1;
% number of layers in the phantom
switch job.layers_opt
    case 0% homogeneous
        for ilayer=2:max(max(max(Y)))
            Yphantom(Yphantom==ilayer)=1;%MG
        end
    case 1% 2 layers
        Yphantom(Yphantom==2)=1;
        for ilayer=3:max(max(max(Y)))
            Yphantom(Yphantom==ilayer)=5;
        end
end

%number of inclusions
switch job.inclusion
    case 1% one in the gray matter
        inc_c1 = 6*ones(5,5,5);% en voxels
        Yphantom(24:28,20:24,30:34) = inc_c1;
    case 2% one in the grey matter and the other in the skin
        inc_c1 = 6*ones(6,5,4);% en voxels
        inc_c5 = 6*ones(4,5,4);
        Yphantom(50:55,46:50,39:42) = inc_c1;%%% sur l'image choisie
        Yphantom(39:42,30:34,30:33) = inc_c5;
end

Vphantom = struct('fname',fullfile(sDtp,['phantom' int2str(job.layers_opt) 'l_' int2str(job.inclusion) 'incl.nii']),...
    'dim',  V.dim,...
    'dt',   V.dt,...
    'pinfo',V.pinfo,...
    'mat',  V.mat);
Vphantom = spm_create_vol(Vphantom);
spm_write_vol(Vphantom, Yphantom);