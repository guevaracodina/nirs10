function out = nirs_MCsegment_rebels2(job)
%--> for rebel voxels
rebel_count = 0;
surnd = 3; %why not 1? 3 allows more accurate attribution of layer

%loading shadow and Ysegmented
Ysegmented = job.Yseg;
Vb = spm_vol(job.Yb);%Yb = spm_read_vols(Vb);
Vh = spm_vol(job.Yh);Yh = spm_read_vols(Vh);

% if isfield(job,'Yb2')
%     Yb2 = job.Yb2;
% end
%
% Y1 = job.Y1;
% Y2 = job.Y2;
% Y3 = job.Y3;
% Y4 = job.Y4;
% Y5m = job.Y5m;
% Y6nus = job.Y6nus;

% while rebel_count_processed>0
spm_progress_bar('Init',size(Ysegmented,1),'Completing MC Segment','Progress');
for i=1+surnd:size(Ysegmented,1)-surnd
    spm_progress_bar('Set',i);
    % non attributed voxels are detected
    [posAlpha,posAleph]=find(squeeze(Yh(i,:,:))>0 & squeeze(Ysegmented(i,:,:))==0);
    
    if ~isempty(posAlpha)
        % on each surrounding, one compute the number of voxels of
        % each layer
        for j=1:size(posAlpha,1)
            
            %             if Y6nus(i,posAlpha(j),posAleph(j))==0
            %                 Ysegmented(i,posAlpha(j),posAleph(j))=6;
            %             else
            nbrhood = Ysegmented(i-surnd:i+surnd,...
                max(1,posAlpha(j)-surnd):min(Vb.dim(2),posAlpha(j)+surnd),...
                max(1,posAleph(j)-surnd):min(Vb.dim(3),posAleph(j)+surnd));
            
            %%%% methode 2 % on interpole et on prend round des vqleurs
            %%%% trouvees
            % % % % % % % % %                 if posAlpha(j,1)<=surnd
            % % % % % % % % %                     vox(2) = posAlpha(j,1);
            % % % % % % % % %                 elseif posAlpha(j,1)>=V5.dim(2)-surnd+1
            % % % % % % % % %                     vox(2) = V5.dim(2)-posAlpha(j,1);
            % % % % % % % % %                 else
            % % % % % % % % %                     vox(2) = surnd+1;
            % % % % % % % % %                 end
            % % % % % % % % %                 if posAleph(j,1)<=surnd
            % % % % % % % % %                     vox(3) = posAleph(j,1);
            % % % % % % % % %                 elseif posAleph(j,1)>=V5.dim(3)-surnd+1
            % % % % % % % % %                     vox(3) = V5.dim(3)-posAleph(j,1);
            % % % % % % % % %                 else
            % % % % % % % % %                     vox(3) = surnd+1;
            % % % % % % % % %                 end
            % % % % % % % % %                 %%%%---> on interpole en 3D
            % % % % % % % % %                 % Create linear spaces for initial and final volume
            % % % % % % % % %                 xi = 1:size(nbrhood,1);
            % % % % % % % % %                 yi = 1:size(nbrhood,2);
            % % % % % % % % %                 zi = 1:size(nbrhood,3);
            % % % % % % % % %                 xf = 1:size(nbrhood,1);
            % % % % % % % % %                 yf = 1:size(nbrhood,2);
            % % % % % % % % %                 zf = 1:size(nbrhood,3);
            % % % % % % % % %
            % % % % % % % % %                 % Make grids, to use with 3-dimensional interpolation
            % % % % % % % % %                 [xi1,yi1,zi1] = meshgrid(xi,yi,zi);
            % % % % % % % % %                 [xf1,yf1,zf1] = meshgrid(xf,yf,zf);
            % % % % % % % % %
            % % % % % % % % %                 nbrhood = permute(nbrhood,[2,1,3]);
            % % % % % % % % %                 nbrhoodI = interp3(xi1,yi1,zi1,nbrhood,xf1,yf1,zf1,'nearest');
            % % % % % % % % %                 nbrhoodI = permute(nbrhoodI,[2,1,3]);
            % % % % % % % % %
            % % % % % % % % %                 Ysegmented(i-surnd:i+surnd,...
            % % % % % % % % %                     max(1,posAlpha(j)-surnd):min(Vb.dim(2),posAlpha(j)+surnd),...
            % % % % % % % % %                     max(1,posAleph(j)-surnd):min(Vb.dim(3),posAleph(j)+surnd)) = nbrhoodI;
            %                 disp(['i : ' int2str(i) ' (vraie position : ' int2str(posAlpha(j,1)) ' ' int2str(posAleph(j,1)) ') et vox : ' int2str(vox)])
            % % %             % % % % % %             layer = nbrhoodI(surnd+1,vox(2),vox(3));
            
            %%%% methode 1 % the most widespread value is the one given to the
            %%%% central voxel
            rebel_value_count = zeros(6,1);
            for i_value =1:6
                rebel_value_count(i_value) = length(find(nbrhood==(i_value-1)));
            end
            [dummy,layer] = max(rebel_value_count);
            layer = layer-1;
            %                 if layer ==0 hence air, then values from the
            %                 NewSegment layers are used
            
            %%%% merging des deux chemins
            if layer==0
                nbrhood1 = Ysegmented(i-surnd:i+surnd,...
                    max(1,posAlpha(j)-5):min(Vb.dim(2),posAlpha(j)+5),...
                    max(1,posAleph(j)-5):min(Vb.dim(3),posAleph(j)+5));
                rebel_value_count = zeros(6,1);
                for i_value =1:6
                    rebel_value_count(i_value) = length(find(nbrhood1==(i_value-1)));
                end
                [dummy,layer] = max(rebel_value_count);
                layer = layer-1;
                Ysegmented(i,posAlpha(j),posAleph(j)) = layer;
                
                %                 v_c1 = Y1(i,posAlpha(j),posAleph(j));
                %                 v_c2 = Y2(i,posAlpha(j),posAleph(j));
                %                 v_c3 = Y3(i,posAlpha(j),posAleph(j));
                %
                %                 v_c4 = Y4(i,posAlpha(j),posAleph(j));
                %                 v_c5 = Y5(i,posAlpha(j),posAleph(j));
                %                 [confidence,layer2] = max([v_c1,v_c2,v_c3,v_c4,v_c5]);
                %                 if layer2~=0 %&& confidence>0.3
                %                     Ysegmented(i,posAlpha(j),posAleph(j)) = layer2;
                %                 else
                %                     Ysegmented(i,posAlpha(j),posAleph(j)) = 6;
                %                 end
            else
                Ysegmented(i,posAlpha(j),posAleph(j)) = layer;
            end
            %             end
            
        end
    end
    
    
    rebel_count = rebel_count +  size(posAlpha,1);
end

if surnd==3
    disp(['Remaining rebel voxels before processing : ',int2str(rebel_count)]);
end
rebel_count_processed = 0;
for i=1+surnd:size(Ysegmented,1)-surnd
    [posAlpha,dummy]=find(squeeze(Yh(i,:,:))>0 & squeeze(Ysegmented(i,:,:))==0);
    rebel_count_processed = rebel_count_processed +  size(posAlpha,1);
end
disp(['Remaining rebel voxels after processing : ',int2str(rebel_count_processed)]);
% % % % % % %     surnd = surnd+2;
% % % % % % % end

% calcul des covariances des voxels (voir Diffuse optical tomography with a priori anatomical information, Guven, 2005)
% on calcule l ecart type pour l'ense;ble des voxels
% appartenant a une couche
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% A CODER %

% save of nifti image :
Vsegmented = struct('fname',fullfile(job.dir,[job.output_prefix,'_segmented_',job.name,'.nii']),...
    'dim',  Vb.dim,...
    'dt',   Vb.dt,...
    'pinfo',Vb.pinfo,...
    'mat',  Vb.mat);
Vsegmented.pinfo(1,1)=1;% pinfo(1,1) doit etre egal a 1...
spm_write_vol(Vsegmented, Ysegmented);

out = Vsegmented.fname;
end