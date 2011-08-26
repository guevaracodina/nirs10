function out = nirs_MCsegment_rebels(job)
%--> for rebel voxels
rebel_count = 0;
rebel_count_processed =1;
surnd = 3;

%loading shadow and Ysegmented
Ysegmented = job.Yseg;
Vb = spm_vol(job.Yb);Yb = spm_read_vols(Vb);
Vh = spm_vol(job.Yh);Yh = spm_read_vols(Vh);

V1 = spm_vol(job.Y1);Y1 = spm_read_vols(V1);
V2 = spm_vol(job.Y2);Y2 = spm_read_vols(V2);
V3 = spm_vol(job.Y3);Y3 = spm_read_vols(V3);
V4 = spm_vol(job.Y4);Y4 = spm_read_vols(V4);
V5 = spm_vol(job.Y5);Y5 = spm_read_vols(V5);

while rebel_count_processed>0
    surnd = surnd+2;
    
    for i=1+surnd:size(Ysegmented,1)-surnd
        % non attributed voxels are detected
        [posAlpha,posAleph]=find(squeeze(Yh(i,:,:))>=0 & squeeze(Ysegmented(i,:,:))==0);
        
        % on each surrounding, one compute the number of voxels of
        % each layer
        for j=1:size(posAlpha,1)
            nbrhood = Ysegmented(i-surnd:i+surnd,...
                max(1,posAlpha(j)-surnd):min(Vb.dim(2),posAlpha(j)+surnd),...
                max(1,posAleph(j)-surnd):min(Vb.dim(3),posAleph(j)+surnd));
            
            if posAlpha(j,1)<=surnd
                vox(2) = posAlpha(j,1);
            elseif posAlpha(j,1)>=V5.dim(2)-surnd+1
                vox(2) = V5.dim(2)-posAlpha(j,1);
            end
            if posAleph(j,1)<=surnd
                vox(3) = posAleph(j,1);
            elseif posAleph(j,1)>=V5.dim(3)-surnd+1
                vox(3) = V5.dim(3)-posAleph(j,1);
            end
            
            %%%% methode 1 % the most widespread value is the one given to the
            %%%% central voxel
%             rebel_value_count = zeros(6,1);
%             for i_value =1:6
%                 rebel_value_count(i_value) = size(find(nbrhood==(i_value-1)),1);
%             end
%             [dummy,layer] = max(rebel_value_count);
%             layer = layer-1;
            %%%% methode 2 % on interpole et on prend round des vqleurs
            %%%% trouvees
            
            %%%%---> on interpole en 3D
            % Create linear spaces for initial and final volume
            xi = 1:size(nbrhood,1);
            yi = 1:size(nbrhood,2);
            zi = 1:size(nbrhood,3);
            xf = 1:size(nbrhood,1);
            yf = 1:size(nbrhood,2);
            zf = 1:size(nbrhood,3);
            
            % Make grids, to use with 3-dimensional interpolation
            [xi1,yi1,zi1] = meshgrid(xi,yi,zi);
            [xf1,yf1,zf1] = meshgrid(xf,yf,zf);
            
            nbrhood = permute(nbrhood,[2,1,3]);
            nbrhoodI = interp3(xi1,yi1,zi1,nbrhood,xf1,yf1,zf1,'nearest');
            nbrhoodI = permute(nbrhoodI,[2,1,3]);
            disp(['i : ' int2str(i) ' et vox ; ' int2str(vox)])
            layer = nbrhoodI(surnd+1,vox(2),vox(3));
            % if layer ==0 hence air, then values from the
            % NewSegment layers are used
            if layer==0
                v_c1 = Y1(i,posAlpha(j),posAleph(j));
                v_c2 = Y2(i,posAlpha(j),posAleph(j));
                v_c3 = Y3(i,posAlpha(j),posAleph(j));
                v_c4 = Y4(i,posAlpha(j),posAleph(j));
                v_c5 = Y5(i,posAlpha(j),posAleph(j));
                [dummy,layer] = max([v_c1,v_c2,v_c3,v_c4,v_c5]);
                Ysegmented(i,posAlpha(j),posAleph(j)) = layer;
            else
                Ysegmented(i,posAlpha(j),posAleph(j)) = layer;
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
end
out = V.fname;
end