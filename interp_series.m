function Yi = interp_series(Y,kk,Q)
%Interpolated data series at point kk on map
%data series Y given as time x channels
%HbO, HbR and HbT are treated separately
if ~isempty(Q)
    if size(kk,2) == 1
        mode = 1;
    else
        mode = 0;
    end
    B = Q.B;
    Yi = zeros(size(Y,1),size(kk,1));
    for i=1:size(kk,1)
        if mode
            rmv = Q.rmask{1};
            cmv = Q.cmask{1};
            B2(:,1) = B(rmv(kk(i)), cmv(kk(i)),:);
        else
            B2(:,1) = B(kk(i,1),kk(i,2),:);
        end
        Yi(:,i) = Y*B2;
    end
    Yi=mean(Yi,2);
else
    if size(kk,2) == 1
        mode = 1;
    else
        mode = 0;
    end
    Yi=zeros(size(kk,1),1);
    for i=1:size(kk,1)
        if mode
            rmv = Q.rmask{1};
            cmv = Q.cmask{1};
            Yi(i) = Y(rmv(kk(i)), cmv(kk(i)));
        else
            Yi(i) = Y(kk(i,1),kk(i,2));
        end
    end
    Yi=mean(Yi);
end
end