function varargout = spatial_LPF(mode,K,data2D)
%Gaussian spatial LPF - based on NIRS_SPM function spm_filter_HPF_LPF_WMDL
%Careful for boundary effects with this filter! near the edges, it will reduce the amplitude by
%as much as 50% if radius is 3, and by much more for larger radius
switch mode
    case 'set'        
        FWHM = 2*K.radius;
        sigma   = FWHM;
        h       = round(4*sigma);
        h       = exp(-(-h:h).^2/(2*sigma^2));
        n       = length(h);
        d       = (1:n) - (n + 1)/2;
        if      n == 1, h = 1; end
        
        k = K.k1;
        L = spdiags(ones(k,1)*h, d, k,k);
        K.K1 = spdiags(1./sum(L')',0,k,k)*L;
        k = K.k2;
        L = spdiags(ones(k,1)*h, d, k,k);
        K.K2 = spdiags(1./sum(L')',0,k,k)*L;
        K.Ks1 = K.K1*K.K1';
        K.Ks2 = K.K2*K.K2';
        varargout{1} = K;
    case 'lpf'
        ndata2D = zeros(size(data2D));
%         for i=1:size(data2D,1)
%             ndata2D(i,:,:) = K.K1 * squeeze(data2D(i,:,:)) * K.K2';
%         end
        for i=1:size(data2D,1)
            ndata2D(i,:,:) = K.Ks1 * squeeze(data2D(i,:,:)) * K.Ks2';
        end

        varargout{1} = ndata2D;
    otherwise
end
end