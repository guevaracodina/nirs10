function out = nirs_criugm_trackpace(varargin)% Action,fft_slab,fft_freq_step

% Action = varargin{1};
% if size(varargin,2)<3
%     outbeattest = varargin{2};
% else
fft_slab = varargin{1};
fft_freq_step = varargin{2};
MinRate = varargin{3};
MaxRate = varargin{4};

%     if size(varargin,2)==4, energie = varargin{4};
%     elseif size(varargin,2)>4
%         energie = varargin{4};
%         heart_possiblepace = varargin{5};
%     end
% end

%eventuellement : energie des ondes de Mayer
%--- detection de la pointe associee aux ondes de Mayer :
% Détection de maximum :
% Mayer Waves
% c_mayer=1;
% for i=floor(0.1/fft_freq_step):floor(2/fft_freq_step)% on cherche entre 0.2 et 1Hz
%     if (fft_slab(i)>fft_slab(i+1) && fft_slab(i)>fft_slab(i-1))% alors on a un pas
%         mayer_possiblepace(c_mayer,1) = i;
%         c_mayer = c_mayer+1;
%     end
% end
%
% if ~isempty(mayer_possiblepace)
%     [mayer_amp mayer_pace] = max(fft_slab(mayer_possiblepace));
%     energie_mayer = sum(energie(mayer_possiblepace(mayer_pace)-6:mayer_possiblepace(mayer_pace)+6));
%     outbeattest.mayer.pace = mayer_possiblepace(mayer_pace)*fft_freq_step;
%     outbeattest.mayer.energie = energie_mayer;
% else
%     outbeattest.mayer.pace = 0;
%     outbeattest.mayer.energie =0;
% end


% mayer_possiblepace=zeros(0);
heart_possiblepace=zeros(0);

energie = fft_slab.^2;

c_heart=1;
for i=floor(MinRate/fft_freq_step):floor(MaxRate/fft_freq_step)% on cherche entre 0.8 et 3Hz - but was at 5Hz???
    if (fft_slab(i)>fft_slab(i+1) && fft_slab(i)>fft_slab(i-1))% alors on a un max
        heart_possiblepace(c_heart,1) = i;
        c_heart = c_heart+1;
    end
end


coeff =0.3;

if ~isempty(heart_possiblepace)
    [heart_amp heart_pace] = max(fft_slab(heart_possiblepace));
    % autour de ce max on cherche les deux minima :
    i_up = heart_possiblepace(heart_pace);
    while fft_slab(i_up)>fft_slab(i_up+1) && i_up<=1024
        i_up = i_up+1;
    end
    i_lo = heart_possiblepace(heart_pace);
    while fft_slab(i_lo)>fft_slab(i_lo-1) && i_lo>=0
        i_lo = i_lo-1;
    end
    % test sur l'energie :
    energie_heart = sum(energie(i_lo:i_up));
    if energie_heart> coeff*sum(energie(floor(MinRate/fft_freq_step):floor(MinRate/fft_freq_step)));
        outbeattest.heart.pace = heart_possiblepace(heart_pace)*fft_freq_step;
        outbeattest.heart.energie = energie_heart;

    else
        outbeattest.heart.pace = 0;
        outbeattest.heart.energie =0;
    end
else
    outbeattest.heart.pace = 0;
    outbeattest.heart.energie =0;
end
    
out = outbeattest;

% function out = CRiugm_trackpace(varargin)% Action,fft_slab,fft_freq_step
%
% Action = varargin{1};
% if size(varargin,2)<3
%     outbeattest = varargin{2};
% else
%     fft_slab = varargin{2};
%     fft_freq_step = varargin{3};
%     if size(varargin,2)==4, energie = varargin{4};
%     elseif size(varargin,2)>4
%         energie = varargin{4};
%         heart_possiblepace = varargin{5};
%     end
% end
%
% %eventuellement : energie des ondes de Mayer
% %--- detection de la pointe associee aux ondes de Mayer :
% %% Détection de maximum :
% % Mayer Waves
% % c_mayer=1;
% % for i=floor(0.1/fft_freq_step):floor(2/fft_freq_step)% on cherche entre 0.2 et 1Hz
% %     if (fft_slab(i)>fft_slab(i+1) && fft_slab(i)>fft_slab(i-1))% alors on a un pas
% %         mayer_possiblepace(c_mayer,1) = i;
% %         c_mayer = c_mayer+1;
% %     end
% % end
% %
% % if ~isempty(mayer_possiblepace)
% %     [mayer_amp mayer_pace] = max(fft_slab(mayer_possiblepace));
% %     energie_mayer = sum(energie(mayer_possiblepace(mayer_pace)-6:mayer_possiblepace(mayer_pace)+6));
% %     outbeattest.mayer.pace = mayer_possiblepace(mayer_pace)*fft_freq_step;
% %     outbeattest.mayer.energie = energie_mayer;
% % else
% %     outbeattest.mayer.pace = 0;
% %     outbeattest.mayer.energie =0;
% % end
%
%
% switch lower(Action)
%     case 'welcome'
%         mayer_possiblepace=zeros(0);
%         heart_possiblepace=zeros(0);
%
%         energie = fft_slab.^2;
%         CRiugm_trackpace('track_max',fft_slab,fft_freq_step,energie);
%
%     case 'track_max'
%         c_heart=1;
%         for i=floor(0.5/fft_freq_step):floor(5/fft_freq_step)% on cherche entre 0.8 et 3Hz
%             if (fft_slab(i)>fft_slab(i+1) && fft_slab(i)>fft_slab(i-1))% alors on a un max
%                 heart_possiblepace(c_heart,1) = i;
%                 c_heart = c_heart+1;
%             end
%         end
%         CRiugm_trackpace('comp_nrg',fft_slab,fft_freq_step,energie,heart_possiblepace);
%
%     case 'comp_nrg'
%         coeff =0.4;
%
%         if ~isempty(heart_possiblepace)
%             [heart_amp heart_pace] = max(fft_slab(heart_possiblepace));
%             % autour de ce max on cherche les deux minima :
%             i_up = heart_possiblepace(heart_pace);
%             while fft_slab(i_up)>fft_slab(i_up+1)
%                 i_up = i_up+1;
%             end
%             i_lo = heart_possiblepace(heart_pace);
%             while fft_slab(i_lo)>fft_slab(i_lo-1)
%                 i_lo = i_lo-1;
%             end
%             % test sur l'energie :
%             energie_heart = sum(energie(heart_possiblepace(heart_pace)-i_lo:heart_possiblepace(heart_pace)+i_up));
%             if energie_heart> coeff*sum(energie);
%                 outbeattest.heart.pace = heart_possiblepace(heart_pace)*fft_freq_step;
%                 outbeattest.heart.energie = energie_heart;
%                 CRiugm_trackpace('save',outbeattest);
%             else
%                 CRiugm_trackpace('track_max',fft_slab,fft_freq_step,energie);
%             end
%         else
%             outbeattest.heart.pace = 0;
%             outbeattest.heart.energie =0;
%         end
% end
% out = 1;