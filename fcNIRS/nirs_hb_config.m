function [inv_exs2, EPF2] = nirs_hb_config(PVF, DPF, nChannels, Cgp, Cwl, wl, nIter)
% Pre-computation of extinction coefficients and effective path length
% SYNTAX
% INPUTS
% OUTPUTS
%_______________________________________________________________________________
% Copyright (C) 2014 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Partial volume correction factor
% PVF = [50 50];
% Differential pathlength factor
% From Duncan et al., Phys. Med. Biol. 1995 40(2):295-304.
% DPF = [6.51 5.86]; % For 690nm and 832nm
% NOTE: Portable NIRS-EEG wavelengths are 735nm and 850nm
DPF = repmat(DPF',[1 nChannels]); % nLambda x nChannels
PVF = repmat(PVF',[1 nChannels]); % nLambda x nChannels
EPF = zeros(1, nChannels); %EPF = L * PPF = L * DPF/PVF at each wavelength
% Source-detector distances
% Cgp = NIRS.Cf.H.C.gp;
% Channels wavelength
% Cwl = NIRS.Cf.H.C.wl;
% Device wavelength
% wl = NIRS.Cf.dev.wl;
% for iChannels = 1:NIRS.Cf.H.C.N
for iChannels = 1:nChannels
    EPF(1,iChannels) = Cgp(1,iChannels)*DPF(Cwl(1,iChannels),iChannels)./PVF(Cwl(1,iChannels),iChannels);
end
% Effective path length
EPF2 = ones(nIter,1)*EPF; %PP used to be EPF2 = EPF *ones(1,size(d,2));
% exs(:,1): HbO for each wavelength; exs(:,2): HbR for each wavelength Alexis's
% choice of extinction coefficients corresponds to case 1 in Homer, which
% appears to be their preferred choice too
[exs, extcoeff_ref] = GetExtinctions(wl,1);
inv_exs = pinv(exs(:,1:2)); % inv_exs has size 2 x #wl
inv_exs2 = kron(inv_exs,eye(nChannels/size(wl,2))); % size #pairs*2 x #pairs*#wl
% (number of channels nChannels = number of pairs x number of wavelengths)


% % Partial volume correction factor
% PVF = [50 50];
% % Differential pathlength factor
% % From Duncan et al., Phys. Med. Biol. 1995 40(2):295-304.
% DPF = [6.51 5.86]; % For 690nm and 832nm
% % NOTE: Portable NIRS-EEG wavelengths are 735nm and 850nm
% DPF = repmat(DPF',[1 nChannels]); % nLambda x nChannels
% PVF = repmat(PVF',[1 nChannels]); % nLambda x nChannels
% EPF = zeros(1, nChannels); %EPF = L * PPF = L * DPF/PVF at each wavelength
% % Source-detector distances
% Cgp = NIRS.Cf.H.C.gp;
% % Channels wavelength
% Cwl = NIRS.Cf.H.C.wl;
% % Device wavelength
% wl = NIRS.Cf.dev.wl;
% % for iChannels = 1:nChannels
% for iChannels = 1:NIRS.Cf.H.C.N
%     EPF(1,iChannels) = Cgp(1,iChannels)*DPF(Cwl(1,iChannels),iChannels)./PVF(Cwl(1,iChannels),iChannels);
% end
% % Effective path length
% EPF2 = ones(nIter,1)*EPF; %PP used to be EPF2 = EPF *ones(1,size(d,2));
% % exs(:,1): HbO for each wavelength; exs(:,2): HbR for each wavelength Alexis's
% % choice of extinction coefficients corresponds to case 1 in Homer, which
% % appears to be their preferred choice too
% [exs, extcoeff_ref] = GetExtinctions(wl,1);
% inv_exs = pinv(exs(:,1:2)); % inv_exs has size 2 x #wl
% inv_exs2 = kron(inv_exs,eye(nChannels/size(wl,2))); % size #pairs*2 x #pairs*#wl
% % (number of channels nChannels = number of pairs x number of wavelengths)

% EOF
