function colormapOut = nirs_get_colormap(map)
% Creates colormaps that are adequate to display images in SS-OCT system. Also
% creates maps adequate for photoacoustic tomography.
% SYNTAX:
% colormapOut = nirs_get_colormap(map)
% INPUTS:
% map           String that describes the colormap to retrieve:
%               'octgold'
%               'fdrainbow'
%               'tdrainbow'
%               'wob'
%               'bow'
%               'flow'
%               'rwbdoppler'
%               'bwrdoppler'
%               'robdoppler'
%               'bordoppler'
%               'redmap'
%               'so2'
%               'redbluecmap'
% OUTPUTS:
% colormapOut   3 columns matrix, which values are in the range from 0 to 1.
%_______________________________________________________________________________
% Copyright (C) 2013 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

ColorMapSize = 128;

switch lower(map)
    case 'octgold'
        % OCT Gold: typical high contrast, high dynamic OCT map
        x = [0 7 60 108 157 190 224 253 254 255] + 1;
        rgb = [ ...
            0.00392157  0.00392157  0.0156863;
            0.0352941   0           0;
            0.376471    0.105882    0;
            0.678431    0.337255    0;
            0.992157    0.701961    0.203922;
            0.992157    0.878431    0.341176;
            0.996078    0.996078    0.482353;
            0.94902     0.996078    0.603922;
            0.905882    0.952941    0.607843;
            0.905882    0.952941    0.607843];
    case 'fdrainbow'
        % Fourier-Domain Rainbow: multi-colour, high dynamic range colour map for
        % novel OCT systems
        x = [0 22 53 68 106 141 168 194 245 255] + 1;
        rgb = [ ...
            0.00392157  0.00392157  0.0156863;
            0.0313725   0.0313725   0.117647;
            0.00392157  0.113725    0.541176;
            0.203922    0.486275    0.207843;
            0.294118    0.780392    0.145098;
            0.847059    0.752941    0.435294;
            1           0.966912    0.0735296;
            0.992157    0.792157    0.411765;
            1           0.3         0.3;
            1           0.984314    0.984314];
    case 'tdrainbow'
        % Time-Domain Rainbow: mimicking the old time domain high dynamic range
        % colour map found in the commercial devices
        x = [32 64 78 138 221 255] + 1;
        rgb = [ ...
            0.00392157  0.00392157  0.0156863;
            0.3         0.3         1;
            0.176471    0.760784    0.254902;
            1           1           0.301961;
            1           0.3         0.3;
            1           0.984314    0.984314];
    case 'wob'
        % grey-scale with transparency for frequency domain OCT data
        x = [0 40 72 136 255] + 1;
        rgb = [...
            0.00392157  0.00392157  0.0156863;
            0           0           0;
            0.301961    0.301961    0.301961;
            0.733333    0.733333    0.733333;
            1           1           1];
    case 'bow'
        % grey-scale with transparency for frequency domain OCT data (not very
        % commonly used by the community, but good for faint signals)
        x = [0 35 96 148 255] + 1;
        rgb = [...
            1           1           1;
            0.890196    0.890196    0.890196;
            0.584314    0.584314    0.584314;
            0.282353    0.282353    0.282353;
            0           0           0];
    case 'flow'
        % High flow velocity contrast 20kHz
        x = [0 69 128 184 255] + 1;
        rgb = [...
            0.133333    0.188235    1;
            0.301961    0.760784    1;
            0.286275    0.913725    0.27451;
            0.972549    1           0.0862745;
            1           0.1         0.1];
    case 'rwbdoppler'
        % Red on blue, with white background for Doppler imaging
        % Also for SO2 contrast in photoacosutics
        minColor    = [0 0 1]; % blue
        medianColor = [1 1 1]; % white   
        maxColor    = [1 0 0]; % red      
       
        int1 = zeros(ColorMapSize,3); 
        int2 = zeros(ColorMapSize,3);
        for k=1:3
            int1(:,k) = linspace(minColor(k), medianColor(k), ColorMapSize);
            int2(:,k) = linspace(medianColor(k), maxColor(k), ColorMapSize);
        end
        colormapOut = [int1(1:end-1,:); int2];
        return
    case 'bwrdoppler'
        % Blue on red, with white background for Doppler imaging
        minColor    = [1 0 0]; % red
        medianColor = [1 1 1]; % white
        maxColor    = [0 0 1]; % blue
        
        int1 = zeros(ColorMapSize,3);
        int2 = zeros(ColorMapSize,3);
        for k=1:3
            int1(:,k) = linspace(minColor(k), medianColor(k), ColorMapSize);
            int2(:,k) = linspace(medianColor(k), maxColor(k), ColorMapSize);
        end
        colormapOut = [int1(1:end-1,:); int2];
        return
    case 'robdoppler'
        % Red on blue, with black background for Doppler imaging
        minColor    = [0 0 1]; % blue
        medianColor = [0 0 0]; % white   
        maxColor    = [1 0 0]; % red      

        int1 = zeros(ColorMapSize,3); 
        int2 = zeros(ColorMapSize,3);
        for k=1:3
            int1(:,k) = linspace(minColor(k), medianColor(k), ColorMapSize);
            int2(:,k) = linspace(medianColor(k), maxColor(k), ColorMapSize);
        end
        colormapOut = [int1(1:end-1,:); int2];
        return
    case 'bordoppler'
        % Blue on red, with black background for Doppler imaging
        minColor    = [1 0 0]; % red
        medianColor = [0 0 0]; % black
        maxColor    = [0 0 1]; % blue
        
        int1 = zeros(ColorMapSize,3);
        int2 = zeros(ColorMapSize,3);
        for k=1:3
            int1(:,k) = linspace(minColor(k), medianColor(k), ColorMapSize);
            int2(:,k) = linspace(medianColor(k), maxColor(k), ColorMapSize);
        end
        colormapOut = [int1(1:end-1,:); int2];
        return
    case 'so2'
        % Red on blue, sO2 PAT imaging
        minColor    = [0 0 1]; % blue
        medianColor = [0.5 0 0.5]; % ??
        maxColor    = [1 0 0]; % red
        
        int1 = zeros(ColorMapSize,3);
        int2 = zeros(ColorMapSize,3);
        for k=1:3
            int1(:,k) = linspace(minColor(k), medianColor(k), ColorMapSize);
            int2(:,k) = linspace(medianColor(k), maxColor(k), ColorMapSize);
        end
        colormapOut = [int1(1:end-1,:); int2];
        return
    case 'redmap'
        % Red map for HbT contrast in PAT, from VisualSonics
        colormapOut = [1 0 0];
        return
    case 'redbluecmap'
        % red represents values above the mean, white represents the mean, and
        % blue represents values below the mean
        colormapOut = redbluecmap;
        return
    otherwise
        % Inverted linear gray colormap
        colormapOut = flipud(colormap(gray(255)));
        return
end

%% Calculate colormap
% ----------------------- Piecewise linear interpolation -----------------------
nSegments           = numel(x) - 1;
samplesPerSegment   = diff(x);
colormapOut         = zeros([sum(samplesPerSegment) 3]);

for iSegments = 1:nSegments,
    for iColors = 1:3,
    colormapOut(x(iSegments):x(iSegments+1),iColors) = linspace(rgb(iSegments,iColors),...
        rgb(iSegments+1,iColors),...
        samplesPerSegment(iSegments)+1);
    end
end

% --------------- uint8 array to export colormaps to LabView -------------------
% colormapOut = uint8(round(colormapOut*255));
% figure; plot(colormapOut)

% ---------------------------- Spline interpolation ----------------------------
% colormapOut(:,1) = spline(x, rgb(:,1), (x(1):x(end))');
% colormapOut(:,2) = spline(x, rgb(:,2), (x(1):x(end))');
% colormapOut(:,3) = spline(x, rgb(:,3), (x(1):x(end))');
% colormapOut(colormapOut<0) = 0;         % Get rid of negative numbers
% colormapOut(colormapOut>1) = 1;         % Truncate to 1
% figure; plot(colormapOut)

%% OCT Gold Alpha channel
% x                   = [0; 22; 70; 123; 251] + 1;
% alphaValues         = [0; 0; 0.0543478; 0.717391; 1];
% nPoints             = numel(x);
% nSegments           = nPoints - 1;
% samplesPerSegment   = diff(x);
% OCTgoldAlphaMap     = zeros([sum(samplesPerSegment) 1]);
% 
% for iSegments = 1:nSegments,
%     OCTgoldAlphaMap(x(iSegments):x(iSegments+1)) = linspace(alphaValues(iSegments),...
%         alphaValues(iSegments+1),...
%         samplesPerSegment(iSegments)+1);
% end

function c = redbluecmap(m,varargin)
%REDBLUECMAP creates a red and blue colormap.
%
%   REDBLUECMAP(M) returns an M-by-3 matrix containing a red and blue
%   diverging color palette. M is the number of different colors in the
%   colormap with a minimun of 3 and a maximun of 11. Low values are dark
%   blue, values in the center of the map are white, and high values are
%   dark red. If M is empty, a default value of 11 will be used.
%
%   Example:
% 
%       % Reset the colormap of the current figure, type
%             colormap(redbluecmap)
%
%   See also CLUSTERGRAM, COLORMAP, COLORMAPEDITOR, REDGREENCMAP.

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2.4.1 $  $Date: 2008/01/23 21:09:34 $

% Reference: 
% http://colorbrewer.org.

%== Setting default
if nargin < 1 || isempty(m) || ~isnumeric(m)
    m = 11;
end

if ~isscalar(m)
    m = m(:);
end

m = max(abs(fix(m)), 3);
m = min(m, 11);

switch (m)
    case 3
        c = [239	138     98;
             247	247     247;
             103	169     207];
    case 4
        c = [202	0       32;
             244	165     130;
             146	197     222;
             5      113     176];
    case 5
        c = [202	0       32;
             244	165     130;
             247	247     247;
             146	197     222;
             5      113     176];
    case 6
        c = [178	24      43;
             239	138     98;
             253	219     199;
             209	229     240;
             103	169     207;
             33     102     172];
    case 7
        c = [178	24      43;
             239	138     98;
             253	219     199;
            247     247     247;
            209     229     240;
            103     169     207;
            33      102     172];
    case 8
        c = [178	24      43;
             214	96      77;
             244	165     130;
             253	219     199;
             209	229     240;
             146	197     222;
             67     147     195;
             33     102     172];
    case 9
        c = [178	24      43;
             214	96      77;
             244	165     130;
             253	219     199;
             247	247     247;
             209	229     240;
             146	197     222;
             67     147     195;
             33     102     172];
    case 10
        c = [103	0       31;
            178     24      43;
            214     96      77;
            244     165     130;
            253     219     199;
            209     229     240;
            146     197     222;
            67	    147     195;
            33      102     172;
            5       48      97];
    case 11
        c = [103    0       31;
            178     24      43;
            214     96      77;
            244     165     130;
            253     219     199;
            247     247     247;
            209     229     240;
            146     197     222;
            67      147     195;
            33      102     172;
            5       48      97];
end
c = flipud(c/255);
x = round(linspace(1,256,11));
% ----------------------- Piecewise linear interpolation -----------------------
nSegments           = numel(x) - 1;
samplesPerSegment   = diff(x);
colormapOut         = zeros([sum(samplesPerSegment) 3]);

for iSegments = 1:nSegments,
    for iColors = 1:3,
    colormapOut(x(iSegments):x(iSegments+1),iColors) = linspace(c(iSegments,iColors),...
        c(iSegments+1,iColors),...
        samplesPerSegment(iSegments)+1);
    end
end
c = colormapOut;
% ==============================================================================
% [EOF]
