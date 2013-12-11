function varargout = nirs_clusterdata(X, varargin)
% nirs_clusterdata Construct clusters from data.
%   [T Z] = nirs_clusterdata(X, CUTOFF) constructs clusters from data X.
%   X is a matrix of size M by N, treated as M observations of N
%   variables.  CUTOFF is a threshold for cutting the hierarchical
%   tree generated by LINKAGE into clusters. When 0 < CUTOFF < 2,
%   clusters are formed when inconsistent values are greater than
%   CUTOFF (see INCONSISTENT). When CUTOFF is an integer and CUTOFF >= 2,
%   then CUTOFF is considered as the maximum number of clusters to
%   keep in the hierarchical tree generated by LINKAGE. The output T is
%   a vector of size M containing a cluster number for each observation.
%   Z is a (m-1)-by-3 matrix, where m is the number of observations in the
%   original data. Columns 1 and 2 of Z contain cluster indices linked in pairs
%   to form a binary tree. The leaf nodes are numbered from 1 to m. Leaf nodes
%   are the singleton clusters from which all higher clusters are built. Each
%   newly-formed cluster, corresponding to row Z(I,:), is assigned the index
%   m+I. Z(I,1:2) contains the indices of the two component clusters that form
%   cluster m+I. There are m-1 higher clusters which correspond to the interior
%   nodes of the clustering tree. Z(I,3) contains the linkage distances between
%   the two clusters merged in row Z(I,:).
%
%   When 0 < CUTOFF < 2, T = CLUSTERDATA(X,CUTOFF) is equivalent to: 
%      Y = pdist(X, 'euclid');
%      Z = linkage(Y, 'single');
%      T = cluster(Z, 'cutoff', CUTOFF);
%
%   When CUTOFF is an integer >= 2, T = CLUSTERDATA(X,CUTOFF) is equivalent
%   to: 
%
%      Y = pdist(X,'euclid'); 
%      Z = linkage(Y,'single'); 
%      T = cluster(Z,'maxclust',CUTOFF)
%
%   [T Z] = CLUSTERDATA(X,'PARAM1',VAL1,'PARAM2',VAL2,...) provides more
%   control over the clustering through a set of parameter/value pairs.
%   Valid parameters are the following:
%
%      Parameter    Value
%      'distance'   Any of the distance metric names allowed by PDIST
%                   (the 'minkowski' option can be followed by the
%                   value of the exponent P).
%      'linkage'    Any of the linkage methods allowed by LINKAGE.
%      'cutoff'     Cutoff for inconsistent or distance measure.
%      'maxclust'   Maximum number of clusters to form.
%      'criterion'  Either 'inconsistent' or 'distance'.
%      'depth'      Depth for computing inconsistent values.
%   
%   See also PDIST, LINKAGE, INCONSISTENT, CLUSTER, DENDROGRAM, KMEANS.

%   Copyright 1993-2008 The MathWorks, Inc. 
%   $Revision: 1.13.4.4 $

if nargin < 2
  error('stats:clusterdata:TooFewInputs','Not enough input arguments.');
end

% Determine which args go where
pdistargs = {};
linkageargs = {};
donelist = [];
if ischar(varargin{1})
   j = 1;
   while(j<nargin-1)
      pname = varargin{j};
      if isequal(pname,'distance')
         pdistargs(end+1) = varargin(j+1);
         donelist(end+1:end+2) = j:j+1;
         if isequal(varargin{j+1},'minkowski')
            if (length(varargin)>=j+2) & isnumeric(varargin{j+2})
               pdistargs(end+1) = varargin(j+2);
               donelist(end+1) = j+2;
               j = j+1;
            end
         end
      elseif isequal(pname,'linkage')
         linkageargs(end+1) = varargin(j+1);
         donelist(end+1:end+2) = j:j+1;
      end
      j = j+2;
   end
   varargin(donelist) = [];
end
if length(varargin)==0
   error('stats:clusterdata:MissingParameter',...
         'Must specify a CUTOFF or MAXCLUST value.')
end

if isempty(pdistargs)
   pdistargs = {'euclid'};
end
if isempty(linkageargs)
   linkageargs = {'single'};
end

% Call calculation routines with their own args.  Note that the
% pdist calculation is done within the linkage function.
Z = linkage(X,linkageargs{1},pdistargs);
T = cluster(Z,varargin{:});

msg = nargoutchk(1, 2, nargout);
if isempty(msg)
    % nout = max(nargout, 1) - 1;
    varargout{1} = T;
    varargout{2} = Z;
else
    disp(msg)
end
    
% EOF
