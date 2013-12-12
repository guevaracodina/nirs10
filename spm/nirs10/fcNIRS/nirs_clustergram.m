function varargout = nirs_clustergram(dataNIRS,varargin)
% nirs_clustergram creates a dendrogram and heat map on the same figure.
%
%   nirs_clustergram(dataNIRS) creates a dendrogram and heat map from DATA using
%   hierarchical clustering with Euclidean distance metric and average linkage
%   used to generate the hierarchical tree. The clustering is performed on the
%   rows of dataNIRS. The rows of dataNIRS are typically channels and the
%   columns are time points.  To cluster the columns instead of the rows,
%   transpose the data using the ' operator.

%% Using clusterdata (hierarchical clustering)
% Number of independent runs
nIter           = 1;
% Number of clusters to form
nClusters       = 12;
% Type of distance metric
distanceType    = 'correlation';
% Linkage method
linkageMethod   = 'weighted';
% FDR correction
doFDR           = true;

% nirs_clusterdata performs all of the necessary steps for you. You do not need
% to execute the pdist, linkage, or cluster functions separately.
T2 = zeros([size(dataNIRS,1) nIter]);
Z2 = zeros([size(dataNIRS,1)-1 3 nIter]);
groupsClusterData = zeros([nIter, nClusters]);

for i = 1:nIter,
    tic
    [T2(:,i), Z2(:,:,i)] = nirs_clusterdata(dataNIRS,'distance', distanceType,...
        'linkage', linkageMethod, 'maxclust', nClusters);
    for j = 1:nClusters,
        groupsClusterData(i,j) = numel(find(T2(:,i)==j));
    end
    fprintf('clusterData iter%d computed in: %s\n', i, datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));
end

%% clustergram (dendrogram + correlation matrix)
% position the figure
figure('units','normalized','Position',[0.5641    0.2407    0.3807    0.6426]);
mainPanel = axes('Position',[.25 .08 .69 .69]);
leftPanel = axes('Position',[.08 .08 .17 .69]);
topPanel =  axes('Position',[.25 .77 .69 .21]);

%% Left dendrogram
axes(leftPanel);
color = Z2(end-nClusters+2,3)-eps;
[H21, T31, PERM1] = dendrogram(Z2, 0, 'colorthreshold', color, 'orient','left');
% h = dendrogram(Z_samples,'orient','left');

%% Top dendrogram
axes(topPanel);
[H22, T32, PERM2] = dendrogram(Z2, 0, 'colorthreshold', color);

%% Correlation matrix
dataNIRSorder = dataNIRS';
dataNIRSorder = dataNIRSorder(:, PERM1);
tic; [corrMat p] = corrcoef(dataNIRSorder); toc
if doFDR
    fprintf('FDR correction done\n');
    p = nirs_fdr(p(:));
    p = reshape(p, size(corrMat));
%     p = nirs_fdr(p);
    % Make all non-significant elements equal to zero
    corrMat(p>0.05) = 0;
end
axes(mainPanel);
imagesc(corrMat,[-1 1]); colormap(nirs_get_colormap('redbluecmap')); 
set(mainPanel,'Xticklabel',[],'yticklabel',[]);
axis xy; 

%% align the x and y axis
set(leftPanel,'ylim',[1 size(corrMat,1)],'Visible','Off');
set(topPanel,'xlim',[1 size(corrMat,2)],'Visible','Off');
axes(mainPanel);axis([1 size(corrMat,2) 1 size(corrMat,1)]);

%% Add annotations
axes(mainPanel); 
set(gca, 'XTickLabel', cellfun(@(x) num2str(x), num2cell(PERM1), 'UniformOutput',false))
xlabel('NIRS channels','Fontsize',14); 
% Color bar
colorbar('Location','northoutside','Position',[ 0.0584    0.8761    0.3082    0.0238]); 
set(leftPanel,'yaxislocation','left');
% set(get(leftPanel,'ylabel'),'string','Samples','Fontsize',14);
set(findall(leftPanel, 'type', 'text'), 'visible', 'on');
set(gcf,'color',[1 1 1],'paperpositionmode','auto');

% EOF
