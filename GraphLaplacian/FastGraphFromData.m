function [vGraph,vNNInfo,vOpts]= FastGraphFromData( cX, cOpts, cF, cVerbose )

%
% function [vGraph,vNNInfo,vOpts]= FastGraphFromData( cX, cOpts, cF, cVerbose )
%
% IN:
%   cOpts         : structure containing:
%                     [Type]    : 'nn','gauss', 'harmonic', 'selftuning'. Default : 'gauss'.
%                     [NNsymm]  : how to symmetrize the weight matrix (empty otherwise)
%                               'ave' for 1/2(W+W*)
%                               'leftadj' for W*W
%                               'rightadj' for WW*
%                               Default: does not symmetrize.
%                      [NNsymm01] : Boolean: true means 0 or 1 weights after symmetrization. NNsymm needs to be set. Default: false.
%                     if Type == 'nn':
%                               [kNN]           : number of nearest neighbors if cType=='nn'. Default : 10
%                               [Delta]         : width of Gaussian or harmonic weight (depending on WeightType)
%                                                   If Delta == Inf, then uses 1 for all edge weights.
%                                                   If Delta == 0, then uses location-depedent Gaussian width,
%                                                       comparable to the distance of the farther kNN-th nearest neighbor                           
%                                                   Default: Inf.
%                               [WeightType]    : 'gauss' or 'harmonic'. Default: 'gauss'
%                     if Type == 'gauss' or 'harmonic'
%                               [Delta]         : width of Gaussian. Default: estimate with FastLaplacianEstimateGoodDelta
%                               [Precision]     : determines the range_search radius when cType=='gauss'. Default: 1e-3.
%                               [Radius]        : determines the range_search radius when cType=='gauss'. If specified it overrides Precision.
%                     if Type =='selftuning'
%                               [kNN]           : number of nearest neighbors if cType=='nn'. Default : 10.
%                               [kNNdelta]      : the kNNsigma-th nearest neighbor is used to set the scale delta. Default: 7.
%                               [alpha]         : factor in front of the scaling \sigma. See self-tuning spectral clustering paper. Default: 1.
%                     [Epsilon]   : no diffusion outside any local metric ball of radius Epsilon (used only by flexible metrics)
%                     [kNNCoarse] : do kNNCoarse fast nearest neighbor with constraints on certain coordinates. Default: [].
%                     [NNCoordinates] : used only if kNNCoarse is specified, is the list of coordinates to which constrain a first kNNCoarse nn search.
%                     [GetWeightsFcn] : pointer to function to compute weights, as in GetWeightsFcn_Template.m
%                     [GetWeightsFcnParam] : an extra parameter if a custom weight function is provided
%                     [UseFastNNSearch] : whether to use a fast nearest neighbor search or not. If not specified, decides which is likely to perform best based on number of points and dimensionality.
%                     [SparsifyW] : structure containing the following fields:
%                               Type : 'absolute','relative','quantile'
%                               Abs  : true or false depending on whether the size or distribution should be of the entries or of the absolute value of the entries of W
%                               Value : a real number describing the upper bound on the entries of W, or the upper bound on the relative size of the entries of W,
%                                       or the percentage of entries of W to be kept
%   cF            : the function whose gradient determines the weights.
%
% OUT:
%   vGraph      : sparse graph obtained from the data. Contains the fields D and W. W is the adjancency matrix (with weights), D is just the vector of row sums of W.
%   vNNInfo     : structure containing the following.
%                   Atria  : tree for fast nn and range searches
%                 If cType == 'nn':
%                   NNIdxs  : M by cOpts.kNN matrix whose (i,j) entry is the j-th nn to the i-th point
%                   NNDist  : M by cOpts.kNN matrix whose (i,j) entry is the distance of the j-th nn to the i-th point (these are sorted in increasing order
%                 If cType == 'gauss':
%                   NNCounts : M column vector whose i-th entry is the number of \epsilon neighbors of the i-th point
%                   NNIdxs   : cell array whose i-th cell contains the indices of the \epsilon neighbors of the i-th point
%                   NNDist   :  cell array whose i-th cell contains the distances of the \epsilon neighbors of the i-th point from the i-th point, NOT sorted
%   vOpts       : version of cOpts actually used.
%
% USES:
%   nn_prepare, nn_search, range_search, FastLaplacianEstimateGoodDelta, bruteforce_nn_search
%
%
% SC:
%   MM      03/13/04    [it's pretty good now, fast and includes range_searches]
%   MM      05/31/04    [added generic weight function also for nn, reorganized code, added epsilon for calls to flexible metrics, updated doc]
%   MM      07/06/04    [added some parameter checking]
%   MM      10/04/05    [added self-tuning]
%
% (c) Copyright Plain Sight Systems Inc., 2004
%

vOpts = cOpts;

% Default arguments
if nargin<2,
    cOpts = [];
end;
if nargin<4,
    cVerbose = false;
end;
if ~isfield(cOpts,'Type'),
    cOpts.Type = 'gauss';
end;
if strcmpi(cOpts.Type,'gauss') | strcmpi(cOpts.Type,'harmonic'),
    % Gaussian graph
    if ~isfield(cOpts,'Delta'),     
        cOpts.Delta = FastLaplacianEstimateGoodDelta( cX );
    end;
    if ~isfield(cOpts,'Precision'),
        cOpts.Precision = 1e-3;
    end;
elseif strcmpi(cOpts.Type,'nn'),
    % k-NN graph
    if ~isfield(cOpts,'Delta'),
        cOpts.Delta = Inf;
    end;
    if ~isfield(cOpts,'kNN'),
        cOpts.kNN = 10;
    end;
    % If the number of points is smaller than the number of neighbors, fix this and warn
    if cOpts.kNN > size(cX,1),
        cOpts.kNN = size(cX,1);
        warning('FastGraphFromData:Number of nearest neighbors set to size of set!');
    end;
    if ~isfield(cOpts,'WeightType'),
        cOpts.WeightType = 'gauss';
    end;
else
    % Self-tuning graph
    cOpts.Type = 'selftuning';
    if ~isfield(cOpts,'kNN'),
        cOpts.kNN = 10;
    end;
    if ~isfield(cOpts,'kNNdelta'),
        cOpts.kNNdelta = min([7,cOpts.kNN]);
    end;
    if ~isfield(cOpts,'alpha'),
        cOpts.alpha = 1;
    end;
end;

if (~isfield(cOpts,'NNsymm01')) | (isempty(cOpts.NNsymm01)),
    cOpts.NNsymm01 = false;
end;


% Number of points and dimensionality
[lNumberOfPoints lDim] = size(cX);

% Check if a weight function is specified
if isfield(cOpts,'GetWeightsFcn'),
    lWeightsFcnSpecified = true;
    if ~isfield(cOpts,'Epsilon')
        cOpts.Epsilon = 0;
    end;
    if ~isfield(cOpts,'GetWeightsFcnParam'),
        cOpts.GetWeightsFcnParam = [];
    end;
else
    lWeightsFcnSpecified = false;
end;

% Check if use of fast nn-search is specified
if ~isfield(cOpts,'UseFastNNSearch'),
    % Decide on whether to use it or not...
    if ( lDim>lNumberOfPoints ) | ( lDim > 300 )
        cOpts.UseFastNNSearch = false;
    else
        cOpts.UseFastNNSearch = true;
    end;
end;

% Pre-process points
cX = full(cX);                                  % TODO: modify code so that it accepts sparse matrices...could be a pain

if cOpts.UseFastNNSearch,
    if cVerbose, fprintf('Preparing for fast nearest neighbors...'); end;
    vNNInfo.Atria = nn_prepare( cX );
    if cVerbose, fprintf('done.\n'); end;
end;
% Check if vector field F is specified
if (nargin<4) | (isempty(cF))
    lFSpecified = false;
else
    lFSpecified = true;
    lDeltaSquared = cOpts.Delta^2;
end;

if (strcmpi(cOpts.Type,'nn') == 1) | (strcmpi(cOpts.Type,'selftuning')),
    if (isfield(cOpts,'kNNCoarse')==1) & (length(cOpts.kNNCoarse)>0),
        if cVerbose, fprintf('Computing constrained nearest neighbors...'); end;
        % Do fast nearest neighbors constrained on certain coordinates
        [vNNInfo.NNIdxs,vNNInfo.NNDist,vNNInfo.NNDist1,vNNInfo.NNDist2]=FastGraphFromDataConditionalNNSearch(cX,cOpts);%FGFDCFNN(cX,cOpts);
    else
        % Perform the nn search
        if cOpts.UseFastNNSearch,
            if cVerbose, fprintf('Computing fast nearest neighbors...'); end;
            [vNNInfo.NNIdxs vNNInfo.NNDist] = nn_search( cX, vNNInfo.Atria, cX, cOpts.kNN );
        else
            if cVerbose, fprintf('Computing nearest neighbors...'); end;
            [vNNInfo.NNIdxs vNNInfo.NNDist] = bruteforce_nn_search( cX, cX, cOpts.kNN );
        end;
    end;
    % If self-tuning weights, compute the \sigma_i's
    if strcmpi(cOpts.Type,'selftuning')==1,
        if (isfield(cOpts,'kNNCoarse')==1) & (length(cOpts.kNNCoarse)>0),
            lSigma1 = vNNInfo.NNDist1(:,cOpts.kNNdelta)';
            lSigma2 = vNNInfo.NNDist2(:,cOpts.kNNdelta)';
        else
            lSigma = vNNInfo.NNDist(:,cOpts.kNNdelta)';
        end;
        lSigma(find(lSigma==0))=1;          % In order to avoid 0/0 below...
    end;
    if cVerbose, fprintf('done.\n'); end;
    % Total number of entries in the W matrix
    lTotalSize = cOpts.kNN*lNumberOfPoints;
    % We need to create the sparse W matrix all at once (Matlab is VERY slow in adding one element to a sparse matrix)
    % Allocate memory for needed indices
    lIdx = 1;lIdxsI=[];lIdxsJ=[];lIdxsI(lTotalSize) = int16(0);lIdxsJ(lTotalSize) = int16(0);lEntries = zeros(1,lTotalSize);
    % Re-organize the stuff returned from the nn search.
    if cVerbose, fprintf('Building matrix of weights...'); end;
    for lk = 1:lNumberOfPoints
        % Build the sparse row and column index sets
        lIdxsI(lIdx:lIdx+cOpts.kNN-1) = lk;
        lIdxsJ(lIdx:lIdx+cOpts.kNN-1) = vNNInfo.NNIdxs(lk,:);
        if ~lFSpecified
            % No vector field specified
            if lWeightsFcnSpecified
                % A weight function is specified
                if isfield(cOpts,'Delta'),
                    lExtraParam.Delta = cOpts.Delta;
                end;
                lExtraParam.Distances = vNNInfo.NNDist(lk,:);
                if isfield(cOpts,'Epsilon'),
                    lExtraParam.Epsilon = cOpts.Epsilon;
                end;
                lExtraParam.Custom = cOpts.GetWeightsFcnParam;
                lExtraParam.CenterIdx = lk;
                % Call the weight function
                try
                    lW = feval( cOpts.GetWeightsFcn,cX,vNNInfo.NNIdxs(lk,:),cX(lk,:), lExtraParam );
                    lEntries(lIdx:lIdx+cOpts.kNN-1) = lW;
                catch
                    warning('FastGraphFromData:','Could not call the GetWeightsFcn, using default Gaussian weights');
                    % No weight function is specified: use standard Gaussian weight
                    if cOpts.Delta == 0,
                        lDelta = vNNInfo.NNDist(lk,cOpts.kNN)/sqrt(-log(0.05));
                        if lDelta == 0,
                            lDelta = 0.0001;
                        end;
                    else
                        lDelta = cOpts.Delta;
                    end;
                    if strcmpi(cOpts.WeightType,'gauss'),
                        lEntries(lIdx:lIdx+cOpts.kNN-1) = exp(-(vNNInfo.NNDist(lk,:)/lDelta).^2);
                    elseif strcmpi(cOpts.WeightType,'harmonic'),
                        lEntries(lIdx:lIdx+cOpts.kNN-1) = lDelta./(lDelta+vNNInfo.NNDist(lk,:));
                    end;
                end;
            else
                % No weight function is specified
                if (strcmpi(cOpts.Type,'nn') == 1),
                    %use standard Gaussian weight
                    if cOpts.Delta == 0,
                        lDelta = vNNInfo.NNDist(lk,cOpts.kNN)/sqrt(-log(0.05));
                        if lDelta == 0,
                            lDelta = 0.0001;
                        end;
                    else
                        lDelta = cOpts.Delta;
                    end;
                    if strcmpi(cOpts.WeightType,'gauss'),
                        lEntries(lIdx:lIdx+cOpts.kNN-1) = exp(-(vNNInfo.NNDist(lk,:)/lDelta).^2);
                    elseif strcmpi(cOpts.WeightType,'harmonic'),
                        lEntries(lIdx:lIdx+cOpts.kNN-1) = lDelta./(lDelta+vNNInfo.NNDist(lk,:));
                    end;
                else
                    % Use self-tuning weight
                    if (isfield(cOpts,'kNNCoarse')==1) & (length(cOpts.kNNCoarse)>0),
                        lEntries(lIdx:lIdx+cOpts.kNN-1) = exp(-(vNNInfo.NNDist1(lk,:).*vNNInfo.NNDist1(lk,:))./(cOpts.alpha*lSigma1(lk)*lSigma1(vNNInfo.NNIdxs(lk,:)))...
                                                              -(vNNInfo.NNDist2(lk,:).*vNNInfo.NNDist2(lk,:))./(cOpts.alpha*lSigma2(lk)*lSigma2(vNNInfo.NNIdxs(lk,:))));
                    else                        
                        lEntries(lIdx:lIdx+cOpts.kNN-1) = exp(-(vNNInfo.NNDist(lk,:).^2)./(cOpts.alpha*lSigma(lk)*lSigma(vNNInfo.NNIdxs(lk,:))));
                    end;
                end;
            end;
        else
            % A vector field is specified: flow along it
            lEntries(lIdx:lIdx+cOpts.kNN-1) = exp(-((vNNInfo.NNDist(lk,:)/cOpts.Delta).^2+abs(cF(vNNInfo.NNIdxs(lk,:))-cF(lk)))/lDeltaSquared);
        end;
        lIdx = lIdx + cOpts.kNN;
    end;
    if cVerbose, fprintf('done.\n'); end;

elseif strcmpi(cOpts.Type,'gauss') | strcmpi(cOpts.Type,'harmonic'),
    % Compute the range for requested precision
    if isfield(cOpts,'Radius'),
        if cOpts.Radius>0,
            lRange = cOpts.Radius;
        end;
    else
        if isfield(cOpts,'Precision')
            if cOpts.Precision <= 0, cOpts.Precision = 1e-4;end;
        else
            cOpts.Precision = 1e-4;
        end;
        % TODO: this is wrong, or at least inefficient, if cOpts.Delta=Inf...
        if strcmpi(cOpts.Type,'gauss'),
            lRange = cOpts.Delta*(log(1/cOpts.Precision))^(0.5);
        else
            lRange = cOpts.Delta*(1-cOpts.Precision)/cOpts.Precision;
        end;
    end;
    % Perform the range search
    if cOpts.UseFastNNSearch,
        [vNNInfo.NNCounts lTemp] = range_search( cX, vNNInfo.Atria, 1:lNumberOfPoints, lRange, -1);
    else
        [vNNInfo.NNCounts lTemp] = bruteforce_range_search( cX, cX, lRange );
    end;
    % Total number of entries in the W matrix
    lTotalSize = sum(vNNInfo.NNCounts);
    % We need to create the sparse W matrix all at once (Matlab is VERY slow in adding one element to a sparse matrix)
    % Allocate memory for needed indices
    lIdx = 1;lIdxsI=[];lIdxsJ=[];lIdxsI(lTotalSize) = int16(0);lIdxsJ(lTotalSize) = int16(0);lEntries = zeros(1,lTotalSize);
    % Prepare the vNNInfo structure as nn-search does
    %    if ~isfield(cOpts,'kNN'),
    for lk = 1:lNumberOfPoints
        vNNInfo.NNIdxs{lk} = lTemp{lk,1};
        vNNInfo.NNDist{lk} = lTemp{lk,2};
    end;
    %    else
    % for lk = 1:lNumberOfPoints
    %     vNNInfo.NNCounts(lk) = min([cOpts.kNN,vNNInfo.NNCounts(lk)]);           % TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %     vNNInfo.NNIdxs{lk} = lTemp{lk,1}(1:vNNInfo.NNCounts(lk));
    %     vNNInfo.NNDist{lk} = lTemp{lk,2}(1:vNNInfo.NNCounts(lk));
    % end;
    %    end;
    % Re-organize the stuff returned from the range search.
    for lk = 1:lNumberOfPoints
        % Build the sparse row and column index sets
        lIdxsI(lIdx:lIdx+vNNInfo.NNCounts(lk)-1) = lk;
        lIdxsJ(lIdx:lIdx+vNNInfo.NNCounts(lk)-1) = vNNInfo.NNIdxs{lk};
        if ~lFSpecified
            % No vector field specified
            if lWeightsFcnSpecified
                % A weight function is specified
                if isfield(cOpts,'Delta'),
                    lExtraParam.Delta = cOpts.Delta;
                end;
                if isfield(cOpts,'Epsilon'),
                    lExtraParam.Epsilon = cOpts.Epsilon;
                end;
                lExtraParam.Distances = vNNInfo.NNDist{lk};
                lExtraParam.Custom = cOpts.GetWeightsFcnParam;
                lExtraParam.CenterIdx = lk;
                % Call the weight function
                lW = feval( cOpts.GetWeightsFcn,cX,vNNInfo.NNIdxs{lk},cX(lk,:), lExtraParam );
                lEntries(lIdx:lIdx+vNNInfo.NNCounts(lk)-1) = lW;
            else
                % No weight function is specified: use standard Gaussian weight
                if cOpts.Delta == Inf
                    % Contant potential
                    lEntries(lIdx:lIdx+vNNInfo.NNCounts(lk)-1) = 1;
                elseif strcmpi(cOpts.Type,'gauss')
                    % Gaussian potential
                    lEntries(lIdx:lIdx+vNNInfo.NNCounts(lk)-1) = exp(-(vNNInfo.NNDist{lk}/cOpts.Delta).^2);
                elseif strcmpi(cOpts.Type,'harmonic')
                    % Harmonic potential
                    lEntries(lIdx:lIdx+vNNInfo.NNCounts(lk)-1) = cOpts.Delta./(cOpts.Delta+vNNinfo.NNDist{lk});
                end;
            end;
        else
            % A vector field is specified: flow along it
            lEntries(lIdx:lIdx+vNNInfo.NNCounts(lk)-1) = exp(-((vNNInfo.NNDist{lk}/cOpts.Delta).^2+abs(cF(vNNInfo.NNIdxs{lk})-cF(lk)))/lDeltaSquared);
        end;
        lIdx = lIdx + vNNInfo.NNCounts(lk);
    end;
end;

clear lTemp;

% Check if needs to be sparsified
lIdxs = [];
if isfield(cOpts,'SparsifyW'),
    fprintf('Sparsifying the graph.... \n'); 
    if (isfield(cOpts.SparsifyW,'Type')) & (isfield(cOpts.SparsifyW,'Abs')) & (isfield(cOpts.SparsifyW,'Value')),
        if cOpts.SparsifyW.Abs,
            lEntriesS = abs(lEntries);
        else
            lEntriesS = lEntries;
        end;
        if strcmpi(cOpts.SparsifyW.Type,'absolute')
            % Absolute threshold
            lIdxs = find( lEntriesS>=cOpts.SparsifyW.Value );
        elseif strcmpi(cOpts.SparsifyW.Type,'relative')
            % Relative threshold
            lIdxs = find( lEntriesS>=max(lEntriesS)*cOpts.SparsifyW.Value );
        elseif strcmpi(cOpts.SparsifyW.Type,'quantile')
            % Quantile threshold
            [lSorted lSortedIdxs] = sort(lEntriesS);
            lIdxs = lSortedIdxs(floor(cOpts.SparsifyW.Value*length(lEntriesS)):length(lEntriesS));
        end;
    end;
end;

% Fill-in the matrix.
if isempty(lIdxs),
    vGraph.W = sparse(lIdxsI,lIdxsJ,lEntries,lNumberOfPoints,lNumberOfPoints);
else
    vGraph.W = sparse(lIdxsI(lIdxs),lIdxsJ(lIdxs),lEntries(lIdxs),lNumberOfPoints,lNumberOfPoints);
end;

lIdxsI = [];lIdxsJ=[];lEntries=[];

% Symmetrize if requested
if isfield(cOpts,'NNsymm'),
    % Symmetrize it
    [lI,lJ,lV] = find(vGraph.W);
    lGraphWT = sparse(lJ,lI,lV,size(vGraph.W,1),size(vGraph.W,2),length(lI));
    lI = [];lJ=[];lV=[];
    try
        if strcmpi(cOpts.NNsymm,'ave')==1,
            vGraph.W = (lGraphWT+vGraph.W);            
            vGraph.W = vGraph.W/2;
        elseif strcmpi(cOpts.NNsymm,'leftadj')==1,
            vGraph.W = lGraphWT*vGraph.W;
        elseif strcmpi(cOpts.NNsymm,'rightadj')==1,
            vGraph.W = vGraph.W*lGraphWT;
        end;
    catch
        % Default is to average with adjoint
        vGraph.W = (lGraphWT+vGraph.W)/2;
        vGraph.W = vGraph.W/2;
    end;
    if (cOpts.NNsymm01),
        vGraph.W = (vGraph.W~=0);
    end;

end;


vGraph.D = sum(vGraph.W,2);
vOpts = cOpts;

return;
