function [vNNIdxs,vNNDists] = bruteforce_nn_search( cPointSet, cPointFromSet, cNumberOfNeighbors, cP )

%
% function [vNNIdxs,vNNDists] = bruteforce_nn_search( cPointSet, cPointFromSet, cNumberOfNeighbors, cP )
%
%   Computes the cNumberOfNeighbors nearest neighbors of each point in cPointFromSet that are in cPointSet.
%   It is done brute force, and its computational complexity is O(N*M*(M')*log(M)), where M is the number of points
%   in cPointSet, M' the number of points in cPointFromSet, and N the dimensionality of the space.
%   NOTE: The algorithm does not randomly breaks ties (which is necessary if used in classification problems).
%
% IN:
%   cPointSet           : M by N matrix of M points in N dimensions
%   cPointFromSet       : M' by N matrix of M' points in N dimensions. If empty [], equals cPointSet
%   cNumberOfNeighbors  : number of nearest neighbors to find
%   [cP]                : which p-norm use as distance, p=2 by default, allowed any value in (0,+\infty)
%
% OUT:
%   vNNIdxs    : M' by cNumberOfNeighbors matrix, the i-th row being the set of indices, into cPointSet of the cNumberOfNeighbors neighbors of the i-th cPointFromSet point
%   vNNDists   : M' by cNumberOfNeighbors matrix, the i-th row being the set of distances to the cP-th power corresponding to the i-th row points in vNNidxs
%
% USES:
%   bruteforce_distances
%
% SC:
%   MM : 08/11/04
%   MM : 08/21/04 : small bug fix (cNumberOfNeighbors>number of points in cPointSet)
%
% NOTES:
%   - Writing a customized sort would reduce the log(M) term to the constant cNumberOfNeighbors.
%   - Could have plug-in function for arbitrary distances
%   
%
% (c) Copyright Plain Sight Systems Inc, 2004
%

if nargin<4,
    cP = 2;
end;

% Have to do this because bruteforce_nn_search expects full matrices: this is terrible: TODO!
cPointSet = full(cPointSet);

if isempty(cPointFromSet)
    cPointFromSet = cPointSet;
else
    cPointFromSet = full(cPointFromSet);
end;

% Get number of points and dimensions
[lNumberOfPointFrom lDimFrom] = size(cPointFromSet);
[lNumberOfPoints lDim] = size(cPointSet);

% Sanity check
if lDimFrom ~= lDim,
    warning('The dimensions of the query points does not match the dimension of the points given');
    vNNIdxs = []; vNNDists = [];
    return;
end;

% Allocate memory for local variables and results
vNNIdxs   = zeros(lNumberOfPointFrom,min([cNumberOfNeighbors,lNumberOfPoints]));
vNNDists  = zeros(lNumberOfPointFrom,min([cNumberOfNeighbors,lNumberOfPoints]));

% Compute all the distances from cPointFromSet(lk,:) to all the points in cPointSet
lDist = bruteforce_distances( cPointSet,cPointFromSet,cP );

% Loop through the query points
for lk = 1:lNumberOfPointFrom,
    % Sort the distances (actually, only need the top cNumberOfNeighbors, could be optimized if needed)
    [lTempDists, lTempNNs] = sort(lDist(lk,:));
    vNNDists(lk,:) = lTempDists(1:min([cNumberOfNeighbors,lNumberOfPoints]));
    vNNIdxs(lk,:) = lTempNNs(1:min([cNumberOfNeighbors,lNumberOfPoints]));
end;

return;