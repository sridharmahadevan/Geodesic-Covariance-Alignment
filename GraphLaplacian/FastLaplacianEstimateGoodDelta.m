function [vDelta,vDists] = FastLaplacianEstimateGoodDelta( cX )

%
% function [vDelta,vDists] = FastLaplacianEstimateGoodDelta( cX )
%
% Returns a good value for delta for Gaussian weighted exp(-(d(x,y)/delta)^2) for the set of points.
%
% IN:
%   cX          : M by N matrix for M points in N dimensions
%
% OUT:
%   vDelta      : suggested value of variance of exp. weights
%   vDists      : distances on a random subset of points of cX, used to compute vDelta
%
% USES:
%
% SC:
%   MM  :   2/13/05
%
% (c) Copyright Plain Sight Systems Inc., 2005
%

lNumberOfPoints = size(cX,1);

% Pick a subset of cX of at most 1000 random points
lN = min(floor(lNumberOfPoints/5),1000);
lRandPerm = randperm(lNumberOfPoints);
% Compute pairwise distances on the random subset
lDists = pdist(cX(lRandPerm(1:lN),:));
% Estimate a delta based on the distances
lDelta = median(lDists)/2;

% Pass return results
vDelta = lDelta;
vDists = lDists;

return;
