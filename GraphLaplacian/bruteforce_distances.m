function vDist = bruteforce_distances( cTo, cFrom, cP );

%
% function vDist =bruteforce_distances( cTo, cFrom, cP );
%   
%   Brute force computation of all cP-norm distances, to the cP-th power, from
%   any point in cFrom to any point in cTo.
%
%   NOTE: for cP<1 these are only quasi-distances
%
% IN:
%   cTo     : M by N matrix of M points in N dimensions
%   cFrom   : M' by N matrix of M' points in N dimensions
%   [cP]    : which p-norm to the p-th power to compute, p any value in (0,+\infty)
%
% OUT:
%   vDist   : M' by M cP-distances to the cP-th power, the (i,j) entry being the distance between the i-th
%             point in cFrom to the j-th point in cTo
%
% EXAMPLE:
%   distances([1.1,0;2.2,0],[2,0;3,0;4,0],1)
%
% USES:
%   distances.c
%
% SC:
%   SL  : initial version
%   MM : 8/22/04 : added p norms
%