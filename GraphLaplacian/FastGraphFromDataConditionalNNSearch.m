function [vNNIdxs,vNNDist,vNNDist1,vNNDist2]=FastGraphFromDataConditionalNNSearch(cX,cOpts);
%function [lIdxsI,lIdxsJ,lEntries]=FGFDCFNN(cX,cOpts);

%
% function [lIdxsI,lIdxsJ,lEntries]=FGFDCFNN(cX,cOpts);
%
% IN:
%   cX          : M by N matrix of M vectors in N dimensions
%   cOpts       : structure containing the following fields:
%                   NNCoordinates : List of coordinates in which to do the preliminary constraining kNNCoarse nn search
%                   kNNCoarse:  Number of nearest neighbors to compute in cX(:,cOpts.NNCoordinates)
%                   kNN: Number of nearest neighbors to compute in high-dimensions after the constrained nn search
%                   Delta: Scaling for computation of the exponential weight on the edges of the graph
%
% OUT:
%   lIdxsI,lIdxsJ,lEntries : sparse matrix of fast constrained nn-search
%
% USES:
%   nn_prepare,nn_search,
%
% SC:
%   ADS     :   9/05
%   MM      :   9/16/05. Just cleaned up, a couple of bug fixes and commented
%


lNumberOfPoints=size(cX,1);

% Do constrained NN neighbors on the specified coordinates
lXProj=cX(:,cOpts.NNCoordinates);
lXProjC =cX(:,setdiff(1:size(cX,2),cOpts.NNCoordinates));              % A bit of wasted memory here, traded for speed below in the loop.
atria=nn_prepare(lXProj);
if lNumberOfPoints<cOpts.kNNCoarse 
    cOpts.kNNCoarse = lNumberOfPoints;
end     
% Perform the nn search
[CIdxs CDist] = nn_search( lXProj, atria, lXProj, cOpts.kNNCoarse );

% Total number of entries in the W matrix
lTotalSize = cOpts.kNN*lNumberOfPoints;
% We need to create the sparse W matrix all at once (Matlab is VERY slow in adding one element to a sparse matrix)
% Allocate memory for variables needed
lIdx = 1;lIdxsI=[];lIdxsJ=[];lIdxsI(lTotalSize) = int16(0);lIdxsJ(lTotalSize) = int16(0);
vNNIdxs = zeros(lNumberOfPoints,cOpts.kNN);
vNNDist = zeros(lNumberOfPoints,cOpts.kNN);
if strcmpi(cOpts.Type,'selftuning')==1,
    vNNDist1 = zeros(lNumberOfPoints,cOpts.kNN);
    vNNDist2 = zeros(lNumberOfPoints,cOpts.kNN);   
else
    vNNDist1 = [];
    vNNDist2 = [];
end;

% Re-organize the stuff returned from the nn search.
for lk = 1:lNumberOfPoints    
    % Build the sparse row and column index sets
    [lNNIdxs lNNDist] = bruteforce_nn_search( cX(CIdxs(lk,:),:), cX(lk,:), cOpts.kNN );   
    vNNIdxs(lk,:) = CIdxs(lk,lNNIdxs);
    vNNDist(lk,:) = lNNDist;
    % No weight function is specified: use standard Gaussian weight or no weight
    if strcmpi(cOpts.Type,'selftuning'),
        vNNDist1(lk,:) = bruteforce_distances( lXProj(vNNIdxs(lk,:),:) , lXProj(lk,:),2 );
        vNNDist2(lk,:) = bruteforce_distances( lXProjC(vNNIdxs(lk,:),:), lXProjC(lk,:),2 ); 
    end;
end;    

return;